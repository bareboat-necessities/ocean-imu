"""High-precision feasibility diagnostic for one magnetically informed H18
service superword. It is not a certificate and cannot become one.

The controlling obligation is a finite-error storage inequality

    V_(j+1) <= rho V_j + c_d ||d||^2,   rho < 1,

with coercivity and every-prefix retention on the same physical history. The
architecture requires a non-promoting feasibility check before any rigorous
enclosure of such a construction: if the worst admissible ratio of the
candidate storage is already at or above one on an actually reached execution,
the formulation has to change instead of unrelated bounds being sharpened.

This module runs exactly that check for the candidate storage

    V(e) = e^T P^(-1) e

built from the shipping covariance itself, on the export produced by
`h18_superword_export.cpp`. Central differences about a reached nonzero error measure a local incremental
surrogate, not the complete finite-error map or its affine supply. The reference
forcing and nonlinear remainder still require separate same-history bounds. Source
uniformity, the supply constant, capture, recurrence, release retention and
arithmetic totality are not addressed here at all, and no result of this module
sets any proof obligation.

All evaluation is done in `decimal` at a precision far beyond the
single-precision estimator that produced the map, so the reported numbers are
limited by the exported shipping map rather than by this evaluation.
"""
from __future__ import annotations

import argparse
import json
import math
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Sequence

ERROR_DIMENSION = 21
ATTITUDE_OFFSET = 0
GYRO_BIAS_OFFSET = 3
ACC_BIAS_OFFSET = 18

HERE = Path(__file__).resolve().parent
CONSTANTS = HERE / "constants.json"

try:  # imported as part of the theorem package
    from .finite_error import HeldBiasSuperwordBlock, held_bias_non_contraction
except ImportError:  # executed directly from the command line
    from finite_error import HeldBiasSuperwordBlock, held_bias_non_contraction

Matrix = list[list[Decimal]]
Vector = list[Decimal]


class SuperwordExportError(ValueError):
    """The export does not describe a usable magnetically informed H18 superword."""


# --------------------------------------------------------------------------
# dense linear algebra at arbitrary precision, without third-party packages
# --------------------------------------------------------------------------

def zeros(rows: int, cols: int) -> Matrix:
    return [[Decimal(0) for _ in range(cols)] for _ in range(rows)]


def matvec(a: Matrix, x: Vector) -> Vector:
    return [sum((row[j] * x[j] for j in range(len(x))), Decimal(0)) for row in a]


def transpose_matvec(a: Matrix, x: Vector) -> Vector:
    cols = len(a[0])
    out = [Decimal(0)] * cols
    for row, xi in zip(a, x):
        if xi == 0:
            continue
        for j in range(cols):
            out[j] += row[j] * xi
    return out


def dot(x: Vector, y: Vector) -> Decimal:
    return sum((xi * yi for xi, yi in zip(x, y)), Decimal(0))


def exact(value: float) -> Decimal:
    """Exact decimal image of a float; the export carries round-trip literals."""
    return Decimal(repr(float(value)))


def symmetric_from_upper(values: Sequence[float], n: int) -> Matrix:
    expected = n * (n + 1) // 2
    if len(values) != expected:
        raise SuperwordExportError(f"expected {expected} upper-triangular entries")
    out = zeros(n, n)
    at = 0
    for i in range(n):
        for j in range(i, n):
            value = exact(values[at])
            out[i][j] = value
            out[j][i] = value
            at += 1
    return out


def square_from_rows(values: Sequence[float], n: int) -> Matrix:
    if len(values) != n * n:
        raise SuperwordExportError(f"expected {n * n} matrix entries")
    return [[exact(values[i * n + j]) for j in range(n)] for i in range(n)]


def ldl_factor(a: Matrix) -> tuple[Matrix, Vector]:
    """Symmetric LDL^T factorization; fails unless the matrix is positive definite."""
    n = len(a)
    lower = zeros(n, n)
    diag = [Decimal(0)] * n
    for j in range(n):
        pivot = a[j][j]
        for k in range(j):
            pivot -= lower[j][k] * lower[j][k] * diag[k]
        if pivot <= 0:
            raise SuperwordExportError(
                f"storage metric is not positive definite at index {j}")
        diag[j] = pivot
        lower[j][j] = Decimal(1)
        for i in range(j + 1, n):
            entry = a[i][j]
            for k in range(j):
                entry -= lower[i][k] * lower[j][k] * diag[k]
            lower[i][j] = entry / pivot
    return lower, diag


def ldl_solve(factor: tuple[Matrix, Vector], b: Vector) -> Vector:
    lower, diag = factor
    n = len(b)
    y = list(b)
    for i in range(n):
        for k in range(i):
            y[i] -= lower[i][k] * y[k]
    for i in range(n):
        y[i] /= diag[i]
    for i in range(n - 1, -1, -1):
        for k in range(i + 1, n):
            y[i] -= lower[k][i] * y[k]
    return y


def inverse(a: Matrix) -> Matrix:
    """Gauss-Jordan inverse with partial pivoting, for the general root map."""
    n = len(a)
    work = [row[:] + [Decimal(1) if i == j else Decimal(0) for j in range(n)]
            for i, row in enumerate(a)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda row: abs(work[row][col]))
        if work[pivot][col] == 0:
            raise SuperwordExportError(f"singular root response at column {col}")
        work[col], work[pivot] = work[pivot], work[col]
        scale = work[col][col]
        work[col] = [value / scale for value in work[col]]
        for row in range(n):
            if row == col:
                continue
            factor = work[row][col]
            if factor == 0:
                continue
            work[row] = [value - factor * other
                         for value, other in zip(work[row], work[col])]
    return [row[n:] for row in work]


def storage(metric: tuple[Matrix, Vector], e: Vector) -> Decimal:
    """V(e) = e^T P^(-1) e for the shipping covariance P behind `metric`."""
    return dot(e, ldl_solve(metric, e))


# --------------------------------------------------------------------------
# storage ratio over one superword
# --------------------------------------------------------------------------

class SuperwordMap:
    """Psi_k = E_k M, with M the inverse of the measured root response.

    The map is kept as its two factors so that both Psi x and Psi^T y stay
    matrix-vector products; no dense product is rebuilt per iterate.
    """

    def __init__(self, difference: Matrix, root_inverse: Matrix):
        self.difference = difference
        self.root_inverse = root_inverse

    def __call__(self, x: Vector) -> Vector:
        return matvec(self.difference, matvec(self.root_inverse, x))

    def transpose(self, y: Vector) -> Vector:
        return transpose_matvec(self.root_inverse, transpose_matvec(self.difference, y))

    def dense(self, n: int) -> Matrix:
        columns = [self([Decimal(1) if i == j else Decimal(0) for i in range(n)])
                   for j in range(n)]
        return [[columns[j][i] for j in range(n)] for i in range(n)]


def _rayleigh_iteration(superword: "SuperwordMap", root_covariance: Matrix,
                        root_metric: tuple[Matrix, Vector],
                        end_metric: tuple[Matrix, Vector], seed: Vector,
                        iterations: int, tolerance: Decimal
                        ) -> tuple[Decimal, Vector, Decimal, int, bool]:
    """One normalized power iteration of P_root A from `seed`.

    A = Psi^T P_end^(-1) Psi and B = P_root^(-1). The operator B^(-1)A = P_root A
    is self-adjoint in the B inner product, so the iteration converges to the
    eigenvalue whose eigenspace the seed has a component in -- the dominant one
    unless the seed happens to be B-orthogonal to it.
    """
    x = list(seed)
    ratio = Decimal(0)
    residual = Decimal(0)
    for step in range(iterations):
        a_x = superword.transpose(ldl_solve(end_metric, superword(x)))
        denominator = storage(root_metric, x)
        if denominator <= 0:
            raise SuperwordExportError("degenerate iterate in the storage metric")
        updated = dot(x, a_x) / denominator
        residual = abs(updated - ratio)
        ratio = updated
        nxt = matvec(root_covariance, a_x)
        scale = storage(root_metric, nxt)
        if scale <= 0:
            # The map sends this seed to zero: its ratio is exactly zero and
            # there is nothing left to iterate on. Other seeds carry the rest.
            return Decimal(0), x, Decimal(0), step + 1, True
        x = [value / scale.sqrt() for value in nxt]
        if step > 2 and residual < tolerance * (Decimal(1) + abs(ratio)):
            return ratio, x, residual, step + 1, True
    return ratio, x, residual, iterations, False


def worst_admissible_ratio(superword: "SuperwordMap", root_covariance: Matrix,
                           root_metric: tuple[Matrix, Vector],
                           end_metric: tuple[Matrix, Vector], n: int,
                           iterations: int = 400,
                           tolerance: Decimal = Decimal("1e-30")
                           ) -> tuple[Decimal, Vector, Decimal, int]:
    """Largest V_end(Psi e)/V_root(e) over nonzero directions, with its argument.

    The generalized problem is A x = lambda B x with A = Psi^T P_end^(-1) Psi
    and B = P_root^(-1), and the ratio is the largest eigenvalue of B^(-1)A.

    A single power iteration is not enough to claim that largest eigenvalue: a
    seed B-orthogonal to the dominant eigenspace converges, cleanly and with a
    small Rayleigh-quotient increment, to a subdominant eigenvalue instead. The
    shipping map has exactly invariant subspaces -- the held accelerometer bias
    is one -- so that case is real here rather than hypothetical, and reporting
    a subdominant value would understate the ratio and could assert contraction
    that does not hold.

    The iteration therefore runs from every coordinate direction. Those seeds
    span the space, so at least one has a component in the dominant eigenspace
    and the maximum over the runs is the largest eigenvalue. A seed that has
    not converged within the budget fails the diagnostic closed rather than
    contributing an unestablished number.
    """
    best_ratio = Decimal(0)
    best_direction: Vector | None = None
    worst_residual = Decimal(0)
    total_steps = 0
    for index in range(n):
        seed = [Decimal(1) if i == index else Decimal(0) for i in range(n)]
        scale = storage(root_metric, seed)
        if scale <= 0:
            raise SuperwordExportError("degenerate storage metric at the root")
        seed = [value / scale.sqrt() for value in seed]
        ratio, direction, residual, steps, converged = _rayleigh_iteration(
            superword, root_covariance, root_metric, end_metric, seed,
            iterations, tolerance)
        if not converged:
            raise SuperwordExportError(
                f"the storage ratio did not converge from direction {index} "
                f"within {iterations} iterations; no ratio is established")
        total_steps += steps
        worst_residual = max(worst_residual, residual)
        if best_direction is None or ratio > best_ratio:
            best_ratio, best_direction = ratio, direction
    if best_direction is None:
        raise SuperwordExportError("no admissible direction was evaluated")
    return best_ratio, best_direction, worst_residual, total_steps


def eigenvalue_residual(superword: "SuperwordMap", root_covariance: Matrix,
                        root_metric: tuple[Matrix, Vector],
                        end_metric: tuple[Matrix, Vector],
                        ratio: Decimal, direction: Vector) -> Decimal:
    """Relative B-norm residual of A x = lambda B x at the reported pair.

    This is an a posteriori check on the returned eigenpair, reported for the
    reader. It is not a bound on anything the diagnostic is measuring.
    """
    a_x = superword.transpose(ldl_solve(end_metric, superword(direction)))
    gap = [value - ratio * component
           for value, component in zip(matvec(root_covariance, a_x), direction)]
    scale = storage(root_metric, direction)
    if scale <= 0:
        raise SuperwordExportError("degenerate limiting direction")
    return (storage(root_metric, gap) / scale).sqrt()


# --------------------------------------------------------------------------
# magnetic information actually applied inside the superword
# --------------------------------------------------------------------------

def applied_magnetic_information(export: dict, maps: Sequence[SuperwordMap],
                                 scales: dict) -> tuple[Decimal, list[dict]]:
    """Smallest eigenvalue of sum G_k^T G_k in normalized root coordinates.

    G_k = W_k H_(m,k) Phi_(k<-s) E_hb with W_k^T W_k = (S_(m,k)^act)^(-1), so
    G_k^T G_k = M_k^T (S_(m,k)^act)^(-1) M_k for M_k = H_(m,k) Phi_(k<-s) E_hb.
    Only corrections the shipping update actually applied appear in the export,
    and nothing here promotes one window to recurring service.
    """
    rotation = square_from_rows(export["root_rotation_world_to_body"], 3)
    axial = [rotation[i][2] for i in range(3)]   # world down axis in body coordinates
    heading_scale = exact(scales["heading_coordinate_scale_rad"])
    bias_scale = exact(scales["axial_gyro_bias_coordinate_scale_rad_s"])

    injection: list[Vector] = []
    for offset, scale in ((ATTITUDE_OFFSET, heading_scale), (GYRO_BIAS_OFFSET, bias_scale)):
        column = [Decimal(0)] * ERROR_DIMENSION
        for i in range(3):
            column[offset + i] = axial[i] * scale
        injection.append(column)

    total = [[Decimal(0), Decimal(0)], [Decimal(0), Decimal(0)]]
    contributions = []
    for event in export["applied_magnetic_corrections"]:
        prefix = int(event["prefix"])
        sensitivity = [exact(v) for v in event["sensitivity_axis"]]
        innovation = ldl_factor(square_from_rows(event["innovation_covariance"], 3))
        # The pre-correction response and pre-correction sensitivity must be
        # paired with the innovation covariance this correction actually used.
        preceding = SuperwordMap(
            square_from_rows(event["pre_error_difference"], ERROR_DIMENSION),
            maps[0].root_inverse)
        transported = [preceding(column) for column in injection]
        # H_(m,k) = -[zhat]x on the attitude block, zero elsewhere.
        sensed = []
        for column in transported:
            attitude = column[ATTITUDE_OFFSET:ATTITUDE_OFFSET + 3]
            sensed.append([
                -(sensitivity[1] * attitude[2] - sensitivity[2] * attitude[1]),
                -(sensitivity[2] * attitude[0] - sensitivity[0] * attitude[2]),
                -(sensitivity[0] * attitude[1] - sensitivity[1] * attitude[0]),
            ])
        block = [[dot(sensed[i], ldl_solve(innovation, sensed[j])) for j in range(2)]
                 for i in range(2)]
        for i in range(2):
            for j in range(2):
                total[i][j] += block[i][j]
        contributions.append({
            "prefix": prefix,
            "time_s": float(event["time_s"]),
            "information_trace": float(block[0][0] + block[1][1]),
        })

    trace = total[0][0] + total[1][1]
    determinant = total[0][0] * total[1][1] - total[0][1] * total[1][0]
    discriminant = trace * trace - Decimal(4) * determinant
    if discriminant < 0:
        discriminant = Decimal(0)
    smallest = (trace - discriminant.sqrt()) / Decimal(2)
    return smallest, contributions


# --------------------------------------------------------------------------
# validation and top-level evaluation
# --------------------------------------------------------------------------

REQUIRED_KEYS = (
    "error_dimension", "superword_samples", "reference_error", "covariance_upper",
    "error_difference", "root_difference_coarse", "endpoint_difference_coarse",
    "applied_magnetic_corrections", "root_rotation_world_to_body",
    "acc_bias_held_through_superword", "attitude_injection_finite", "inherited_state",
)


def validate(export: dict) -> None:
    for key in REQUIRED_KEYS:
        if key not in export:
            raise SuperwordExportError(f"missing export field: {key}")
    if export["error_dimension"] != ERROR_DIMENSION:
        raise SuperwordExportError("unexpected shipping error dimension")
    samples = int(export["superword_samples"])
    if samples != export["superword_samples"]:
        raise SuperwordExportError("superword_samples must be an integer")
    if samples < 1:
        raise SuperwordExportError("an empty superword proves nothing")
    for key in ("reference_error", "covariance_upper", "error_difference"):
        if len(export[key]) != samples + 1:
            raise SuperwordExportError(
                f"{key} must carry the root and every prefix of the superword")
    if not export["applied_magnetic_corrections"]:
        raise SuperwordExportError(
            "a magnetically informed superword needs an actually applied correction")
    previous_prefix = 0
    for event in export["applied_magnetic_corrections"]:
        prefix = int(event["prefix"])
        if prefix != event["prefix"] or not previous_prefix < prefix <= samples:
            raise SuperwordExportError(
                "applied corrections need distinct ordered integer prefixes in the superword")
        previous_prefix = prefix
        if event.get("phase") != "pre_correction":
            raise SuperwordExportError("magnetic sensitivity/transport must precede correction")
        for field, size in (("pre_error_difference", ERROR_DIMENSION**2),
                            ("sensitivity_axis", 3), ("innovation_covariance", 9)):
            values = event.get(field, [])
            if len(values) != size or any(not math.isfinite(float(v)) for v in values):
                raise SuperwordExportError(f"invalid magnetic {field}")
    for flag in ("acc_bias_held_through_superword", "attitude_injection_finite",
                 "inherited_state"):
        if export[flag] is not True:
            raise SuperwordExportError(f"the export does not satisfy {flag}")
    for key, size in (("reference_error", ERROR_DIMENSION),
                      ("covariance_upper", ERROR_DIMENSION*(ERROR_DIMENSION+1)//2),
                      ("error_difference", ERROR_DIMENSION**2)):
        for row in export[key]:
            if len(row) != size:
                raise SuperwordExportError(f"invalid row size in {key}")
            for value in row:
                if not math.isfinite(float(value)):
                    raise SuperwordExportError(f"nonfinite entry in {key}")
    for key, size in (("root_difference_coarse", ERROR_DIMENSION**2),
                      ("endpoint_difference_coarse", ERROR_DIMENSION**2),
                      ("root_rotation_world_to_body", 9)):
        if len(export[key]) != size or any(not math.isfinite(float(v)) for v in export[key]):
            raise SuperwordExportError(f"invalid entries in {key}")


def _relative_map_discrepancy(fine: Matrix, coarse: Matrix) -> Decimal:
    """Largest entrywise gap between the two difference scales, relative to scale."""
    worst = Decimal(0)
    magnitude = Decimal(0)
    for row_fine, row_coarse in zip(fine, coarse):
        for a, b in zip(row_fine, row_coarse):
            worst = max(worst, abs(a - b))
            magnitude = max(magnitude, abs(a), abs(b))
    if magnitude == 0:
        return Decimal(0)
    return worst / magnitude


def evaluate(export: dict, *, precision: int = 60, prefix_stride: int = 1) -> dict:
    """Measure the candidate storage over one superword. Nothing is promoted."""
    validate(export)
    if prefix_stride < 1:
        raise SuperwordExportError("prefix stride must be positive")
    constants = json.loads(CONSTANTS.read_text(encoding="utf-8"))
    service = constants["magnetic_service"]
    history_id = (f"h18-root-{float(export.get('root_time_s', 0.0)):.3f}s"
                  f"-window-{int(export['superword_samples'])}")

    with localcontext() as context:
        context.prec = precision
        n = ERROR_DIMENSION
        prefixes = len(export["error_difference"])
        differences = [square_from_rows(row, n) for row in export["error_difference"]]
        covariances = [symmetric_from_upper(row, n) for row in export["covariance_upper"]]
        root_inverse = inverse(differences[0])
        maps = [SuperwordMap(difference, root_inverse) for difference in differences]

        root_metric = ldl_factor(covariances[0])
        end_metric = ldl_factor(covariances[-1])
        ratio, direction, residual, iterations = worst_admissible_ratio(
            maps[-1], covariances[0], root_metric, end_metric, n)
        eigenpair_residual = eigenvalue_residual(
            maps[-1], covariances[0], root_metric, end_metric, ratio, direction)

        # The same map measured in a frozen metric separates amplification of
        # the shipping map from drift of the covariance that weights it.
        frozen_ratio, _, _, _ = worst_admissible_ratio(
            maps[-1], covariances[0], root_metric, root_metric, n)

        root_storage = storage(root_metric, direction)
        profile = []
        retention = Decimal(0)
        checkpoints = sorted(set(list(range(0, prefixes, prefix_stride)) + [prefixes - 1]))
        for prefix in checkpoints:
            metric = root_metric if prefix == 0 else ldl_factor(covariances[prefix])
            value = storage(metric, maps[prefix](direction)) / root_storage
            retention = max(retention, value)
            profile.append((prefix, value))

        reference = []
        for prefix in (0, prefixes - 1):
            metric = root_metric if prefix == 0 else end_metric
            error = [exact(v) for v in export["reference_error"][prefix]]
            reference.append(float(storage(metric, error)))

        information, events = applied_magnetic_information(export, maps, service)

        coarse_root = square_from_rows(export["root_difference_coarse"], n)
        coarse_end = square_from_rows(export["endpoint_difference_coarse"], n)
        coarse_map = SuperwordMap(coarse_end, inverse(coarse_root))
        coarse_ratio, _, _, _ = worst_admissible_ratio(
            coarse_map, covariances[0], root_metric, end_metric, n)
        endpoint = maps[-1].dense(n)
        conditioning = _relative_map_discrepancy(endpoint, coarse_map.dense(n))

        # The held accelerometer-bias block of the shipping map, measured rather
        # than assumed: in H18 the estimator applies no bias dynamics and freezes
        # the bias rows of the gain, so this block is expected to be the identity.
        held_block_deviation = Decimal(0)
        for i in range(ACC_BIAS_OFFSET, n):
            for j in range(ACC_BIAS_OFFSET, n):
                unit = Decimal(1) if i == j else Decimal(0)
                held_block_deviation = max(held_block_deviation,
                                           abs(endpoint[i][j] - unit))

        # The two remaining literal inputs of the held-bias obstruction: the
        # bias cross-covariances the hold zeroed, and the frozen bias block.
        cross_covariance = Decimal(0)
        for covariance in (covariances[0], covariances[-1]):
            for i in range(ACC_BIAS_OFFSET, n):
                for j in range(ACC_BIAS_OFFSET):
                    cross_covariance = max(cross_covariance, abs(covariance[i][j]))
        block_change = Decimal(0)
        for i in range(ACC_BIAS_OFFSET, n):
            for j in range(ACC_BIAS_OFFSET, n):
                block_change = max(block_change,
                                   abs(covariances[-1][i][j] - covariances[0][i][j]))
        obstruction = held_bias_non_contraction(HeldBiasSuperwordBlock(
            history_id=history_id,
            map_block=tuple(tuple(float(endpoint[ACC_BIAS_OFFSET + i][ACC_BIAS_OFFSET + j])
                                  for j in range(3)) for i in range(3)),
            cross_covariance_max=float(cross_covariance),
            covariance_block_change=float(block_change)))

        floor = exact(service["mu_M"])
        report = {
            "qualification": "OU3_H18_SERVICE_SUPERWORD_FEASIBILITY_V2",
            "history_id": history_id,
            "role": ("non-promoting feasibility diagnostic of one candidate storage on "
                     "one reached execution; not a certificate and not source uniform"),
            "storage": "V(e) = e^T P^-1 e with the shipping covariance P",
            "evaluation_precision_digits": precision,
            "superword_samples": int(export["superword_samples"]),
            "prefix_checkpoints": len(profile),
            "worst_admissible_ratio": float(ratio),
            "worst_admissible_ratio_coarse_scale": float(coarse_ratio),
            "worst_admissible_ratio_frozen_metric": float(frozen_ratio),
            "held_bias_block_deviation_from_identity": float(held_block_deviation),
            "held_bias_obstruction": obstruction,
            "map_scale_discrepancy": float(conditioning),
            "power_iteration_residual": float(residual),
            "power_iteration_steps": iterations,
            "power_iteration_seeds": ERROR_DIMENSION,
            "eigenpair_relative_residual": float(eigenpair_residual),
            "strict_contraction_observed": bool(ratio < 1),
            "endpoint_direction_prefix_ratio_max": float(retention),
            "endpoint_direction_prefix_nonexpansive": bool(retention <= 1),
            "limiting_direction": [float(v) for v in direction],
            "limiting_direction_blocks": _direction_blocks(direction),
            "endpoint_direction_prefix_profile": [{"prefix": prefix, "storage_ratio": float(value)}
                               for prefix, value in profile],
            "reference_storage_root": reference[0],
            "reference_storage_end": reference[1],
            "applied_magnetic_corrections": len(events),
            "magnetic_information_min_eigenvalue": float(information),
            "magnetic_information_floor": float(floor),
            "magnetic_information_meets_floor": bool(information >= floor),
            "magnetic_event_contributions": events,
            "local_incremental_map_only": True,
            "finite_error_remainder_bounded": False,
            "all_direction_prefix_retention_evaluated": False,
            "supply_constant_evaluated": False,
            "source_uniform": False,
            "certificate_complete": False,
            "obligation_discharged": False,
        }
    return report


BLOCK_NAMES = ("attitude", "gyro_bias", "velocity", "position", "integral",
               "wave_acceleration", "accelerometer_bias")


def _direction_blocks(direction: Vector) -> dict:
    """Block energy shares of the limiting direction, for the failure report."""
    total = sum((value * value for value in direction), Decimal(0))
    out = {}
    for index, name in enumerate(BLOCK_NAMES):
        block = direction[3 * index:3 * index + 3]
        share = sum((value * value for value in block), Decimal(0))
        out[name] = float(share / total) if total > 0 else 0.0
    return out


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export", type=Path, required=True,
                        help="JSON produced by h18_superword_export")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--precision", type=int, default=60)
    parser.add_argument("--prefix-stride", type=int, default=1)
    args = parser.parse_args(argv)

    export = json.loads(args.export.read_text(encoding="utf-8"))
    report = evaluate(export, precision=args.precision, prefix_stride=args.prefix_stride)
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    summary = {key: report[key] for key in (
        "worst_admissible_ratio", "strict_contraction_observed",
        "endpoint_direction_prefix_ratio_max", "magnetic_information_min_eigenvalue",
        "magnetic_information_meets_floor", "certificate_complete")}
    print(json.dumps(summary, sort_keys=True))
    # The diagnostic reports; it never gates on its own numbers.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
