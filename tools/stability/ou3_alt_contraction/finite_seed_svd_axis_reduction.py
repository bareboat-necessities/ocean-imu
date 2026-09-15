"""Source-specific reduction of the near-antiparallel Eigen seed obligation.

For Eigen 3.4.0's real 2x3 JacobiSVD with ComputeFullV, the requested axis
is exactly column 2 of the 3x2 QR preconditioner's Q on every returning run.
Jacobi sweeps and sorting can change only columns 0 and 1.  This is a memory
write/loop-index invariant, so it introduces no rounding error, eigenvalue gap,
or SVD-axis uniqueness assumption.  It does not prove that the run returns.

The reviewed headers are the ones shipped in Arduino Eigen 0.3.2.  Header
hashes bind the deduction to that source; target compiler/library arithmetic,
QR totality, the two-reflector output, and loop termination remain explicit.
"""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path

QUALIFICATION = 'OU3_ALT_EIGEN_340_SVD_UNUSED_COLUMN_REDUCTION_V1'
ARDUINO_ARCHIVE_SHA256 = '35648d8077ebc91471af88cb1f5397a0de4269041579b15fd2adc847f2287069'
SOURCES = {
    'Eigen/src/Geometry/Quaternion.h': '8e281fe29f4971828e790f334653fcaaafe34a7853b0eece4224a63a1ccf1fa8',
    'Eigen/src/SVD/JacobiSVD.h': '8409c413f2193d51ea9669c2f0b7d280319da95265a35c323fabc65f8a407108',
    'Eigen/src/QR/ColPivHouseholderQR.h': '191a654a935f1fdbb58468fe50eff4012ac9b2a530d57f3bdf31f295ff7c32c8',
    'Eigen/src/Householder/Householder.h': '4a6f35c282bd96fad7f18717a8b4ad66a7681f2594c4ef5594b3b313688e8b2d',
    'Eigen/src/Householder/HouseholderSequence.h': 'a89c7bbc9a8e7f0ba4d9b382c13fea480e2e411039b832b4ee584c3a624d4ce8',
    'Eigen/src/misc/RealSvd2x2.h': '5fe0795f45cd92a09b48c02259ef743bacbbe86f890cb0c0725c0c03fb56886a',
    'Eigen/src/Jacobi/Jacobi.h': '7aaafa09ea25c02f6cf78c9e08594617f31d4afced1aa0ebf0e2e74fa5c3a74f',
    'Eigen/src/Core/Redux.h': 'dcdb68827c7f9063af745c9a0e1e85cf53c0a69b67d3d4880134a2ef8ed1f3f6',
}


def audit_headers(include_root):
    """Accept the reviewed source, never an Eigen version label alone."""
    root = Path(include_root)
    for name, digest in SOURCES.items():
        path = root/name
        if not path.is_file() or sha256(path.read_bytes()).hexdigest()!=digest:
            raise ValueError('unreviewed Eigen startup source: '+name)
    return dict(SOURCES)


def build():
    # JacobiSVD::allocate uses min(rows,cols), not V.cols().  These enumerate
    # every possible index in the actual nested loops, including all possible
    # outcomes of maxCoeff during singular-value sorting.
    rows, cols = 2, 3
    diag = min(rows, cols)
    rotations = tuple((p,q) for p in range(1,diag) for q in range(p))
    swaps = tuple((i,i+pos) for i in range(diag) for pos in range(diag-i) if pos)
    written = frozenset(i for pair in (*rotations,*swaps) for i in pair)
    assert written == frozenset((0,1)) and 2 not in written
    # The real sign-correction loop changes U, never V.  Singular-value
    # scaling changes the singular-value vector, never V.  Therefore one
    # sweep and every finalization branch preserve V.col(2); induction gives
    # any finite number of sweeps.  No scalar arithmetic identity is assumed.
    return {
        'qualification': QUALIFICATION,
        'Arduino_Eigen_package_version': '0.3.2',
        'Eigen_source_version': '3.4.0',
        'Arduino_archive_sha256': ARDUINO_ARCHIVE_SHA256,
        'reviewed_header_hashes': dict(SOURCES),
        'matrix_shape': (rows,cols),
        'diag_size': diag,
        'all_Jacobi_rotation_column_pairs': rotations,
        'all_sorting_swap_column_pairs': swaps,
        'post_QR_V_write_columns': tuple(sorted(written)),
        'axis_column': 2,
        'QR_input_shape': (3,2),
        'QR_householder_count': 2,
        'returning_solver_axis_equals_QR_Q_column_2': True,
        'Jacobi_rotation_roundoff_accumulates_in_axis_column': False,
        'rank_1_and_rank_2_same_column_preservation': True,
        'source_uniform_QR_axis_totality_and_bounds_closed': False,
        'source_uniform_Jacobi_loop_termination_closed': False,
        'target_compiler_arithmetic_correspondence_closed': False,
        'near_antiparallel_JacobiSVD_solver_correspondence_qualified': False,
        'storage_search_allowed': False,
    }


def validate(report):
    return [key+' differs from reviewed Eigen column-invariance theorem'
            for key,value in build().items() if report.get(key)!=value]
