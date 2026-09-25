"""Source text of each sea-state orchestrator, including its shared headers.

The OU-II, OU-III and TFG wrappers keep their estimator-specific layer in
their own header and share the marine front end, the adaptation mechanics,
the vibration conditioning, the magnetic startup primitives and (for the two
OU families) the proxy-startup wrapper through src/kalman_common/.  A
source-level contract about "what the OU-III implementation does" therefore
has to read the wrapper together with the shared headers it is built from.

wrapper_source(family) returns exactly that unit.  estimator_layer(family)
returns the wrapper alone, for contracts that are specifically about the
estimator-specific layer (the regularizer law, the state layout, ...).
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

WRAPPERS = {
    "ou2": "src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h",
    "ou3": "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h",
    "tfg": "src/kalman_tfg/SeaStateFusionFilter_TFG.h",
}

_COMMON = "src/kalman_common/"
_OU_SHARED = (
    _COMMON + "SeaStateOUFamilyDefaults.h",
    _COMMON + "SeaStateFusionDefaults.h",
    _COMMON + "AccelVibrationConditioning.h",
    _COMMON + "MarineWaveFrontEnd.h",
    _COMMON + "SeaStateAdaptationCommon.h",
    _COMMON + "ProxyStartupFusion.h",
    _COMMON + "MagneticStartupCommon.h",
    _COMMON + "SeaStateFusionFilterCommon.h",
)
SHARED = {
    "ou2": _OU_SHARED,
    "ou3": _OU_SHARED,
    "tfg": (
        _COMMON + "SeaStateFusionDefaults.h",
        _COMMON + "AccelVibrationConditioning.h",
        _COMMON + "MagneticStartupCommon.h",
        _COMMON + "SeaStateFusionFilterCommon.h",
    ),
}


def _text(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def estimator_layer(family: str) -> str:
    return _text(WRAPPERS[family])


def wrapper_source(family: str) -> str:
    parts = [estimator_layer(family)]
    parts.extend(_text(path) for path in SHARED[family])
    return "\n".join(parts)


def source_paths(family: str) -> tuple[str, ...]:
    return (WRAPPERS[family],) + tuple(SHARED[family])


_NUMBER = r"(-?[0-9]+(?:\.[0-9]*)?(?:[eE][-+]?[0-9]+)?)f?"
_RHS = rf"(?:{_NUMBER}\s*/\s*{_NUMBER}|{_NUMBER}|[A-Za-z_:][A-Za-z0-9_:]*)"


def _evaluate(rhs: str) -> float | str:
    import re

    quotient = re.fullmatch(rf"{_NUMBER}\s*/\s*{_NUMBER}", rhs)
    if quotient:
        return float(quotient.group(1)) / float(quotient.group(2))
    literal = re.fullmatch(_NUMBER, rhs)
    if literal:
        return float(literal.group(1))
    return rhs.split("::")[-1]


def resolved_value(source: str, name: str) -> float | None:
    """Value of the definition of `name` in source, following named defaults.

    Declarations (`float x = ...;`, `constexpr ... X = ...;`, a config field
    initializer) are preferred over later assignments such as a setter's
    `x = c;`.  The initializer may be a numeric literal, a literal quotient
    such as `1.0f / 200.0f`, or a (possibly qualified) name of another
    constant defined in the same source unit, e.g.
    `proxy_two_kp = defaults::STARTUP_PROXY_TWO_KP;`, which is followed to the
    literal it is defined as.  Returns None when no definition resolves.
    """
    import re

    def resolve(current: str, seen: frozenset[str]) -> float | None:
        if current in seen:
            return None
        pattern = rf"\b{re.escape(current)}\s*=\s*({_RHS})\s*;"
        declared = rf"(?:\bfloat|\bint|\bdouble|\bauto)\s+{pattern}"
        matches = [m.group(1) for m in re.finditer(declared, source)]
        matches += [m.group(1) for m in re.finditer(pattern, source)]
        for rhs in matches:
            value = _evaluate(rhs)
            if isinstance(value, float):
                return value
            resolved = resolve(value, seen | {current})
            if resolved is not None:
                return resolved
        return None

    return resolve(name, frozenset())
