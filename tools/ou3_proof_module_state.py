#!/usr/bin/env python3
"""Scoped restoration for process-global OU-III proof backends.

Several historical proof producers install validated covariance/rotation backends
by rebinding attributes on already imported ``ou3_*`` modules.  That is harmless
for the intended standalone producer process, but it is unsafe when a newer proof
route invokes those producers inside one long-lived unittest/CI interpreter.

This helper snapshots every currently loaded ``ou3_*`` module binding, allows a
nested proof calculation to install whatever backend it needs, and restores the
original module dictionaries in ``finally``.  Newly imported ``ou3_*`` modules
are removed from ``sys.modules`` so a later consumer gets a clean import rather
than a half-installed backend.

The snapshot is deliberately about *bindings*.  Current backend installers
rebind functions/modules/scalar audit globals; they do not rely on in-place
mutation of shared containers.  If a future backend mutates shared mutable state
in place, that backend must either stop doing so or add its own explicit reset.
"""
from __future__ import annotations

from contextlib import contextmanager
import sys
from types import ModuleType
from typing import Iterator


def _tracked_modules(prefix: str) -> dict[str, ModuleType]:
    return {
        name: module
        for name, module in tuple(sys.modules.items())
        if name.startswith(prefix) and isinstance(module, ModuleType)
    }


@contextmanager
def preserve_module_bindings(prefix: str = "ou3_") -> Iterator[None]:
    """Restore all loaded proof-module attribute bindings after the body.

    This is intentionally broad: a route should not need to know which legacy
    V2/V3 installer happened to touch ``_initial_covariance``, ``SIGNED``, a
    source-structure helper, or an audit scalar.  The complete loaded-module
    namespace is the isolation boundary.
    """
    before = _tracked_modules(prefix)
    snapshots = {name: dict(module.__dict__) for name, module in before.items()}
    try:
        yield
    finally:
        # Restore pre-existing modules in place because other modules may retain
        # references to these module objects.
        for name, module in before.items():
            snap = snapshots[name]
            namespace = module.__dict__
            for key in tuple(namespace):
                if key not in snap:
                    namespace.pop(key, None)
            namespace.update(snap)

        # Modules imported only inside the scoped computation may themselves
        # contain installed global backends.  Drop them so a later import starts
        # from their source definition rather than inheriting scoped state.
        for name in tuple(sys.modules):
            if name.startswith(prefix) and name not in before:
                sys.modules.pop(name, None)
