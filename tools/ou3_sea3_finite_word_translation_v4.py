#!/usr/bin/env python3
"""Diagnostic wrapper: use eleven conservative commit segments.

V3 deliberately admits a 20-sample commit gap to absorb binary32 clock-edge
uncertainty.  Eleven such segments make the recurrent word strictly longer
than the 1 s PE recurrence even with the outward dt lower endpoint.  Once the
finite-word construction closes this wrapper is folded back into the canonical
producer; it is not a second theorem route.
"""
import ou3_sea3_finite_word_translation_v3 as BASE

BASE.WORD_SEGMENTS = 11

if __name__ == "__main__":
    raise SystemExit(BASE.main())
