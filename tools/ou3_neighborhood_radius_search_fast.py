#!/usr/bin/env python3
"""Prefix-optimized entry point for the OU-III sampled basin search."""
from __future__ import annotations

import ou3_neighborhood_radius_search as search
from pathlib import Path

search.DRIVER = Path(__file__).resolve().with_name("ou3_neighborhood_prefix_driver.py")

if __name__ == "__main__":
    raise SystemExit(search.main())
