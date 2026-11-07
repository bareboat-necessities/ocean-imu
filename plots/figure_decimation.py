#!/usr/bin/env python3
"""Uniform decimation of long plotted lines in every saved figure.

Committed SVG figures and the PGF figures the articles input are vector
files: every vertex of every line is written out.  A 20-minute 200 Hz trace
is 240 000 vertices per line, megabytes of SVG and pages that crawl in a PDF
viewer, while a page-wide axis resolves a few thousand points at most.

``install()`` wraps ``matplotlib.figure.Figure.savefig`` so that, before any
figure is saved in any format, every 2-D line drawn with a line style and
more than ``max_line_points()`` vertices is reduced to that budget:

* lines with monotonic x (time series, spectra) keep, in each of
  ``budget / 2`` equal index buckets, the samples at the bucket's minimum and
  maximum of y, in time order -- the envelope and every peak survive.  The
  budget is spent on the part inside the axis x-limits (a zoomed window of a
  long record stays sharp); outside it only the endpoints, the extremes and
  the samples next to the window edges are kept, so data limits and
  autoscaling are unchanged;
* other lines (trajectories) keep every k-th vertex, endpoints included.

Marker-only lines (scatter-like data) and 3-D lines are not touched; bands
drawn with ``fill_between`` are decimated by the script with ``stride_xy()``.
The budget is ``PLOT_MAX_POINTS`` (environment) or 2000 per line.  Decimation
happens at save time in the generating script, so a re-run study produces
decimated figures by construction; the data files are never changed.
"""

from __future__ import annotations

import os

import numpy as np

DEFAULT_MAX_LINE_POINTS = 2000
_INSTALLED_ATTR = "_ocean_imu_line_decimation"


def max_line_points() -> int:
    raw = os.getenv("PLOT_MAX_POINTS", "").strip()
    try:
        value = int(raw) if raw else DEFAULT_MAX_LINE_POINTS
    except ValueError:
        value = DEFAULT_MAX_LINE_POINTS
    return max(16, value)


def minmax_indices(y: np.ndarray, max_points: int) -> np.ndarray:
    """Indices of each bucket's minimum and maximum of ``y`` (NaN-aware), in order."""
    n = len(y)
    buckets = max(1, max_points // 2)
    edges = np.linspace(0, n, buckets + 1).astype(int)
    keep = []
    for a, b in zip(edges[:-1], edges[1:]):
        if b <= a:
            continue
        seg = y[a:b]
        finite = np.isfinite(seg)
        if not finite.any():
            keep.append(a)            # keep the gap visible
            continue
        lo = int(np.nanargmin(np.where(finite, seg, np.nan)))
        hi = int(np.nanargmax(np.where(finite, seg, np.nan)))
        keep.extend(sorted({a + lo, a + hi}))
    keep = sorted(set(keep) | {0, n - 1})
    return np.asarray(keep, dtype=int)


def stride_indices(n: int, max_points: int) -> np.ndarray:
    step = int(np.ceil(n / max_points))
    idx = np.arange(0, n, step)
    if idx[-1] != n - 1:
        idx = np.append(idx, n - 1)
    return idx


def stride_xy(x, *ys, max_points: int | None = None):
    """Uniform-stride decimation of aligned arrays to the line budget, for
    data drawn as polygons (``fill_between`` bands), which ``install()`` does
    not touch."""
    budget = max_line_points() if max_points is None else max_points
    x = np.asarray(x)
    if len(x) <= budget:
        return (x, *(np.asarray(y) for y in ys))
    idx = stride_indices(len(x), budget)
    return (x[idx], *(np.asarray(y)[idx] for y in ys))


def _as_float(a) -> np.ndarray | None:
    try:
        arr = np.ma.filled(np.ma.asarray(a, dtype=float), np.nan)
    except (TypeError, ValueError):
        return None
    return np.asarray(arr, dtype=float).reshape(-1)


def windowed_minmax_indices(x: np.ndarray, y: np.ndarray, max_points: int,
                            view: tuple[float, float] | None) -> np.ndarray:
    """``minmax_indices`` with the budget spent inside ``view`` (x-limits).

    ``x`` is monotonic.  Outside the view only each side's endpoints and
    extremes survive, plus the sample next to each window edge so the segment
    crossing the edge is exact; the data limits are those of the full line."""
    n = len(y)
    inside = np.zeros(n, dtype=bool)
    if view is not None:
        lo, hi = min(view), max(view)
        with np.errstate(invalid="ignore"):
            inside = (x >= lo) & (x <= hi)
    if view is None or inside.all():
        return minmax_indices(y, max_points)
    keep = [0, n - 1]
    pos = np.flatnonzero(inside)
    if len(pos):
        a, b = max(int(pos[0]) - 1, 0), min(int(pos[-1]) + 1, n - 1)
        keep.extend(a + minmax_indices(y[a:b + 1], max_points))
    else:
        a = b = n                     # no sample in view: one outside part
    for s0, s1 in ((0, a), (b + 1, n)):
        if s1 - s0 > 0:
            keep.extend(s0 + minmax_indices(y[s0:s1], 4))
    return np.asarray(sorted(set(keep)), dtype=int)


def _x_view(line) -> tuple[float, float] | None:
    """The axis x-limits when the line's x is in data coordinates."""
    ax = line.axes
    if ax is None:
        return None
    try:
        x_is_data = line.get_transform().contains_branch_seperately(ax.transData)[0]
    except (AttributeError, ValueError):
        return None
    return tuple(ax.get_xlim()) if x_is_data else None


def decimate_line(line, max_points: int) -> bool:
    """Decimate one ``Line2D`` in place; return whether it changed."""
    if hasattr(line, "get_data_3d"):
        return False
    if line.get_linestyle() in ("None", "none", "", " "):
        return False
    x = _as_float(line.get_xdata(orig=True))
    y = _as_float(line.get_ydata(orig=True))
    if x is None or y is None or len(x) != len(y) or len(x) <= max_points:
        return False
    dx = np.diff(x[np.isfinite(x)])
    if len(dx) and (np.all(dx >= 0) or np.all(dx <= 0)):
        idx = windowed_minmax_indices(x, y, max_points, _x_view(line))
    else:
        idx = stride_indices(len(x), max_points)
    if len(idx) >= len(x):
        return False
    line.set_data(x[idx], y[idx])
    return True


def decimate_figure(figure, max_points: int | None = None) -> int:
    """Decimate every eligible line of ``figure``; return how many changed."""
    budget = max_line_points() if max_points is None else max_points
    changed = 0
    for ax in figure.get_axes():
        for line in ax.get_lines():
            changed += decimate_line(line, budget)
    return changed


def install(matplotlib_module=None) -> None:
    """Decimate every figure before ``savefig`` (idempotent).  Without
    Matplotlib there is nothing to save, so tools that plot lazily still import."""
    try:
        import matplotlib.figure as mfig
    except ImportError:
        return

    if getattr(mfig.Figure.savefig, _INSTALLED_ATTR, False):
        return
    original = mfig.Figure.savefig

    def savefig(self, *args, **kwargs):
        decimate_figure(self)
        return original(self, *args, **kwargs)

    setattr(savefig, _INSTALLED_ATTR, True)
    savefig.__doc__ = original.__doc__
    mfig.Figure.savefig = savefig
