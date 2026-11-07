"""Bound dense chart geometry before export, independently of the file format.

Only the temporary rendering geometry is sampled. Source arrays, statistics,
FFTs, axis limits and the figure's original artists remain unchanged.
"""
from __future__ import annotations

import numpy as np

MAX_PLOT_POINTS = 2000


def envelope_indices(coordinates, max_points=MAX_PLOT_POINTS):
    """Ordered bucket extrema, endpoints and nonfinite breaks within a budget.

    Include extrema of every coordinate so trajectories retain turns as well as
    time-series peaks. This is a visual envelope, not signal-processing input.
    """
    data = np.column_stack(coordinates)
    if max_points < 2 * data.shape[1] + 3:
        raise ValueError("point budget too small for coordinate extrema")
    if len(data) <= max_points:
        return np.arange(len(data))
    per_bucket = 2 * data.shape[1] + 1
    buckets = (max_points - 2) // per_bucket
    selected = [0, len(data) - 1]
    for indices in np.array_split(np.arange(1, len(data) - 1), buckets):
        block = data[indices]
        valid = np.all(np.isfinite(block), axis=1)
        if not valid.all():
            selected.append(indices[np.flatnonzero(~valid)[0]])
        if valid.any():
            finite_indices = indices[valid]
            finite = block[valid]
            for column in finite.T:
                selected.extend((finite_indices[column.argmin()], finite_indices[column.argmax()]))
    return np.unique(selected)


def save_figure(figure, *args, **kwargs):
    """Export bounded line and filled-band geometry to SVG, PGF, PDF or raster.

    Restore artists even if a backend fails, so repeated exports use identical
    source data and analysis after saving still sees the full resolution.
    """
    from figure_decimation import max_line_points

    budget = max_line_points()
    restore = []
    try:
        for axes in figure.axes:
            # Several overlaid filters share the same screen width. Share the
            # budget among dense curves rather than multiplying it per filter.
            dense = sum(len(line.get_xdata()) > 256 for line in axes.lines)
            line_budget = max(min(256, budget), budget // max(1, dense))
            for line in axes.lines:
                if line.get_linestyle() in ("None", "none", "", " "):
                    continue
                is_3d = hasattr(line, "get_data_3d")
                original = line.get_data_3d() if is_3d else line.get_data()
                if len(original[0]) <= line_budget:
                    continue
                numeric = original if is_3d else (line.get_xdata(orig=False), line.get_ydata(orig=False))
                indices = envelope_indices(numeric, max_points=line_budget)
                setter = line.set_data_3d if is_3d else line.set_data
                restore.append((setter, original))
                setter(*(np.asanyarray(values)[indices] for values in numeric))
            # fill_between bands are polygon paths, not Line2D artists.
            from matplotlib.collections import PolyCollection
            for collection in axes.collections:
                if not isinstance(collection, PolyCollection):
                    continue
                for path in collection.get_paths():
                    if len(path.vertices) <= budget:
                        continue
                    vertices, codes = path.vertices, path.codes
                    indices = envelope_indices(vertices.T, max_points=budget)
                    restore.append((lambda v, c, p=path: (setattr(p, "vertices", v), setattr(p, "codes", c)),
                                    (vertices, codes)))
                    path.vertices = vertices[indices]
                    path.codes = codes[indices] if codes is not None else None
        return figure.savefig(*args, **kwargs)
    finally:
        for setter, original in reversed(restore):
            setter(*original)
