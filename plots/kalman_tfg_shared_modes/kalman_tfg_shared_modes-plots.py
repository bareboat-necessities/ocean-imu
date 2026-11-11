#!/usr/bin/env python3
"""The common CARMA replay gallery for the 27-state TFG shared-mode estimator."""
import importlib.util
from functools import partial
from pathlib import Path

source = Path(__file__).resolve().parents[1] / "kalman_carma" / "kalman_carma-plots.py"
spec = importlib.util.spec_from_file_location("carma_replay_plots", source)
plots = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plots)

if __name__ == "__main__":
    plots.main(family="TFG_CARMA_SHARED_MODES_CSM", prefix="w3d_tfg_csm")
    # Reuse the gallery's drawing/statistics functions at publication-column
    # dimensions so the standalone method does not shrink full-page labels.
    representative = list(Path('.').glob('w3d_jonswap_H1.500_*_fusion_TFG_CARMA_SHARED_MODES_CSM.csv'))
    if len(representative) != 1:
        raise RuntimeError('Expected the unique medium JONSWAP replay for the method figures')
    df = plots.pd.read_csv(representative[0])
    scored = plots.tail(df, plots.SCORED_WINDOW_S)
    win = plots.tail(df, plots.PLOT_WINDOW_S).iloc[::plots.get_decimation_step(plots.BASE_SAMPLE_RATE_HZ)]
    full = df.iloc[::plots.get_decimation_step(plots.BASE_SAMPLE_RATE_HZ, plots.FULL_RUN_RATE_HZ)]
    plots.make_subplots = partial(plots.make_subplots, width=3.45, row_height=1.15)
    with plots.matplotlib.rc_context({'font.size': 8, 'axes.labelsize': 8,
                                     'legend.fontsize': 6, 'xtick.labelsize': 7, 'ytick.labelsize': 7}):
        base = 'w3d_tfg_csm_method_jonswap_medium'
        plots.plot_attitude(win, scored, base)
        plots.plot_vertical(win, scored, base)
        plots.plot_horizontal(win, scored, base)
        plots.plot_bias(full, base, 'acc')
        plots.plot_bias(full, base, 'gyro')
