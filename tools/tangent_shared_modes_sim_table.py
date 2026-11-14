#!/usr/bin/env python3
"""Format the experimental tangent replay with the common CARMA table writer."""
from carma_sim_table import main

if __name__ == '__main__':
    raise SystemExit(main(
        table={
            'tool': 'tools/tangent_shared_modes_sim_table.py',
            'environment': 'table*',
            'variant': '27-state reduced tangent-bias shared modes',
            'caption': 'Experimental tangent-bias replay (last \\SI{900}{s} RMS window). Single realization per record; unchanged TFG quality gates.',
            'label': 'tab:tangent_native_replay'},
        binary_default='tests/kalman_tangent_shared_modes/tangent-replay',
        data_default='tests/kalman_tangent_shared_modes',
        output_default='doc/kalman_tangent_shared_modes/tangent-replay-generated.tex-part'))
