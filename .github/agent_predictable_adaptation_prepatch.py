from pathlib import Path


def replace_once(path, old, new):
    p = Path(path)
    text = p.read_text()
    if old not in text:
        raise SystemExit(f'prepatch pattern not found in {path}: {old[:100]!r}')
    p.write_text(text.replace(old, new, 1))


# Normalize the compact OU-II setter so the main patch can handle both
# orchestrators with the same replacement.
replace_once(
    'src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h',
    '    void enableTuner(bool flag = true) { enable_tuner_ = flag; }\n',
    '    void enableTuner(bool flag = true) {\n        enable_tuner_ = flag;\n    }\n',
)

# Do not let the current sample's newly smoothed drift-correction noise bypass
# the one-sample staging boundary. All adaptive coefficients are committed by
# apply_pending_online_tune_() at the beginning of the next updateTime().
replace_once(
    'src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h',
    '        // Keep linear-block R_p0 / R_v0 tuning responsive in Live mode instead of\n'
    '        // waiting for slow adaptation cadence.\n'
    '        if (startup_stage_ == StartupStage::Live && enable_linear_block_) {\n'
    '            apply_R_p0_tune_();\n'
    '            apply_R_v0_tune_();\n'
    '        }\n\n',
    '        // R_p0/R_v0 are committed with the rest of the staged online\n'
    '        // schedule at the beginning of the next IMU sample.\n\n',
)
replace_once(
    'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h',
    '        // Keep linear-block R_S tuning responsive in Live mode instead of\n'
    '        // waiting for slow adaptation cadence.\n'
    '        if (startup_stage_ == StartupStage::Live && enable_linear_block_) {\n'
    '            apply_RS_tune_();\n'
    '        }\n\n',
    '        // R_S is committed with the rest of the staged online schedule at\n'
    '        // the beginning of the next IMU sample.\n\n',
)
