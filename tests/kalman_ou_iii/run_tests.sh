#!/bin/bash -e

./kalman_ou_iii-sim
./kalman_ou_common-test
./aw_covariance_policy-test
./acc_bias_ou-test
./channel_freeze-test
./wave_period-test
./mag_hard_iron-test
./continuous_mag_hard_iron-test
./tuner_coupling-test
./tuner_schedule-test
./wave_band_sigma-test
./iss_contract-test
./rs_law-test
./startup_init-test

# PR #356 retunes the statistical horizons after removing the layered
# period/frequency and variance smoothing.  A new branch-only workflow created
# through the API did not instantiate on GitHub, whereas this test job is part
# of the repository's established pull-request workflow and is guaranteed to
# run on every synchronize event.  Keep the expensive screen branch-local and
# print the reports into the Actions log so they can be inspected before any
# defaults are changed.
if [[ "${GITHUB_ACTIONS:-}" == "true" && \
      "${GITHUB_HEAD_REF:-}" == "agent/period-log-smoothing-retune" ]]; then
    results="${GITHUB_WORKSPACE}/reports/results/period_statistics_retune"
    data="${GITHUB_WORKSPACE}/plots/kalman_ou_iii"
    mkdir -p "$results"

    echo "=== PR #356 coarse wave-period moment-horizon sweep ==="
    python3 -u "${GITHUB_WORKSPACE}/tools/ou_period_statistics_retune.py" \
      --axis moment --values 4,6,8,10 --baseline 8 --seeds 1 --jobs 4 \
      --data-dir "$data" \
      --csv "$results/coarse_moment_runs.csv" \
      --report "$results/coarse_moment.md"
    cat "$results/coarse_moment.md"

    echo "=== PR #356 coarse log-period smoothing sweep ==="
    python3 -u "${GITHUB_WORKSPACE}/tools/ou_period_statistics_retune.py" \
      --axis log --values 0,0.10,0.25,0.50,1.0 --baseline 0.50 --seeds 1 --jobs 4 \
      --data-dir "$data" \
      --csv "$results/coarse_log_runs.csv" \
      --report "$results/coarse_log.md"
    cat "$results/coarse_log.md"

    echo "=== PR #356 coarse acceleration-variance horizon sweep ==="
    python3 -u "${GITHUB_WORKSPACE}/tools/ou_period_statistics_retune.py" \
      --axis sigma_k --values 1,2,3,4,6 --baseline 2 --seeds 1 --jobs 4 \
      --data-dir "$data" \
      --csv "$results/coarse_sigma_k_runs.csv" \
      --report "$results/coarse_sigma_k.md"
    cat "$results/coarse_sigma_k.md"
fi
