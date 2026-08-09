#!/usr/bin/env python3
from pathlib import Path
import re

# One-shot helper: removed after the validated source commit is published.
ROOT = Path(__file__).resolve().parents[2]

# The hard accelerometer-bias projection was removed because changing the
# nominal state without conditioning/transporting P breaks the posterior
# covariance interpretation. Remove its now-obsolete public setter too.
path = ROOT / "src/kalman_tfg/Kalman3D_Wave_TFG.h"
text = path.read_text()
old = "    void set_accel_bias_limit(T limit) { acc_bias_limit_ = std::abs(limit); }\n"
if old in text:
    text = text.replace(old, "", 1)
path.write_text(text)

# Replace the old test whose purpose was to validate periodic posterior
# covariance forcing. The correct invariant is the opposite: Sigma_aw is a
# process-model parameter and changing it must not silently rewrite P. The
# explicit reset remains legal as initialization of an inactive/newly enabled
# state block and is tested separately here.
test_path = ROOT / "tests/kalman_tfg/tfg_orchestrator-test.cpp"
t = test_path.read_text()
pattern = re.compile(
    r"// The congruent sync must replace the a_w marginal while preserving every\n"
    r"// correlation coefficient -- that is the whole point of it over a reset\.\n"
    r"void test_aw_covariance_sync_preserves_correlations\(\) \{.*?\n\}\n\n",
    re.S,
)
replacement = r'''// Sigma_aw parameterizes the OU process noise; it is not the posterior
// covariance of rho_aw. Changing the process model must therefore leave P
// untouched. An explicit reset is kept only as a deliberate prior
// initialization operation for the inactive/newly-enabled world block.
void test_aw_stationary_covariance_is_model_only() {
    using Mekf = ocean_imu::kalman::Kalman3D_Wave_TFG<float>;
    Mekf m;
    m.initialize_identity();

    auto& P = m.covariance_full();
    P.setIdentity();
    P *= 0.05f;
    for (int k = 0; k < 3; ++k) {
        P(Mekf::OFF_AW + k, Mekf::OFF_AW + k) = 0.9f;
        P(Mekf::OFF_AW + k, k) = 0.06f;
        P(k, Mekf::OFF_AW + k) = 0.06f;
    }
    const auto P_before = P;

    m.set_aw_stationary_std(Vector3f(0.4f, 0.4f, 0.25f));
    m.set_aw_time_constant(2.3f);
    check((P - P_before).cwiseAbs().maxCoeff() == 0.0f,
          "changing OU process parameters rewrote the posterior covariance");

    // Explicit prior initialization is intentionally different: it installs
    // the configured stationary a_w marginal and clears its old correlations.
    m.reset_aw_covariance_to_stationary();
    check(std::fabs(P(Mekf::OFF_AW, Mekf::OFF_AW) - 0.4f * 0.4f) < 1e-6f,
          "explicit a_w prior reset did not install Sigma_aw");
    check(std::fabs(P(Mekf::OFF_AW + 2, Mekf::OFF_AW + 2) - 0.25f * 0.25f) < 1e-6f,
          "explicit a_w prior reset did not install anisotropic Sigma_aw");
    check(std::fabs(P(Mekf::OFF_AW, 0)) == 0.0f,
          "explicit a_w prior reset did not clear stale cross-covariance");
}

'''
t2, n = pattern.subn(replacement, t, count=1)
if n != 1:
    raise RuntimeError(f"expected one obsolete covariance-sync test, found {n}")
t2 = t2.replace(
    "    test_aw_covariance_sync_preserves_correlations();\n",
    "    test_aw_stationary_covariance_is_model_only();\n",
    1,
)
test_path.write_text(t2)

print("Removed obsolete bias-limit API and replaced covariance-sync test")
