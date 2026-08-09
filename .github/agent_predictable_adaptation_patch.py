from pathlib import Path


def replace_once(path, old, new):
    p = Path(path)
    text = p.read_text()
    if old not in text:
        raise SystemExit(f"pattern not found in {path}: {old[:100]!r}")
    if text.count(old) != 1:
        raise SystemExit(f"pattern not unique in {path}: {text.count(old)} matches")
    p.write_text(text.replace(old, new, 1))


headers = [
    "src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h",
    "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h",
]
for path in headers:
    replace_once(
        path,
        "        if (!mekf_) return;\n        if (!(dt > 0.0f) || !std::isfinite(dt)) return;\n",
        "        if (!mekf_) return;\n        if (!(dt > 0.0f) || !std::isfinite(dt)) return;\n\n"
        "        // Predictable online schedule: commit only coefficients staged\n"
        "        // after the previous IMU sample. The current sample can update\n"
        "        // the tuner later in this function, but it cannot choose the\n"
        "        // model/gain schedule used by its own Kalman innovation.\n"
        "        apply_pending_online_tune_();\n",
    )
    replace_once(
        path,
        "    bool enable_clamp_ = true;\n    bool enable_tuner_ = true;\n",
        "    bool enable_clamp_ = true;\n    bool enable_tuner_ = true;\n    bool online_tune_apply_pending_ = false;\n",
    )
    replace_once(
        path,
        "    void enableTuner(bool flag = true) {\n        enable_tuner_ = flag;\n    }\n",
        "    void enableTuner(bool flag = true) {\n"
        "        enable_tuner_ = flag;\n"
        "        if (!flag) online_tune_apply_pending_ = false;\n"
        "    }\n",
    )
    replace_once(
        path,
        "        enable_tuner_ = false;\n",
        "        enable_tuner_ = false;\n        online_tune_apply_pending_ = false;\n",
    )
    replace_once(
        path,
        "        last_adapt_time_sec_ = time_;\n        last_aw_cov_sync_sec_ = time_;\n",
        "        last_adapt_time_sec_ = time_;\n        last_aw_cov_sync_sec_ = time_;\n        online_tune_apply_pending_ = false;\n",
    )

replace_once(
    "src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h",
    "    // sync_covariance is set only by discrete reconfiguration events. The\n",
    "    // Apply the online tuner output only at the next IMU-sample boundary.\n"
    "    // adapt_mekf() may consume y_k and smooth its candidate during step k,\n"
    "    // but this function runs before y_{k+1} reaches the MEKF. Therefore the\n"
    "    // active schedule at step k+1 is measurable with respect to data through k.\n"
    "    void apply_pending_online_tune_() {\n"
    "        if (!online_tune_apply_pending_ || !mekf_) return;\n"
    "        apply_ou_tune_(false);\n"
    "        if (startup_stage_ == StartupStage::Live && enable_linear_block_) {\n"
    "            apply_R_p0_tune_();\n"
    "            apply_R_v0_tune_();\n"
    "        }\n"
    "        online_tune_apply_pending_ = false;\n"
    "    }\n\n"
    "    // sync_covariance is set only by discrete reconfiguration events. The\n",
)
replace_once(
    "src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h",
    "        if (time_ - last_adapt_time_sec_ > adapt_every_secs_) {\n"
    "            if (tuner_.isFreqReady()) {\n"
    "                apply_ou_tune_(false);\n"
    "            }\n"
    "            if (startup_stage_ == StartupStage::Live && enable_linear_block_) {\n"
    "                apply_R_p0_tune_();\n"
    "                apply_R_v0_tune_();\n"
    "            }\n"
    "            last_adapt_time_sec_ = time_;\n"
    "        }\n",
    "        if (time_ - last_adapt_time_sec_ > adapt_every_secs_) {\n"
    "            // y_k may change the smoothed candidate here, but the MEKF keeps\n"
    "            // the schedule that was active before y_k arrived. Commit this\n"
    "            // candidate at the beginning of updateTime(k+1).\n"
    "            online_tune_apply_pending_ = true;\n"
    "            last_adapt_time_sec_ = time_;\n"
    "        }\n",
)

replace_once(
    "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h",
    "    // sync_covariance is set only by discrete reconfiguration events. The\n",
    "    // Apply the online tuner output only at the next IMU-sample boundary.\n"
    "    // adapt_mekf() may consume y_k and smooth its candidate during step k,\n"
    "    // but this function runs before y_{k+1} reaches the MEKF. Therefore the\n"
    "    // active schedule at step k+1 is measurable with respect to data through k.\n"
    "    void apply_pending_online_tune_() {\n"
    "        if (!online_tune_apply_pending_ || !mekf_) return;\n"
    "        apply_ou_tune_(false);\n"
    "        if (startup_stage_ == StartupStage::Live && enable_linear_block_) {\n"
    "            apply_RS_tune_();\n"
    "        }\n"
    "        online_tune_apply_pending_ = false;\n"
    "    }\n\n"
    "    // sync_covariance is set only by discrete reconfiguration events. The\n",
)
replace_once(
    "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h",
    "        if (time_ - last_adapt_time_sec_ > adapt_every_secs_) {\n"
    "            if (tuner_.isFreqReady()) {\n"
    "                apply_ou_tune_(false);\n"
    "            }\n"
    "            if (startup_stage_ == StartupStage::Live && enable_linear_block_) {\n"
    "                apply_RS_tune_();\n"
    "            }\n"
    "            last_adapt_time_sec_ = time_;\n"
    "        }\n",
    "        if (time_ - last_adapt_time_sec_ > adapt_every_secs_) {\n"
    "            // y_k may change the smoothed candidate here, but the MEKF keeps\n"
    "            // the schedule that was active before y_k arrived. Commit this\n"
    "            // candidate at the beginning of updateTime(k+1).\n"
    "            online_tune_apply_pending_ = true;\n"
    "            last_adapt_time_sec_ = time_;\n"
    "        }\n",
)

replace_once(
    "doc/kalman_ou_iii/w3d-iss-stability.tex-part",
    "configuration analyzed here.\n\nAssume a Live interval on which\n",
    "configuration analyzed here.\n\n"
    "The implementation also makes the online schedule predictable with respect to\n"
    "the current Kalman innovation. Let $\\mathcal F_k$ denote the sigma-algebra\n"
    "generated by measurements through IMU sample $k$. At the beginning of step\n"
    "$k$, the filter commits only the tuning candidate staged after step $k-1$;\n"
    "the current sample $y_k$ is then processed by the MEKF and only afterward is\n"
    "it admitted to the tuner. Hence the active schedule satisfies\n"
    "\\begin{equation}\n"
    "\\vartheta_k\\in\\mathcal F_{k-1}.\n"
    "\\label{eq:iss-predictable-schedule}\n"
    "\\end{equation}\n"
    "The tuner may compute a new candidate from $y_k$, but that candidate is staged\n"
    "for step $k+1$. This one-sample scheduling delay is approximately\n"
    "\\SI{5}{ms} at the nominal \\SI{200}{Hz} IMU rate, negligible relative to the\n"
    "adaptation horizons, and supplies the usual predictability condition convenient\n"
    "for the strict stochastic Kalman argument. The deterministic UES/ISS result\n"
    "needs only the exogenous bounded schedule stated above.\n\n"
    "Assume a Live interval on which\n",
)
replace_once(
    "doc/kalman_ou_iii/w3d-fus-methods.tex-part",
    "through the filter's cross-covariances; see Sec.~\\ref{app:iss-tuner}.\n\nThe tuner input frequency is",
    "through the filter's cross-covariances; see Sec.~\\ref{app:iss-tuner}.\n\n"
    "Online coefficients are additionally staged by one IMU sample. At step $k$ the\n"
    "MEKF uses only $\\vartheta_k$ committed before the current innovation; after\n"
    "$y_k$ has been processed, the tuner may update its smoothed candidate, but that\n"
    "candidate becomes active only at the beginning of step $k+1$. Therefore\n"
    "$\\vartheta_k\\in\\mathcal F_{k-1}$ for the measurement filtration\n"
    "$\\mathcal F_k$. The delay is one sample (about \\SI{5}{ms} at\n"
    "\\SI{200}{Hz}), far below the adaptation time scales.\n\n"
    "The tuner input frequency is",
)
replace_once(
    "doc/kalman_ou_ii/kalman_ou_ii-charts.tex",
    "two need not agree, and in a developed sea they do not.\n\n% PM Stokes\n",
    "two need not agree, and in a developed sea they do not.\n\n"
    "The online operating point is applied with a one-IMU-sample scheduling delay.\n"
    "At step $k$ the OU--II Kalman recursion uses only the coefficients committed\n"
    "before the current innovation. The tuner may consume the current IMU sample\n"
    "and update its smoothed candidate afterward, but that candidate becomes active\n"
    "only at step $k+1$. Thus, for the measurement filtration $\\mathcal F_k$,\n"
    "the active schedule obeys $\\vartheta_k\\in\\mathcal F_{k-1}$. At the nominal\n"
    "\\SI{200}{Hz} IMU rate the added delay is about \\SI{5}{ms}, negligible beside\n"
    "the adaptation horizons.\n\n% PM Stokes\n",
)
replace_once(
    "doc/kalman_ou_ii/kalman_ou_ii-charts.tex",
    "\\usepackage{microtype}\n",
    "\\usepackage{microtype}\n\\usepackage{siunitx}\n",
)

Path("tests/kalman_ou_ii/tuner_schedule-test.cpp").write_text(r'''#define EIGEN_NON_ARDUINO
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#define private public
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"
#undef private

const float g_std = 9.80665f;
using Filter = SeaStateFusionFilter_OU_II<TrackerType::KALMANF>;

int main() {
    Filter f(false);
    f.initialize(Eigen::Vector3f::Constant(0.0148f),
                 Eigen::Vector3f::Constant(0.00157f),
                 Eigen::Vector3f::Constant(0.25f));
    const Eigen::Vector3f gyro = Eigen::Vector3f::Zero();
    const Eigen::Vector3f acc(0.0f, 0.0f, -g_std);
    f.initialize_from_acc(acc);
    f.startup_stage_ = Filter::StartupStage::Live;
    f.mekf_->set_linear_block_enabled(true);
    f.time_ = 1.0;
    f.last_adapt_time_sec_ = 0.0;
    f.adapt_every_secs_ = 0.05f;

    const float before = f.getTauApplied();
    f.adapt_mekf(0.1f, 3.0f, 0.8f, 1.2f, 0.7f);
    const float staged = f.tune_.tau_applied;
    if (!f.online_tune_apply_pending_) {
        std::cerr << "FAIL: adaptation was not staged\n";
        return 1;
    }
    if (std::fabs(f.getTauApplied() - before) > 1e-6f) {
        std::cerr << "FAIL: current-sample tuner output changed the active OU-II schedule\n";
        return 1;
    }
    if (std::fabs(staged - before) < 1e-5f) {
        std::cerr << "FAIL: test did not create a distinct staged candidate\n";
        return 1;
    }

    f.enable_tuner_ = false;
    f.updateTime(0.005f, gyro, acc);
    if (std::fabs(f.getTauApplied() - staged) > 1e-6f) {
        std::cerr << "FAIL: staged OU-II schedule was not committed at next sample\n";
        return 1;
    }
    if (f.online_tune_apply_pending_) {
        std::cerr << "FAIL: pending flag survived the commit\n";
        return 1;
    }
    std::cout << "OU-II predictable tuner scheduling passed\n";
    return 0;
}
''')

Path("tests/kalman_ou_iii/tuner_schedule-test.cpp").write_text(r'''#define EIGEN_NON_ARDUINO
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private

const float g_std = 9.80665f;
using Filter = SeaStateFusionFilter_OU_III<TrackerType::KALMANF>;

int main() {
    Filter f(false);
    f.initialize(Eigen::Vector3f::Constant(0.0148f),
                 Eigen::Vector3f::Constant(0.00157f),
                 Eigen::Vector3f::Constant(0.25f));
    const Eigen::Vector3f gyro = Eigen::Vector3f::Zero();
    const Eigen::Vector3f acc(0.0f, 0.0f, -g_std);
    f.initialize_from_acc(acc);
    f.startup_stage_ = Filter::StartupStage::Live;
    f.mekf_->set_linear_block_enabled(true);
    f.time_ = 1.0;
    f.last_adapt_time_sec_ = 0.0;
    f.adapt_every_secs_ = 0.05f;

    const float before = f.getTauApplied();
    f.adapt_mekf(0.1f, 3.0f, 0.8f, 1.2f);
    const float staged = f.tune_.tau_applied;
    if (!f.online_tune_apply_pending_) {
        std::cerr << "FAIL: adaptation was not staged\n";
        return 1;
    }
    if (std::fabs(f.getTauApplied() - before) > 1e-6f) {
        std::cerr << "FAIL: current-sample tuner output changed the active OU-III schedule\n";
        return 1;
    }
    if (std::fabs(staged - before) < 1e-5f) {
        std::cerr << "FAIL: test did not create a distinct staged candidate\n";
        return 1;
    }

    f.enable_tuner_ = false;
    f.updateTime(0.005f, gyro, acc);
    if (std::fabs(f.getTauApplied() - staged) > 1e-6f) {
        std::cerr << "FAIL: staged OU-III schedule was not committed at next sample\n";
        return 1;
    }
    if (f.online_tune_apply_pending_) {
        std::cerr << "FAIL: pending flag survived the commit\n";
        return 1;
    }
    std::cout << "OU-III predictable tuner scheduling passed\n";
    return 0;
}
''')

replace_once(
    "tests/kalman_ou_ii/Makefile",
    "TARGETS = kalman_ou_ii-sim\n",
    "TARGETS = kalman_ou_ii-sim tuner_schedule-test\n",
)
replace_once(
    "tests/kalman_ou_ii/Makefile",
    "kalman_ou_ii-sim: kalman_ou_ii-sim.o W3dSimCommon.o\n\t$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)\n\n",
    "kalman_ou_ii-sim: kalman_ou_ii-sim.o W3dSimCommon.o\n\t$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)\n\n"
    "tuner_schedule-test: tuner_schedule-test.o\n\t$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)\n\n",
)
replace_once(
    "tests/kalman_ou_ii/run_tests.sh",
    "./kalman_ou_ii-sim\n",
    "./kalman_ou_ii-sim\n./tuner_schedule-test\n",
)
replace_once(
    "tests/kalman_ou_iii/Makefile",
    "TARGETS = kalman_ou_iii-sim kalman_ou_common-test aw_covariance_policy-test acc_bias_ou-test channel_freeze-test wave_period-test mag_hard_iron-test tuner_coupling-test\n",
    "TARGETS = kalman_ou_iii-sim kalman_ou_common-test aw_covariance_policy-test acc_bias_ou-test channel_freeze-test wave_period-test mag_hard_iron-test tuner_coupling-test tuner_schedule-test\n",
)
replace_once(
    "tests/kalman_ou_iii/Makefile",
    "mag_hard_iron-test: mag_hard_iron-test.o\n\t$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)\n\n",
    "mag_hard_iron-test: mag_hard_iron-test.o\n\t$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)\n\n"
    "tuner_schedule-test: tuner_schedule-test.o\n\t$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)\n\n",
)
replace_once(
    "tests/kalman_ou_iii/run_tests.sh",
    "./tuner_coupling-test\n",
    "./tuner_coupling-test\n./tuner_schedule-test\n",
)
replace_once(
    "tests/kalman_ou_iii/tuner_coupling-test.cpp",
    "// three implementation facts that statement depends on, so none of them can\n// drift silently.\n",
    "// implementation facts that statement depends on, so none of them can\n"
    "// drift silently. State exogeneity is checked here; the separate\n"
    "// tuner_schedule-test pins the one-sample predictability boundary.\n",
)
