// Non-promoting consistency check between the legal complete-SEA3 shipping
// point experiment and frozen conditional P3 in the zero-radius limit.
//
// Reuse the exact H18 execution machinery from the complete-word feasibility
// diagnostic.  This file adds no source, schedule, covariance, or filter path.
// It merely rescales the same covariance-whitened H18 basis perturbations and
// reports the worst tested rho as the radius tends to zero.
//
// H18 is intentionally used here: it requires no held->active hybrid release,
// so the test cannot be contaminated by the separate A21 gate timing.  The
// 0.1 s magnetic packet spacing in the included point harness is already an
// admissible asynchronous PE realization inside the declared <=1 s recurrence
// window.  A rho that remains above one as radius -> 0 would therefore signal
// a coordinate/attachment mismatch in the point diagnostic, not a nonlinear
// basin failure.

#define main ou3_complete_word_original_main_for_reuse
#include "ou3-sea3-complete-word-feasibility-test.cpp"
#undef main

#include <array>

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            std::cerr << "usage: ou3-sea3-h18-tangent-radius-sweep-test <validated-specific-force-z.txt>\n";
            return 2;
        }
        std::ifstream in(argv[1]);
        if (!in) throw std::runtime_error("cannot open validated complete-SEA3 source trace");
        float x;
        while (in >> x) {
            if (!std::isfinite(x)) throw std::runtime_error("nonfinite source input");
            sfz.push_back(x);
        }
        if (sfz.size() != static_cast<std::size_t>(TOTAL)) {
            throw std::runtime_error("complete-word trace length mismatch");
        }

        constexpr std::array<float,5> radii{{0.1f, 0.03f, 0.01f, 0.003f, 0.001f}};
        std::cout << '{'
                  << "\"qualification\":\"OU3_COMPLETE_SEA3_H18_TANGENT_RADIUS_SWEEP_V1\","
                  << "\"canonical_source\":\"COMPLETE_SEA3_NORMAL_LIVE_WORD\","
                  << "\"same_shipping_H18_word\":true,"
                  << "\"trajectory_replay_used\":false,"
                  << "\"finite_harmonic_source_used\":false,"
                  << "\"P3_changed\":false,"
                  << "\"P4_promoted\":false,"
                  << "\"rows\":[";

        bool first = true;
        for (const float radius : radii) {
            auto probe = prepare_mode('H');
            const Eigen::MatrixXf P0 = covariance_mode(*probe, 18);
            Result best;
            best.rho = -std::numeric_limits<double>::infinity();
            for (int i = 0; i < 18; ++i) {
                Eigen::VectorXf d = local_delta(P0, i, 18);
                d *= radius / LOCAL_RADIUS;
                Result r = evaluate('H', d, "basis_" + std::to_string(i));
                if (r.rho > best.rho) best = std::move(r);
            }
            if (!first) std::cout << ',';
            first = false;
            std::cout << '{'
                      << "\"radius\":" << std::setprecision(9) << radius << ','
                      << "\"V0\":" << std::setprecision(17) << best.V0 << ','
                      << "\"V1\":" << best.V1 << ','
                      << "\"rho\":" << best.rho << ','
                      << "\"distance_to_one\":" << (1.0-best.rho) << ','
                      << "\"worst_basis\":\"" << best.label << "\","
                      << "\"direction\":" << vector_json(best.delta) << ','
                      << "\"imu_samples\":" << best.counts.imu << ','
                      << "\"mag_events\":" << best.counts.mag << ','
                      << "\"due_S_events\":" << best.counts.due_s << ','
                      << "\"effective_process_sigma_min\":" << best.counts.process_sigma_min
                      << '}';
        }
        std::cout << "]}\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "H18 tangent radius sweep failure: " << e.what() << '\n';
        return 1;
    }
}
