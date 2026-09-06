// Non-promoting finite-angle point sweep on the legal complete-SEA3 shipping word.
//
// The checked-in complete-word evaluator is itself bound to shipping's qref
// left-error coordinate. This harness only exposes that same evaluator
// at the declared 15/20/25/30 degree attitude candidates, for both H18 and A21.
// It introduces no source, schedule, covariance, or filter path.

#define main ou3_complete_word_original_main_for_reuse
#include "ou3-sea3-complete-word-feasibility-test.cpp"
#undef main

#include <array>

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            std::cerr << "usage: ou3-sea3-candidate-angle-sweep-test <validated-specific-force-z.txt>\n";
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

        constexpr std::array<float,4> angles{{15.0f,20.0f,25.0f,30.0f}};
        std::cout << '{'
                  << "\"qualification\":\"OU3_COMPLETE_SEA3_CANDIDATE_ANGLE_POINT_SWEEP_V1\","
                  << "\"canonical_source\":\"COMPLETE_SEA3_NORMAL_LIVE_WORD\","
                  << "\"same_continuum_source_across_H_release_A\":true,"
                  << "\"shipping_qref_left_error_coordinate\":true,"
                  << "\"trajectory_replay_used\":false,"
                  << "\"finite_harmonic_source_used\":false,"
                  << "\"independent_schedule_used\":false,"
                  << "\"P3_changed\":false,"
                  << "\"P4_promoted\":false,"
                  << "\"rows\":[";
        bool first = true;
        for (char mode : {'H','A'}) {
            const int dim = mode == 'H' ? 18 : 21;
            for (float deg : angles) {
                double worst = -std::numeric_limits<double>::infinity();
                std::string worst_label;
                double best = std::numeric_limits<double>::infinity();
                for (int axis = 0; axis < 3; ++axis) {
                    for (float sign : {-1.0f,1.0f}) {
                        const float signed_deg = sign * deg;
                        const std::string label = std::string("attitude_axis_")
                            + std::to_string(axis) + "_deg_" + std::to_string(signed_deg);
                        Result r = evaluate(mode, attitude_delta(dim, axis, signed_deg), label);
                        if (r.rho > worst) { worst = r.rho; worst_label = r.label; }
                        best = std::min(best, r.rho);
                    }
                }
                if (!first) std::cout << ',';
                first = false;
                std::cout << '{'
                          << "\"mode\":\"" << mode << "\","
                          << "\"angle_deg\":" << std::setprecision(9) << deg << ','
                          << "\"worst_rho\":" << std::setprecision(17) << worst << ','
                          << "\"best_rho\":" << best << ','
                          << "\"distance_to_one\":" << (1.0-worst) << ','
                          << "\"worst_direction\":\"" << worst_label << "\""
                          << '}';
            }
        }
        std::cout << "]}\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "candidate-angle sweep failure: " << e.what() << '\n';
        return 1;
    }
}
