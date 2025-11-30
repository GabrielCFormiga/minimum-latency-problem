#include "MLP.hpp"

Solution MLP::GILS_RVND(const size_t max_iterations, const size_t max_ils_iterations) {
    Solution global_best;
    global_best.objective = INFINITY;

    std::uniform_int_distribution<size_t> distrib(0, alpha_values.size() - 1);

    for (size_t i = 0; i < max_iterations; ++i) {
        Solution s = randomized(alpha_values[distrib(m_rng)]);
        Solution local_best = s;

        size_t ils_iterations = 0;

        while (ils_iterations <= max_ils_iterations) {
            RVND(s);

            if (s.objective + EPS < local_best.objective) {
                local_best = s;
                ils_iterations = 0;
            }
            
            if (local_best.objective + EPS < global_best.objective) {
                global_best = local_best;
            }

            s = local_best;
            // double_bridge(s);
            ++ils_iterations;
        }
    }
    
    return global_best;
}