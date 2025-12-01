#include "MLP.hpp"
#include <algorithm>

Solution MLP::randomized(double alpha) {
    Solution solution(m_instance);
        
    solution.sequence.push_back(1);

    std::vector<size_t> candidates;

    for (size_t i = 2; i <= m_instance.get_dimension(); ++i) {
        candidates.push_back(i);
    }

    size_t prev = 1;
    double prefix = 0.0;

    while (!candidates.empty()) {
        sort(candidates.begin(), candidates.end(), [&](size_t a, size_t b) {
            return m_instance.get_distance(prev, a) < m_instance.get_distance(prev, b);
        });

        size_t limit = std::min(static_cast<size_t>(floor(100 * alpha * candidates.size())) + 1, candidates.size());
        std::uniform_int_distribution<size_t> distrib(0, limit - 1);
        size_t ind = distrib(m_rng);

        prefix += m_instance.get_distance(prev, candidates[ind]);
        solution.objective += prefix;
        solution.sequence.push_back(candidates[ind]);
        prev = candidates[ind];
        std::swap(candidates[ind], candidates[candidates.size() - 1]);
        candidates.pop_back();
    }

    solution.sequence.push_back(1);
    prefix += m_instance.get_distance(prev, 1);
    solution.objective += prefix;
    update_all_subsequences(solution);

    assert(solution.test_feasibility(m_instance));
    assert(test_subsequences_feasibility(solution));
    
    return solution;
}