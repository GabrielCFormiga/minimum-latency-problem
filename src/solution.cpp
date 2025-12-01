#include "solution.hpp"

#include <iomanip>
#include <iostream>

Solution::Solution(const Instance &instance) {
    sequence.reserve(instance.get_dimension() + 1);
    subseq_matrix.resize(instance.get_dimension() + 1, std::vector<Subsequence>(instance.get_dimension() + 1));
    objective = 0.0;
}

void Solution::print_subsequence_matrix() const {
    for (size_t i = 0; i < subseq_matrix.size(); ++i) {
        for (size_t j = 0; j < subseq_matrix[i].size(); ++j) {
            std::cout << std::setw(4) << subseq_matrix[i][j].acumulated_cost << ' ';
        }
        std::cout << '\n';
    }
}

void Solution::print_sequence() const {
    for (size_t i = 0; i < sequence.size(); ++i) {
        std::cout << sequence[i] << (i == sequence.size() - 1 ? "\n" : " -> "); 
    }
}

bool Solution::test_feasibility(const Instance &instance) {
    if (sequence.size() != instance.get_dimension() + 1) {
        return false;
    }

    if (sequence[0] != sequence[sequence.size() - 1]) return false;

    if (sequence[0] != 1) return false; // Only for MLP with fixed deposit/leader position

    double new_objective = 0;
    double prefix = 0.0;
    std::vector<bool> used(instance.get_dimension() + 1, false);

    for (size_t i = 0; i < sequence.size() - 1; i++) {
        if (sequence[i] < 1 || sequence[i] > instance.get_dimension()) return false;
        if (used[sequence[i]] == true) return false;
        used[sequence[i]] = true;
        prefix += instance.get_distance(sequence[i], sequence[i + 1]);
        new_objective += prefix;
    }
    
    return objective == new_objective;
}