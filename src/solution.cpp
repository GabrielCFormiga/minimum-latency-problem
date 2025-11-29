#include "solution.hpp"

#include <iostream>

Solution::Solution(const Instance &instance) {
    sequence.reserve(instance.get_dimension() + 1);
    objective = 0.0;
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

    if (sequence[0] != 1 || sequence[sequence.size() - 1] != 1) return false;

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