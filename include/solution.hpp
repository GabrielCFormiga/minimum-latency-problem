#ifndef SOLUTION_HPP
#define SOLUTION_HPP

#include <vector>

#include "instance.hpp"
#include "subsequence.hpp"

struct Solution {
    std::vector<size_t> sequence;
  
    std::vector<std::vector<Subsequence>> subseq_matrix;

    double objective;

    Solution() = default;
    Solution(const Instance &instance);

    bool test_feasibility(const Instance &instance);
  
    void print_sequence() const;
};

#endif