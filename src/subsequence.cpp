#include "MLP.hpp"

void MLP::update_all_subsequences(Solution &solution) {
    // Single element subsequences
    for (size_t i = 0; i < solution.sequence.size(); ++i) {
        solution.subseq_matrix[i][i].length = (i > 0);
        solution.subseq_matrix[i][i].cost = 0.0;
        solution.subseq_matrix[i][i].acumulated_cost = 0.0;
        solution.subseq_matrix[i][i].first = solution.sequence[i];
        solution.subseq_matrix[i][i].last = solution.sequence[i];
    }

    // -> Subsequences
    for (size_t i = 0; i < solution.sequence.size() - 1; ++i) {
        for (size_t j = i + 1; j < solution.sequence.size(); ++j) {
            solution.subseq_matrix[i][j] = concatenate_subsequences(
                solution.subseq_matrix[i][j - 1],
                solution.subseq_matrix[j][j]
            );
        }
    }

    // <- Subsequences
    for (size_t i = solution.sequence.size() - 1; i > 0; --i) {
        for (size_t j = i - 1; j >= 0; --j) {
            solution.subseq_matrix[i][j] = concatenate_subsequences(
                solution.subseq_matrix[i][j + 1],
                solution.subseq_matrix[j][j]
            );
            if (j == 0) break;
        }
    }
}

Subsequence MLP::concatenate_subsequences(const Subsequence &a, const Subsequence &b) {
    Subsequence ret;
    ret.length = a.length + b.length;
    ret.cost = a.cost + m_instance.get_distance(a.last, b.first) + b.cost;
    ret.acumulated_cost = a.acumulated_cost + b.length * (a.cost + m_instance.get_distance(a.last, b.first)) + b.acumulated_cost;
    ret.first = a.first;
    ret.last = b.last;
    return ret;
}