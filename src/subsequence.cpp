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
            concatenate_subsequences_inplace(
                solution.subseq_matrix[i][j],
                solution.subseq_matrix[i][j - 1],
                solution.subseq_matrix[j][j]
            );
        }
    }

    // <- Subsequences
    for (size_t i = solution.sequence.size() - 1; i > 0; --i) {
        for (size_t j = i - 1; j >= 0; --j) {
            concatenate_subsequences_inplace(
                solution.subseq_matrix[i][j],
                solution.subseq_matrix[i][j + 1],
                solution.subseq_matrix[j][j]
            );
            if (j == 0) break;
        }
    }
}

void MLP::update_interval_subsequences(Solution &solution, size_t l, size_t r) {
    if (l > r) std::swap(l, r);

    // Single element subsequences
    for (size_t i = l; i <= r; ++i) {
        solution.subseq_matrix[i][i].length = (i > 0);
        solution.subseq_matrix[i][i].cost = 0.0;
        solution.subseq_matrix[i][i].acumulated_cost = 0.0;
        solution.subseq_matrix[i][i].first = solution.sequence[i];
        solution.subseq_matrix[i][i].last = solution.sequence[i];
    }

    // -> Subsequences
    for (size_t i = 0; i <= r; ++i) {
        for (size_t j = std::max(l, i + 1); j < solution.sequence.size(); ++j) {
            concatenate_subsequences_inplace(
                solution.subseq_matrix[i][j],
                solution.subseq_matrix[i][j - 1],
                solution.subseq_matrix[j][j]
            );
        }
    }

    // <- Subsequences
    for (size_t i = solution.sequence.size() - 1; i >= l; --i) {
        for (size_t j = std::min(i - 1, r); j >= 0; --j) {
            concatenate_subsequences_inplace(
                solution.subseq_matrix[i][j],
                solution.subseq_matrix[i][j + 1],
                solution.subseq_matrix[j][j]
            );
            if (j == 0) break;
        }
    }
}

bool MLP::test_subsequences_feasibility(const Solution &solution) {
    std::vector<std::vector<Subsequence>> new_subseq_matrix(m_instance.get_dimension() + 1, std::vector<Subsequence>(m_instance.get_dimension() + 1));

    // Single element subsequences
    for (size_t i = 0; i < solution.sequence.size(); ++i) {
        new_subseq_matrix[i][i].length = (i > 0);
        new_subseq_matrix[i][i].cost = 0.0;
        new_subseq_matrix[i][i].acumulated_cost = 0.0;
        new_subseq_matrix[i][i].first = solution.sequence[i];
        new_subseq_matrix[i][i].last = solution.sequence[i];
    }

    // -> Subsequences
    for (size_t i = 0; i < solution.sequence.size() - 1; ++i) {
        for (size_t j = i + 1; j < solution.sequence.size(); ++j) {
            new_subseq_matrix[i][j] = concatenate_subsequences(
                new_subseq_matrix[i][j - 1],
                new_subseq_matrix[j][j]
            );
        }
    }

    // <- Subsequences
    for (size_t i = solution.sequence.size() - 1; i > 0; --i) {
        for (size_t j = i - 1; j >= 0; --j) {
            new_subseq_matrix[i][j] = concatenate_subsequences(
                new_subseq_matrix[i][j + 1],
                new_subseq_matrix[j][j]
            );
            if (j == 0) break;
        }
    }

    for (size_t i = 0; i < solution.sequence.size(); ++i) {
        for (size_t j = 0; j < solution.sequence.size(); ++j) {
            if (new_subseq_matrix[i][j].length != solution.subseq_matrix[i][j].length) return false;
            if (std::abs(new_subseq_matrix[i][j].cost - solution.subseq_matrix[i][j].cost) > EPS) return false;
            if (std::abs(new_subseq_matrix[i][j].acumulated_cost - solution.subseq_matrix[i][j].acumulated_cost) > EPS) return false;
            if (new_subseq_matrix[i][j].first != solution.subseq_matrix[i][j].first) return false;
            if (new_subseq_matrix[i][j].last != solution.subseq_matrix[i][j].last) return false;
        }
    }

    return true;
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

void MLP::concatenate_subsequences_inplace(Subsequence &a, const Subsequence &b) {
    a.acumulated_cost = a.acumulated_cost + b.length * (a.cost + m_instance.get_distance(a.last, b.first)) + b.acumulated_cost;
    a.cost = a.cost + m_instance.get_distance(a.last, b.first) + b.cost;
    a.length = a.length + b.length;
    a.last = b.last;
}

void MLP::concatenate_subsequences_inplace(Subsequence &ret, Subsequence &a, const Subsequence &b) {
    ret.acumulated_cost = a.acumulated_cost + b.length * (a.cost + m_instance.get_distance(a.last, b.first)) + b.acumulated_cost;
    ret.cost = a.cost + m_instance.get_distance(a.last, b.first) + b.cost;
    ret.length = a.length + b.length;
    ret.first = a.first;
    ret.last = b.last;
}