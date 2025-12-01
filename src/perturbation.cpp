#include "MLP.hpp"

#include <stack>

void MLP::double_bridge(Solution &solution) {
    size_t limit = (m_instance.get_dimension() + 9) / 10;
    
    std::uniform_int_distribution<size_t> distrib(1, limit);
    size_t first_size = distrib(m_rng);
    size_t second_size = distrib(m_rng);

    distrib = std::uniform_int_distribution<size_t>(1, m_instance.get_dimension() - first_size - second_size);
    size_t first_l = distrib(m_rng);
    size_t first_r = first_l + first_size - 1;

    distrib = std::uniform_int_distribution<size_t>(first_r + 1, m_instance.get_dimension() - second_size);
    size_t second_l = distrib(m_rng);
    size_t second_r = second_l + second_size - 1;

    
    double delta = 0.0;
    Subsequence new_sequence;
    
    if (first_r + 1 == second_l) {
        concatenate_subsequences_inplace(
            new_sequence,
            solution.subseq_matrix[0][first_l - 1],
            solution.subseq_matrix[second_l][second_r]
        ); 

        concatenate_subsequences_inplace(
            new_sequence,
            solution.subseq_matrix[first_l][first_r]
        );

        concatenate_subsequences_inplace(
            new_sequence,
            solution.subseq_matrix[second_r + 1][solution.sequence.size() - 1]
        );

        delta = new_sequence.acumulated_cost - solution.objective;
    } else {
        concatenate_subsequences_inplace(
            new_sequence,
            solution.subseq_matrix[0][first_l - 1],
            solution.subseq_matrix[second_l][second_r]
        ); 

        concatenate_subsequences_inplace(
            new_sequence,
            solution.subseq_matrix[first_r + 1][second_l - 1]
        );

        concatenate_subsequences_inplace(
            new_sequence,
            solution.subseq_matrix[first_l][first_r]
        );

        concatenate_subsequences_inplace(
            new_sequence,
            solution.subseq_matrix[second_r + 1][solution.sequence.size() - 1]
        );

        delta = new_sequence.acumulated_cost - solution.objective;
    }
    
    solution.objective += delta;

    std::stack<size_t> temp;

    for (size_t i = 0; i < first_l; ++i) {
        temp.push(solution.sequence[i]);
    }
    
    for (size_t i = second_l; i <= second_r; ++i) {
        temp.push(solution.sequence[i]);
    }

    for (size_t i = first_r + 1; i < second_l; ++i) {
        temp.push(solution.sequence[i]);
    }

    for (size_t i = first_l; i <= first_r; ++i) {
        temp.push(solution.sequence[i]);
    }

    for (size_t i = second_r + 1; i <= m_instance.get_dimension(); ++i) {
        temp.push(solution.sequence[i]);
    }

    for (size_t i = 0; i <= m_instance.get_dimension(); ++i) {
        solution.sequence[m_instance.get_dimension() - i] = temp.top();
        temp.pop();
    }

    update_interval_subsequences(solution, first_l, second_r);
    assert(solution.test_feasibility(m_instance));
    assert(test_subsequences_feasibility(solution));
}