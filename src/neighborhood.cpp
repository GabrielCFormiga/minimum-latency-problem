#include "MLP.hpp"

#include <algorithm>

bool MLP::best_improvement_swap(Solution &solution) {
    double best_delta = 0.0;
    size_t best_i = 0, best_j = 0;

    /* Only for MLP with variable deposit/leader position 
    // New leader swaps
    for (size_t j = 1; j < solution.sequence.size() - 1; ++j) {
        double delta = 0.0;

        if (j == 1) {
            // Remove edges
            delta -= m_instance.get_distance(solution.sequence[solution.sequence.size() - 2], solution.sequence[solution.sequence.size() - 1]);
            delta -= m_instance.get_distance(solution.sequence[j], solution.sequence[j + 1]) * (solution.sequence.size() - (j + 1));
        
            // Add edges
            delta += m_instance.get_distance(solution.sequence[0], solution.sequence[j + 1]) * (solution.sequence.size() - (j + 1));
            delta += m_instance.get_distance(solution.sequence[solution.sequence.size() - 2], solution.sequence[j]);
        } else if (j == solution.sequence.size() - 2) {
            // Remove edges
            delta -= m_instance.get_distance(solution.sequence[0], solution.sequence[1]) * (solution.sequence.size() - 1);
            delta -= m_instance.get_distance(solution.sequence[j - 1], solution.sequence[j]) * (solution.sequence.size() - j);
        
            // Add edges
            delta += m_instance.get_distance(solution.sequence[j], solution.sequence[1]) * (solution.sequence.size() - 1);
            delta += m_instance.get_distance(solution.sequence[j - 1], solution.sequence[0]) * (solution.sequence.size() - j);
        } else {
            // Remove edges
            delta -= m_instance.get_distance(solution.sequence[0], solution.sequence[1]) * (solution.sequence.size() - 1);
            delta -= m_instance.get_distance(solution.sequence[j - 1], solution.sequence[j]) * (solution.sequence.size() - j);
            delta -= m_instance.get_distance(solution.sequence[j], solution.sequence[j + 1]) * (solution.sequence.size() - (j + 1));
            delta -= m_instance.get_distance(solution.sequence[solution.sequence.size() - 2], solution.sequence[solution.sequence.size() - 1]); 

            // Add edges
            delta += m_instance.get_distance(solution.sequence[j], solution.sequence[1]) * (solution.sequence.size() - 1);
            delta += m_instance.get_distance(solution.sequence[j - 1], solution.sequence[0]) * (solution.sequence.size() - j);
            delta += m_instance.get_distance(solution.sequence[0], solution.sequence[j + 1]) * (solution.sequence.size() - (j + 1));
            delta += m_instance.get_distance(solution.sequence[solution.sequence.size() - 2], solution.sequence[j]);
        }

        if (delta + EPS < best_delta) {
            best_delta = delta;
            best_i = 0;
            best_j = j;
        }
    } */

    // Other swaps
    for (size_t i = 1; i < solution.sequence.size() - 2; ++i) {
        for (size_t j = i + 1; j < solution.sequence.size() - 1; ++j) {
            double delta = 0;

            if (j == i + 1) {
                // Remove edges
                delta -= m_instance.get_distance(solution.sequence[i - 1], solution.sequence[i]) * (solution.sequence.size() - i);
                delta -= m_instance.get_distance(solution.sequence[j], solution.sequence[j + 1]) * (solution.sequence.size() - (j + 1));

                // Add edges
                delta += m_instance.get_distance(solution.sequence[i - 1], solution.sequence[j]) * (solution.sequence.size() - i);
                delta += m_instance.get_distance(solution.sequence[i], solution.sequence[j + 1]) * (solution.sequence.size() - (j + 1));
            } else {
                // Remove edges
                delta -= m_instance.get_distance(solution.sequence[i - 1], solution.sequence[i]) * (solution.sequence.size() - i);
                delta -= m_instance.get_distance(solution.sequence[i], solution.sequence[i + 1]) * (solution.sequence.size() - (i + 1));
                delta -= m_instance.get_distance(solution.sequence[j - 1], solution.sequence[j]) * (solution.sequence.size() - j);
                delta -= m_instance.get_distance(solution.sequence[j], solution.sequence[j + 1]) * (solution.sequence.size() - (j + 1));

                // Add edges
                delta += m_instance.get_distance(solution.sequence[i - 1], solution.sequence[j]) * (solution.sequence.size() - i);
                delta += m_instance.get_distance(solution.sequence[j], solution.sequence[i + 1]) * (solution.sequence.size() - (i + 1));
                delta += m_instance.get_distance(solution.sequence[j - 1], solution.sequence[i]) * (solution.sequence.size() - j);
                delta += m_instance.get_distance(solution.sequence[i], solution.sequence[j + 1]) * (solution.sequence.size() - (j + 1));
            }

            if (delta + EPS < best_delta) {
                best_delta = delta;
                best_i = i;
                best_j = j;
            }
        }
    }

    if (best_delta + EPS < 0) {
        /* Only for MLP with variable deposit/leader position 
        if (best_i == 0) {
            solution.sequence[solution.sequence.size() - 1] = solution.sequence[best_j];
        } */
        std::swap(solution.sequence[best_i], solution.sequence[best_j]);
        solution.objective += best_delta;
        update_interval_subsequences(solution, best_i, best_j);
        assert(solution.test_feasibility(m_instance));
        assert(test_subsequences_feasibility(solution));
    }

    return best_delta + EPS < 0;
}

bool MLP::best_improvement_2_opt(Solution &solution) {
    double best_delta = 0.0;
    size_t best_i = 0, best_j = 0;

    for (size_t i = 1; i < solution.sequence.size() - 3; ++i) {
        for (size_t j = i + 2; j < solution.sequence.size() - 1; ++j) {
            Subsequence new_sequence = concatenate_subsequences(
                solution.subseq_matrix[0][i - 1], 
                solution.subseq_matrix[j][i]
            );
            new_sequence = concatenate_subsequences(
                new_sequence, 
                solution.subseq_matrix[j + 1][solution.sequence.size() - 1]
            );

            double delta = new_sequence.acumulated_cost - solution.objective;

            if (delta + EPS < best_delta) {
                best_delta = delta;
                best_i = i;
                best_j = j;
            }
        }
    }

    if (best_delta + EPS < 0) {
        reverse(solution.sequence.begin() + best_i, solution.sequence.begin() + best_j + 1);
        solution.objective += best_delta;
        update_interval_subsequences(solution, best_i, best_j);
        assert(solution.test_feasibility(m_instance));
        assert(test_subsequences_feasibility(solution));
    }

    return best_delta + EPS < 0;
}

bool MLP::best_improvement_or_opt(Solution &solution, size_t segment_size) {
    double best_delta = 0.0;
    size_t best_i = 0, best_j = 0;

    for (size_t i = 1; i < solution.sequence.size() - segment_size; ++i) {
        for (size_t j = 1; j < solution.sequence.size(); ++j) {
            if (j >= i && j <= i + segment_size) continue;

            double delta = 0.0;

            if (j < i) {
                Subsequence new_sequence = concatenate_subsequences(
                    solution.subseq_matrix[0][j - 1],
                    solution.subseq_matrix[i][i + segment_size - 1]
                );
    
                new_sequence = concatenate_subsequences(
                    new_sequence,
                    solution.subseq_matrix[j][i - 1]
                );

                new_sequence = concatenate_subsequences(
                    new_sequence,
                    solution.subseq_matrix[i + segment_size][solution.sequence.size() - 1]
                );

                delta = new_sequence.acumulated_cost - solution.objective;
            } else {
                Subsequence new_sequence = concatenate_subsequences(
                    solution.subseq_matrix[0][i - 1],
                    solution.subseq_matrix[i + segment_size][j - 1]
                );

                new_sequence = concatenate_subsequences(
                    new_sequence,
                    solution.subseq_matrix[i][i + segment_size - 1]
                );

                new_sequence = concatenate_subsequences(
                    new_sequence,
                    solution.subseq_matrix[j][solution.sequence.size() - 1]
                );

                delta = new_sequence.acumulated_cost - solution.objective;
            }

            if (delta + EPS < best_delta) {
                best_delta = delta;
                best_i = i;
                best_j = j;
            }
        }
    }

    if (best_delta + EPS < 0) {
        std::vector<size_t> segment(segment_size);

        for (size_t k = 0; k < segment_size; ++k) {
            segment[k] = solution.sequence[best_i + k];
        }

        if (best_i < best_j) {
            rotate(solution.sequence.begin() + best_i, solution.sequence.begin() + best_i + segment_size, solution.sequence.begin() + best_j);
        } else {
            rotate(solution.sequence.begin() + best_j, solution.sequence.begin() + best_i, solution.sequence.begin() + best_i + segment_size);
        }

        solution.objective += best_delta;
        update_interval_subsequences(solution, std::min(best_i, best_j), std::max(best_i + segment_size - 1, best_j));
        assert(solution.test_feasibility(m_instance));
        assert(test_subsequences_feasibility(solution));
    }

    return best_delta + EPS < 0;
}
