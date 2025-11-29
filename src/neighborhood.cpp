#include "MLP.hpp"

#include <algorithm>

// Always has sequence[0] == sequence.back() == 1

bool MLP::best_improvement_swap(Solution &solution) {
    double best_delta = 0.0;
    size_t best_i = 0, best_j = 0;

    for (size_t i = 1; i < solution.sequence.size() - 2; ++i) {
        for (size_t j = i + 1; j < solution.sequence.size() - 1; ++j) {
            double delta = 0;

            if (j == i + 1) {
                // Remove edges
                delta -= m_instance.get_distance(solution.sequence[i - 1], solution.sequence[i]);
                delta -= m_instance.get_distance(solution.sequence[j], solution.sequence[j + 1]);

                // Add edges
                delta += m_instance.get_distance(solution.sequence[i - 1], solution.sequence[j]);
                delta += m_instance.get_distance(solution.sequence[i], solution.sequence[j + 1]);
            } else {
                // Remove edges
                delta -= m_instance.get_distance(solution.sequence[i - 1], solution.sequence[i]);
                delta -= m_instance.get_distance(solution.sequence[i], solution.sequence[i + 1]);
                delta -= m_instance.get_distance(solution.sequence[j - 1], solution.sequence[j]);
                delta -= m_instance.get_distance(solution.sequence[j], solution.sequence[j + 1]);

                // Add edges
                delta += m_instance.get_distance(solution.sequence[i - 1], solution.sequence[j]);
                delta += m_instance.get_distance(solution.sequence[j], solution.sequence[i + 1]);
                delta += m_instance.get_distance(solution.sequence[j - 1], solution.sequence[i]);
                delta += m_instance.get_distance(solution.sequence[i], solution.sequence[j + 1]);
            }

            if (delta + EPS < best_delta) {
                best_delta = delta;
                best_i = i;
                best_j = j;
            }
        }
    }

    if (best_delta + EPS < 0) {
        std::swap(solution.sequence[best_i], solution.sequence[best_j]);
        solution.objective += best_delta;
        assert(solution.test_feasibility(m_instance));
    }

    return best_delta + EPS < 0;
}

bool MLP::best_improvement_2_opt(Solution &solution) {
    double best_delta = 0.0;
    size_t best_i = 0, best_j = 0;

    for (size_t i = 1; i < solution.sequence.size() - 3; ++i) {
        for (size_t j = i + 2; j < solution.sequence.size() - 1; ++j) {
            double delta = 0;

            // Remove edges
            delta -= m_instance.get_distance(solution.sequence[i - 1], solution.sequence[i]);
            delta -= m_instance.get_distance(solution.sequence[j], solution.sequence[j + 1]);

            // Add edges
            delta += m_instance.get_distance(solution.sequence[i - 1], solution.sequence[j]);
            delta += m_instance.get_distance(solution.sequence[i], solution.sequence[j + 1]);

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
        assert(solution.test_feasibility(m_instance));
    }

    return best_delta + EPS < 0;
}

bool MLP::best_improvement_or_opt(Solution &solution, size_t segment_size) {
    double best_delta = 0.0;
    size_t best_i = 0, best_j = 0;

    for (size_t i = 1; i < solution.sequence.size() - segment_size; ++i) {
        for (size_t j = 1; j < solution.sequence.size(); ++j) {
            if (j >= i && j <= i + segment_size) continue;

            double delta = 0;

            // Remove edges
            delta -= m_instance.get_distance(solution.sequence[i - 1], solution.sequence[i]);
            delta -= m_instance.get_distance(solution.sequence[i + segment_size - 1], solution.sequence[i + segment_size]);
            delta -= m_instance.get_distance(solution.sequence[j - 1], solution.sequence[j]);

            // Add edges
            delta += m_instance.get_distance(solution.sequence[i - 1], solution.sequence[i + segment_size]);
            delta += m_instance.get_distance(solution.sequence[j - 1], solution.sequence[i]);
            delta += m_instance.get_distance(solution.sequence[i + segment_size - 1], solution.sequence[j]);

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
        assert(solution.test_feasibility(m_instance));
    }

    return best_delta + EPS < 0;
}
