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
    
    if (first_r + 1 == second_l) {
        // Remove edges
        delta -= m_instance.get_distance(solution.sequence[first_l - 1], solution.sequence[first_l]);
        delta -= m_instance.get_distance(solution.sequence[first_r], solution.sequence[second_l]);
        delta -= m_instance.get_distance(solution.sequence[second_r], solution.sequence[second_r + 1]);

        // Add edges
        delta += m_instance.get_distance(solution.sequence[first_l - 1], solution.sequence[second_l]);
        delta += m_instance.get_distance(solution.sequence[second_r], solution.sequence[first_l]);
        delta += m_instance.get_distance(solution.sequence[first_r], solution.sequence[second_r + 1]);


    } else {
        // Remove edges
        delta -= m_instance.get_distance(solution.sequence[first_l - 1], solution.sequence[first_l]);
        delta -= m_instance.get_distance(solution.sequence[first_r], solution.sequence[first_r + 1]);
        delta -= m_instance.get_distance(solution.sequence[second_l - 1], solution.sequence[second_l]);
        delta -= m_instance.get_distance(solution.sequence[second_r], solution.sequence[second_r + 1]);
        
        // Add edges
        delta += m_instance.get_distance(solution.sequence[first_l - 1], solution.sequence[second_l]);
        delta += m_instance.get_distance(solution.sequence[second_r], solution.sequence[first_r + 1]);
        delta += m_instance.get_distance(solution.sequence[second_l - 1], solution.sequence[first_l]);
        delta += m_instance.get_distance(solution.sequence[first_r], solution.sequence[second_r + 1]);
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

    assert(solution.test_feasibility(m_instance));
}