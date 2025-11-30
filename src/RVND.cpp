#include "MLP.hpp"

void MLP::RVND(Solution &solution) {
    std::vector<MLP::Neighborhood> neighborhoods = {
        MLP::Neighborhood::SWAP,
        // MLP::Neighborhood::_2_OPT,
        // MLP::Neighborhood::REINSERTION, 
        // MLP::Neighborhood::OR_OPT_2, 
        // MLP::Neighborhood::OR_OPT_3
    };

    bool improved;
    size_t neighborhood;

    while (!neighborhoods.empty()) {
        std::uniform_int_distribution<size_t> distrib(0, neighborhoods.size() - 1);
        neighborhood = distrib(m_rng);

        switch (neighborhoods[neighborhood]) {
            case MLP::Neighborhood::SWAP:
                improved = best_improvement_swap(solution);
                break;
            case MLP::Neighborhood::_2_OPT:
                improved = best_improvement_2_opt(solution);
                break;
            case MLP::Neighborhood::REINSERTION:
                improved = best_improvement_or_opt(solution, 1);
                break;
            case MLP::Neighborhood::OR_OPT_2:
                improved = best_improvement_or_opt(solution, 2);
                break;
            case MLP::Neighborhood::OR_OPT_3:
                improved = best_improvement_or_opt(solution, 3);
                break;
        }

        if (improved) {
            neighborhoods = {
                MLP::Neighborhood::SWAP,
                // MLP::Neighborhood::_2_OPT,
                // MLP::Neighborhood::REINSERTION,
                // MLP::Neighborhood::OR_OPT_2,
                // MLP::Neighborhood::OR_OPT_3
            };
        } else {
            neighborhoods.erase(neighborhoods.begin() + neighborhood);
        }
    }
    
}