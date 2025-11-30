#ifndef MLP_HPP
#define MLP_HPP

#include <instance.hpp>
#include <solution.hpp>
#include <subsequence.hpp>
#include <random>
#include <cassert>

#define EPS 0.0001

class MLP {
    private:
        Instance m_instance;
        std::vector<double> alpha_values;
        std::mt19937_64 m_rng;

    public:
        enum class Neighborhood : uint8_t { SWAP, _2_OPT, REINSERTION, OR_OPT_2, OR_OPT_3 };

        MLP(Instance &instance, uint64_t seed = 0);
        
        // Constructive heuristics
        Solution randomized(double alpha);

        // Local Search Procedures
        void RVND(Solution &solution);
        
        // Neighborhoods
        bool best_improvement_swap(Solution &solution);
        bool best_improvement_2_opt(Solution &solution);
        bool best_improvement_or_opt(Solution &solution, size_t segment_size);
        
        // Methaheuristics
        Solution GILS_RVND(const size_t max_iterations, const size_t max_ils_iterations);
        
        // Perturbations
        void double_bridge(Solution &solution);

        // Subsequences management
        void update_all_subsequences(Solution &solution);

        Subsequence concatenate_subsequences(const Subsequence &a, const Subsequence &b);
};

#endif