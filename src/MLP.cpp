#include "MLP.hpp"
#include <chrono>

MLP::MLP(Instance &instance, uint64_t seed) : m_instance(instance) {
    if (seed == 0) {
        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }
    
    for (double a = 0.0; a <= 0.25 + EPS; a += 0.01) {
        alpha_values.push_back(a);
    }
    
    m_rng.seed(seed);

    // std::cout << "RNG seed: " << seed << std::endl;
}