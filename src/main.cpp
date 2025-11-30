#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <chrono>

#include "MLP.hpp"
#include "instance.hpp"
#include "solution.hpp"

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <instance_file> [seed] [max_iterations] [max_ils_iterations]" << '\n';
        std::cerr << "  instance_file: MLP instance file path" << '\n';
        std::cerr << "  seed: Random seed (default: 0 = use time-based seed)" << '\n';
        std::cerr << "  max_iterations: GRASP iterations (default: 50)" << '\n';
        std::cerr << "  max_ils_iterations: ILS iterations (default: (|V| >= 150) ? |V|/2 : |V|)" << '\n';
        return 1;
    }

    // Read instance
    Instance instance(argc, argv[1]);
    instance.read();

    // Argparse
    uint64_t seed = (argc > 2) ? std::stoull(argv[2]) : 0;
    size_t max_iterations = (argc > 3) ? std::stoull(argv[3]) : 50;
    size_t max_ils_iterations = (argc > 4) ? std::stoull(argv[4]) : (instance.get_dimension() >= 150 ? instance.get_dimension() / 2 : instance.get_dimension());

    // Display configuration
    std::cout << std::left;
    std::cout << std::string(40, '=') << '\n';
    std::cout << "GILS-RVND Configuration:" << '\n';
    std::cout << std::setw(20) << "Instance:" << instance.get_name() << '\n';
    std::cout << std::setw(20) << "GRASP iterations:" << max_iterations << '\n';
    std::cout << std::setw(20) << "ILS iterations:" << max_ils_iterations << '\n';
    std::cout << std::string(40, '=') << '\n';

    MLP MLP(instance, seed);

    Solution s = MLP.GILS_RVND(max_iterations, max_ils_iterations);

    // Display solution
    std::cout << std::string(40, '=') << '\n';
    std::cout << "Solution:" << '\n';   
    std::cout << "Objective: " << s.objective << '\n';
    std::cout << "Sequence: " << '\n';
    s.print_sequence();
    std::cout << std::string(40, '=') << '\n';

    return 0;
}
