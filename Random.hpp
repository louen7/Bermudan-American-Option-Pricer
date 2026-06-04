// ============================================================================
//  Random.hpp - Random Number Generator & Normal CDF
//  Mersenne Twister engine with uniform/normal sampling
//  Abramowitz-Stegun approximation for the standard normal CDF
// ============================================================================

#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <random>

class Random {
private:
    std::mt19937 generator;
    std::uniform_real_distribution<double> uniform_dist;

public:
    Random(unsigned int seed = 12345);

    // Draw a uniform random number in [0, 1)
    double uniforme();

    // Draw a standard normal random number (Box-Muller transform)
    double normale();

    // Standard normal cumulative distribution function (CDF)
    static double abramowitz_stegun(double x);
};

#endif // RANDOM_HPP
