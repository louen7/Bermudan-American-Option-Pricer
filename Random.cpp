// ============================================================================
//  Random.cpp - Random Number Generator & Normal CDF Implementation
// ============================================================================

#include "Random.hpp"
#include <cmath>

// ---------------------------------------------------------------------------
//  Constructor - Initialize Mersenne Twister with a given seed
// ---------------------------------------------------------------------------
Random::Random(unsigned int seed)
    : generator(seed), uniform_dist(0.0, 1.0) {}

// ---------------------------------------------------------------------------
//  uniforme() - Draw a uniform random number in [0, 1)
// ---------------------------------------------------------------------------
double Random::uniforme() {
    return uniform_dist(generator);
}

// ---------------------------------------------------------------------------
//  normale() - Draw a standard normal random number via Box-Muller transform
// ---------------------------------------------------------------------------
double Random::normale() {
    double u1 = uniforme();
    while (u1 <= 1e-15) u1 = uniforme();  // Avoid log(0)
    double u2 = uniforme();

    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
}

// ---------------------------------------------------------------------------
//  abramowitz_stegun() - Standard normal CDF approximation
//  Reference: Abramowitz & Stegun, Handbook of Mathematical Functions (1964)
//  Maximum absolute error: 7.5e-8
// ---------------------------------------------------------------------------
double Random::abramowitz_stegun(double x) {
    if (x < 0.0) {
        return 1.0 - abramowitz_stegun(-x);
    }

    // Rational approximation coefficients
    const double b0 =  0.2316419;
    const double b1 =  0.319381530;
    const double b2 = -0.356563782;
    const double b3 =  1.781477937;
    const double b4 = -1.8211255978;
    const double b5 =  1.330274429;

    double t = 1.0 / (1.0 + b0 * x);
    double poly = b1 * t + b2 * t * t + b3 * t * t * t
                + b4 * t * t * t * t + b5 * t * t * t * t * t;

    return 1.0 - std::exp(-x * x * 0.5) * poly / std::sqrt(2.0 * M_PI);
}
