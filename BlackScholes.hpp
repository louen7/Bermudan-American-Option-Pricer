// ============================================================================
//  BlackScholes.hpp - Black-Scholes Model for European Option Pricing
//  Provides: analytical pricing, asset simulation, and Monte Carlo estimation
// ============================================================================

#ifndef BLACKSCHOLES_HPP
#define BLACKSCHOLES_HPP

#include "Random.hpp"

// Black-Scholes model parameters
struct ParametresBS {
    double S0;     // Spot price
    double r;      // Risk-free interest rate
    double sigma;  // Volatility
    double q;      // Continuous dividend yield (0 if no dividends)
};

class BlackScholes {
public:
    // Simulate asset price S(T) from a standard normal draw Z
    static double simuler_ST(const ParametresBS& p, double T, double Z);

    // Analytical price of a European Put
    static double prix_analytique_put(const ParametresBS& p, double T, double K);

    // Analytical price of a European Call (via put-call parity)
    static double prix_analytique_call(const ParametresBS& p, double T, double K);

    // Monte Carlo price of a European Put with 90% confidence interval
    static void monte_carlo_europeen(const ParametresBS& p, double T, double K,
                                     long long M, Random& rng,
                                     double& prix_mc, double& ic_inf, double& ic_sup);
};

#endif // BLACKSCHOLES_HPP
