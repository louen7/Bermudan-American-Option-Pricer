// ============================================================================
//  BlackScholes.cpp - Black-Scholes Model Implementation
// ============================================================================

#include "BlackScholes.hpp"
#include <cmath>
#include <algorithm>

// ---------------------------------------------------------------------------
//  simuler_ST() - Simulate the asset price at maturity T
//  S(T) = S0 * exp((r - q - sigma^2/2)*T + sigma*Z*sqrt(T))
// ---------------------------------------------------------------------------
double BlackScholes::simuler_ST(const ParametresBS& p, double T, double Z) {
    return p.S0 * std::exp((p.r - p.q - p.sigma * p.sigma / 2.0) * T
                           + p.sigma * Z * std::sqrt(T));
}

// ---------------------------------------------------------------------------
//  prix_analytique_put() - Closed-form European Put price
//  P = K*exp(-rT)*N(-d2) - S0*exp(-qT)*N(-d1)
// ---------------------------------------------------------------------------
double BlackScholes::prix_analytique_put(const ParametresBS& p, double T, double K) {
    double d1 = (std::log(p.S0 / K) + (p.r - p.q + p.sigma * p.sigma * 0.5) * T)
                / (p.sigma * std::sqrt(T));
    double d2 = d1 - p.sigma * std::sqrt(T);

    double N_minus_d1 = Random::abramowitz_stegun(-d1);
    double N_minus_d2 = Random::abramowitz_stegun(-d2);

    return K * std::exp(-p.r * T) * N_minus_d2
         - p.S0 * std::exp(-p.q * T) * N_minus_d1;
}

// ---------------------------------------------------------------------------
//  prix_analytique_call() - Closed-form European Call price
//  C = S0*exp(-qT)*N(d1) - K*exp(-rT)*N(d2)
// ---------------------------------------------------------------------------
double BlackScholes::prix_analytique_call(const ParametresBS& p, double T, double K) {
    double d1 = (std::log(p.S0 / K) + (p.r - p.q + p.sigma * p.sigma * 0.5) * T)
                / (p.sigma * std::sqrt(T));
    double d2 = d1 - p.sigma * std::sqrt(T);

    double N_plus_d1 = Random::abramowitz_stegun(d1);
    double N_plus_d2 = Random::abramowitz_stegun(d2);

    return p.S0 * std::exp(-p.q * T) * N_plus_d1
         - K * std::exp(-p.r * T) * N_plus_d2;
}

// ---------------------------------------------------------------------------
//  monte_carlo_europeen() - Monte Carlo pricing for a European Put
//  Returns: price estimate and 90% confidence interval [ic_inf, ic_sup]
// ---------------------------------------------------------------------------
void BlackScholes::monte_carlo_europeen(const ParametresBS& p, double T, double K,
                                        long long M, Random& rng,
                                        double& prix_mc, double& ic_inf, double& ic_sup) {
    double sum_payoff     = 0.0;
    double sum_payoff_sq  = 0.0;
    double discount       = std::exp(-p.r * T);

    for (long long i = 0; i < M; ++i) {
        double Z      = rng.normale();
        double ST     = simuler_ST(p, T, Z);
        double payoff = std::max(K - ST, 0.0) * discount;

        sum_payoff    += payoff;
        sum_payoff_sq += payoff * payoff;
    }

    prix_mc = sum_payoff / M;

    // Variance and standard error of the estimator
    double variance   = (sum_payoff_sq / M) - (prix_mc * prix_mc);
    double std_error  = std::sqrt(variance / M);
    double z_90       = 1.645;

    ic_inf = prix_mc - z_90 * std_error;
    ic_sup = prix_mc + z_90 * std_error;
}
