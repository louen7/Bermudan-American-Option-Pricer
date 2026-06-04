// ============================================================================
//  Bermudes2Dates.cpp - 2-Dates Bermudean Option Implementation
// ============================================================================

#include "Bermudes2Dates.hpp"
#include <cmath>
#include <algorithm>

// ---------------------------------------------------------------------------
//  Constructor
// ---------------------------------------------------------------------------
Bermude2Dates::Bermude2Dates(double t1, double t2, double strike)
    : T1(t1), T2(t2), K(strike) {}

// ---------------------------------------------------------------------------
//  equation_K_bar() - Arbitrage condition at T1
//  Returns: immediate_payoff(T1) - continuation_value(T1)
//  The unique root of f(S) = 0 defines the critical boundary K_bar
// ---------------------------------------------------------------------------
double Bermude2Dates::equation_K_bar(double S_T1, const ParametresBS& p) {
    double immediate_payoff = std::max(K - S_T1, 0.0);
    ParametresBS p_T1 = { S_T1, p.r, p.sigma };
    double continuation_value = BlackScholes::prix_analytique_put(p_T1, T2 - T1, K);
    return immediate_payoff - continuation_value;
}

// ---------------------------------------------------------------------------
//  calculer_K_bar() - Find the critical boundary by bisection
//  The boundary separates the exercise region (S < K_bar) from the
//  continuation region (S >= K_bar)
// ---------------------------------------------------------------------------
double Bermude2Dates::calculer_K_bar(const ParametresBS& p) {
    double lo  = 0.001;
    double hi  = 2.0 * K;
    double eps = 1e-6;

    while (hi - lo > eps) {
        double mid = (lo + hi) / 2.0;
        if (equation_K_bar(mid, p) > 0.0) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return (lo + hi) / 2.0;
}

// ---------------------------------------------------------------------------
//  monte_carlo_bermude() - 2-dates Bermudean Monte Carlo pricing
//
//  Optional antithetic variance reduction:
//    When enabled, each iteration generates a pair of correlated paths
//    using Z and -Z for the first period, reducing the total variance
// ---------------------------------------------------------------------------
void Bermude2Dates::monte_carlo_bermude(const ParametresBS& p, long long M,
                                        Random& rng, bool utiliser_antithetique,
                                        double& prix_mc, double& ic_inf,
                                        double& ic_sup) {

    // Step 1: compute the critical boundary K_bar
    double K_bar = calculer_K_bar(p);

    double sum_payoff    = 0.0;
    double sum_payoff_sq = 0.0;

    // With antithetic variates, each iteration produces a pair (Z, -Z)
    // so we halve the loop count to maintain M total paths
    long long simulations = utiliser_antithetique ? M / 2 : M;

    // Discount factors for cash-flows at T1 and T2
    double disc_T1 = std::exp(-p.r * T1);
    double disc_T2 = std::exp(-p.r * T2);

    // Step 2: simulation loop
    for (long long i = 0; i < simulations; ++i) {

        // --- PATH 1: standard or first of the antithetic pair ---
        double Z1     = rng.normale();
        double S_T1_1 = BlackScholes::simuler_ST(p, T1, Z1);
        double payoff_1 = 0.0;

        if (S_T1_1 < K_bar) {
            // Below critical boundary: exercise immediately at T1
            payoff_1 = std::max(K - S_T1_1, 0.0) * disc_T1;
        } else {
            // Above critical boundary: continue to T2
            double Z2     = rng.normale();
            ParametresBS p_T1 = { S_T1_1, p.r, p.sigma };
            double S_T2_1 = BlackScholes::simuler_ST(p_T1, T2 - T1, Z2);
            payoff_1 = std::max(K - S_T2_1, 0.0) * disc_T2;
        }

        // --- PATH 2: antithetic mirror (only if enabled) ---
        if (utiliser_antithetique) {
            double S_T1_2 = BlackScholes::simuler_ST(p, T1, -Z1);
            double payoff_2 = 0.0;

            if (S_T1_2 < K_bar) {
                payoff_2 = std::max(K - S_T1_2, 0.0) * disc_T1;
            } else {
                double Z2_anti = rng.normale();
                ParametresBS p_T1_2 = { S_T1_2, p.r, p.sigma };
                double S_T2_2 = BlackScholes::simuler_ST(p_T1_2, T2 - T1, Z2_anti);
                payoff_2 = std::max(K - S_T2_2, 0.0) * disc_T2;
            }

            // Average the antithetic pair
            double pair_avg = (payoff_1 + payoff_2) / 2.0;
            sum_payoff    += pair_avg;
            sum_payoff_sq += pair_avg * pair_avg;

        } else {
            sum_payoff    += payoff_1;
            sum_payoff_sq += payoff_1 * payoff_1;
        }
    }

    // Step 3: compute price estimate and 90% confidence interval
    prix_mc = sum_payoff / simulations;

    double variance  = (sum_payoff_sq / simulations) - (prix_mc * prix_mc);
    double std_error = std::sqrt(variance / simulations);
    double z_90      = 1.645;

    ic_inf = prix_mc - z_90 * std_error;
    ic_sup = prix_mc + z_90 * std_error;
}
