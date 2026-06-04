// ============================================================================
//  Bermudes2Dates.hpp - 2-Dates Bermudean Option (Analytical + Monte Carlo)
//  Computes the critical exercise boundary and prices via MC with antithetics
// ============================================================================

#ifndef BERMUDE2DATES_HPP
#define BERMUDE2DATES_HPP

#include "BlackScholes.hpp"
#include "Random.hpp"

class Bermude2Dates {
private:
    double T1;  // First exercise date
    double T2;  // Second exercise date (final maturity)
    double K;   // Strike price

    // Equation for the critical boundary: Payoff(T1) - ContinuationValue(T1) = 0
    // The unique root of this equation defines the critical boundary K_bar
    double equation_K_bar(double S_T1, const ParametresBS& p);

public:
    Bermude2Dates(double t1, double t2, double strike);

    // Find the critical exercise boundary K_bar via bisection
    double calculer_K_bar(const ParametresBS& p);

    // Monte Carlo pricing with optional antithetic variance reduction
    void monte_carlo_bermude(const ParametresBS& p, long long M, Random& rng,
                             bool utiliser_antithetique,
                             double& prix_mc, double& ic_inf, double& ic_sup);
};

#endif // BERMUDE2DATES_HPP
