// ============================================================================
//  BermudeNDatesLS.hpp - Longstaff-Schwartz Method for N-Dates Bermudean Options
//  Uses least-squares regression to estimate continuation values
// ============================================================================

#ifndef BERMUDENDATESLS_HPP
#define BERMUDENDATESLS_HPP

#include <vector>
#include "BlackScholes.hpp"
#include "Random.hpp"

enum OptionType { PUT, CALL };

class BermudeNDatesLS {
private:
    std::vector<double> dates;  // Exercise dates
    double K;                   // Strike price
    OptionType type;            // PUT or CALL

    // Solve least-squares regression (quadratic polynomial) with data centering
    bool resoudre_moindres_carres(const std::vector<double>& X,
                                  const std::vector<double>& Y,
                                  std::vector<double>& coeff);

    // Evaluate a quadratic polynomial at x
    double evaluer_polynome(double x, const std::vector<double>& coeff);

    // Check if the option is in-the-money
    bool est_itm(double S) const;

    // Compute the option payoff
    double payoff(double S) const;

public:
    BermudeNDatesLS(const std::vector<double>& dates_exercice,
                    double strike, OptionType t = PUT);

    // Longstaff-Schwartz Monte Carlo pricing with 90% confidence interval
    void monte_carlo_longstaff_schwartz(const ParametresBS& p, long long M,
                                        Random& rng,
                                        double& prix_mc, double& ic_inf, double& ic_sup);
};

#endif // BERMUDENDATESLS_HPP
