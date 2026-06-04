#ifndef BERMUDE3DATESLS_HPP
#define BERMUDE3DATESLS_HPP

#include <vector>
#include "BlackScholes.hpp"
#include "Random.hpp"

class Bermude3DatesLS {
private:
    double T1;
    double T2;
    double T3;
    double K;

    bool resoudre_moindres_carres(const std::vector<double>& X,
                                  const std::vector<double>& Y,
                                  std::vector<double>& coeff);

    double evaluer_polynome(double x, const std::vector<double>& coeff);

public:
    Bermude3DatesLS(double t1, double t2, double t3, double strike);

    void monte_carlo_longstaff_schwartz(const ParametresBS& p, long long M, Random& rng,
                                        double& prix_mc, double& ic_inf, double& ic_sup);
};

#endif // BERMUDE3DATESLS_HPP