// ============================================================================
//  PricerArbreBinomial.hpp - Cox-Ross-Rubinstein Binomial Tree Pricer
//  Prices American options via backward induction with early exercise
// ============================================================================

#ifndef PRICERARBREBINOMIAL_HPP
#define PRICERARBREBINOMIAL_HPP

#include "BlackScholes.hpp"

class PricerArbreBinomial {
private:
    int    N;    // Number of time steps in the tree
    double K;    // Strike price
    double T;    // Maturity

public:
    PricerArbreBinomial(int nb_pas, double strike, double maturite);

    // Price of an American Put option
    double prix_put_americain(const ParametresBS& p);

    // Price of an American Call option
    double prix_call_americain(const ParametresBS& p);
};

#endif // PRICERARBREBINOMIAL_HPP
