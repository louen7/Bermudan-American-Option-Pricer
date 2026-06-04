// ============================================================================
//  PricerArbreBinomial.cpp - Cox-Ross-Rubinstein Binomial Tree Implementation
// ============================================================================

#include "PricerArbreBinomial.hpp"
#include <cmath>
#include <algorithm>
#include <vector>

// ---------------------------------------------------------------------------
//  Constructor
// ---------------------------------------------------------------------------
PricerArbreBinomial::PricerArbreBinomial(int nb_pas, double strike, double maturite)
    : N(nb_pas), K(strike), T(maturite) {}

// ---------------------------------------------------------------------------
//  prix_put_americain() - American Put via CRR binomial tree
//
//  Algorithm:
//    1. Build forward tree of stock prices using Cox-Ross-Rubinstein parameters
//    2. Initialize option values at maturity with the put payoff
//    3. Backward induction: at each node, take max(exercise, continuation)
// ---------------------------------------------------------------------------
double PricerArbreBinomial::prix_put_americain(const ParametresBS& p) {

    // CRR tree parameters
    double dt     = T / N;
    double u      = std::exp(p.sigma * std::sqrt(dt));
    double d      = 1.0 / u;
    double discount = std::exp(-p.r * dt);
    double q_prob = (std::exp((p.r - p.q) * dt) - d) / (u - d);

    // Allocate triangular grids for stock prices and option values
    std::vector<std::vector<double>> S(N + 1);
    std::vector<std::vector<double>> V(N + 1);
    for (int t = 0; t <= N; ++t) {
        S[t].resize(t + 1);
        V[t].resize(t + 1);
    }

    // Forward phase: build the stock price tree
    S[0][0] = p.S0;
    for (int t = 1; t <= N; ++t) {
        for (int i = 0; i <= t; ++i) {
            S[t][i] = p.S0 * std::pow(u, i) * std::pow(d, t - i);
        }
    }

    // Initialize option values at maturity
    for (int i = 0; i <= N; ++i) {
        V[N][i] = std::max(K - S[N][i], 0.0);
    }

    // Backward induction with American early exercise condition
    for (int t = N - 1; t >= 0; --t) {
        for (int i = 0; i <= t; ++i) {
            double continuation = discount * (q_prob * V[t + 1][i + 1]
                                            + (1.0 - q_prob) * V[t + 1][i]);
            double exercise     = std::max(K - S[t][i], 0.0);
            V[t][i] = std::max(exercise, continuation);
        }
    }

    return V[0][0];
}

// ---------------------------------------------------------------------------
//  prix_call_americain() - American Call via CRR binomial tree
// ---------------------------------------------------------------------------
double PricerArbreBinomial::prix_call_americain(const ParametresBS& p) {

    double dt       = T / N;
    double u        = std::exp(p.sigma * std::sqrt(dt));
    double d        = 1.0 / u;
    double discount = std::exp(-p.r * dt);
    double q_prob   = (std::exp((p.r - p.q) * dt) - d) / (u - d);

    std::vector<std::vector<double>> S(N + 1);
    std::vector<std::vector<double>> V(N + 1);
    for (int t = 0; t <= N; ++t) {
        S[t].resize(t + 1);
        V[t].resize(t + 1);
    }

    S[0][0] = p.S0;
    for (int t = 1; t <= N; ++t) {
        for (int i = 0; i <= t; ++i) {
            S[t][i] = p.S0 * std::pow(u, i) * std::pow(d, t - i);
        }
    }

    for (int i = 0; i <= N; ++i) {
        V[N][i] = std::max(S[N][i] - K, 0.0);
    }

    for (int t = N - 1; t >= 0; --t) {
        for (int i = 0; i <= t; ++i) {
            double continuation = discount * (q_prob * V[t + 1][i + 1]
                                            + (1.0 - q_prob) * V[t + 1][i]);
            double exercise     = std::max(S[t][i] - K, 0.0);
            V[t][i] = std::max(exercise, continuation);
        }
    }

    return V[0][0];
}
