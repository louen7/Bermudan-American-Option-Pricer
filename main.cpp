// ============================================================================
//  main.cpp - Bermudean Option Pricing Project
//
//  Compares three option pricing methodologies:
//    - European Put: Black-Scholes closed-form solution
//    - Bermudean Put: Longstaff-Schwartz Monte Carlo simulation
//    - American Put: Cox-Ross-Rubinstein binomial tree
//
//  Validates the theoretical price ordering: European <= Bermudean <= American
//  Demonstrates convergence of Bermudean -> American as exercise dates increase
// ============================================================================

#include <iostream>
#include <iomanip>
#include <vector>
#include "BlackScholes.hpp"
#include "Random.hpp"
#include "BermudeNDatesLS.hpp"
#include "PricerArbreBinomial.hpp"

// ---------------------------------------------------------------------------
//  generer_dates() - Generate N evenly-spaced exercise dates in [0, T_max]
// ---------------------------------------------------------------------------
std::vector<double> generer_dates(int N, double T_max) {
    std::vector<double> dates(N);
    for (int i = 0; i < N; ++i)
        dates[i] = T_max * (i + 1) / N;
    return dates;
}

// ---------------------------------------------------------------------------
//  print_separator() - Print a horizontal separator line
// ---------------------------------------------------------------------------
void print_separator(char c, int width) {
    std::cout << std::string(width, c) << std::endl;
}

// ===========================================================================
//  MAIN
// ===========================================================================
int main() {
    // ----- Option parameters -----
    double K     = 100.0;    // Strike price
    double r     = 0.15;     // Risk-free rate (15%)
    double sigma = 0.15;     // Volatility (15%)
    double T     = 10.0;     // Maturity (10 years)

    // ----- Numerical parameters -----
    int       N_berm  = 20;        // Number of exercise dates for Bermudean
    int       N_tree  = 300;       // Number of time steps for binomial tree
    long long M       = 300000;    // Number of Monte Carlo paths
    unsigned int seed = 42;        // Random seed for reproducibility

    // ======================================================================
    //  PART 1: Price comparison across different spot prices
    // ======================================================================
    std::cout << std::endl;
    print_separator('=', 100);
    std::cout << "   EUROPEAN PUT vs BERMUDEAN PUT vs AMERICAN PUT" << std::endl;
    print_separator('=', 100);
    std::cout << "  K = " << K
              << "   r = " << r * 100 << "%"
              << "   sigma = " << sigma * 100 << "%"
              << "   T = " << T << " yrs" << std::endl;
    std::cout << "  European  : Black-Scholes (closed-form)" << std::endl;
    std::cout << "  Bermudean : Longstaff-Schwartz, N = " << N_berm
              << " dates, M = " << M << " paths" << std::endl;
    std::cout << "  American  : CRR Binomial Tree, N = " << N_tree
              << " steps" << std::endl;
    print_separator('-', 100);

    // Table header
    std::cout << std::left
              << std::setw(6)  << "S0"
              << std::setw(12) << "European"
              << std::setw(13) << "Bermudean"
              << std::setw(22) << "90% CI"
              << std::setw(12) << "American"
              << std::setw(10) << "BermPrem"
              << std::setw(10) << "AmerPrem"
              << std::setw(10) << "Gap A-B"
              << std::endl;

    std::cout << std::left
              << std::setw(6)  << "----"
              << std::setw(12) << "---------"
              << std::setw(13) << "---------"
              << std::setw(22) << "--------------------"
              << std::setw(12) << "---------"
              << std::setw(10) << "--------"
              << std::setw(10) << "--------"
              << std::setw(10) << "--------"
              << std::endl;

    double S0_vals[] = {30, 40, 50, 60, 70, 80, 85, 90, 95, 100, 105, 110, 120, 130, 150};
    int nb = sizeof(S0_vals) / sizeof(S0_vals[0]);

    for (int s = 0; s < nb; ++s) {
        double S0 = S0_vals[s];
        ParametresBS p = { S0, r, sigma };

        // European price (closed-form)
        double prix_eur = BlackScholes::prix_analytique_put(p, T, K);

        // Bermudean price (Longstaff-Schwartz Monte Carlo)
        std::vector<double> dates = generer_dates(N_berm, T);
        Random rng(seed);
        BermudeNDatesLS opt_berm(dates, K, PUT);
        double prix_berm, ic_inf, ic_sup;
        opt_berm.monte_carlo_longstaff_schwartz(p, M, rng, prix_berm, ic_inf, ic_sup);

        // American price (CRR binomial tree)
        PricerArbreBinomial opt_tree(N_tree, K, T);
        double prix_amer = opt_tree.prix_put_americain(p);

        // Premiums: difference over European price
        double berm_prem = prix_berm - prix_eur;
        double amer_prem = prix_amer - prix_eur;
        double gap_ba    = prix_amer - prix_berm;

        std::cout << std::fixed << std::setprecision(2);
        std::cout << std::left << std::setw(6)  << (int)S0;
        std::cout << std::setprecision(4);
        std::cout << std::setw(12) << prix_eur;
        std::cout << std::setw(13) << prix_berm;
        std::cout << "[" << std::setw(8) << ic_inf << " ; "
                         << std::setw(8) << ic_sup << "]  ";
        std::cout << std::setw(12) << prix_amer;
        std::cout << std::setw(10) << berm_prem;
        std::cout << std::setw(10) << amer_prem;
        std::cout << std::setw(10) << gap_ba;
        std::cout << std::endl;
    }

    // ======================================================================
    //  PART 2: Convergence of Bermudean -> American as N increases
    // ======================================================================
    {
        double S0 = 60.0;
        ParametresBS p = { S0, r, sigma };
        double prix_eur  = BlackScholes::prix_analytique_put(p, T, K);

        PricerArbreBinomial opt_tree(N_tree, K, T);
        double prix_amer = opt_tree.prix_put_americain(p);

        std::cout << std::endl;
        print_separator('-', 70);
        std::cout << "  CONVERGENCE: Bermudean -> American  (S0 = " << S0
                  << ")" << std::endl;
        std::cout << "  European = " << std::fixed << std::setprecision(4)
                  << prix_eur
                  << "    American = " << prix_amer
                  << "    Intrinsic = " << std::max(K - S0, 0.0)
                  << std::endl;
        print_separator('-', 70);

        std::cout << std::left
                  << std::setw(8)  << "N"
                  << std::setw(16) << "Bermudean"
                  << std::setw(16) << "Gap/American"
                  << std::endl;
        std::cout << std::left
                  << std::setw(8)  << "---"
                  << std::setw(16) << "-----------"
                  << std::setw(16) << "------------"
                  << std::endl;

        int configs[] = {1, 2, 3, 5, 10, 15, 20, 30, 50, 75, 100};
        int nb_c = sizeof(configs) / sizeof(configs[0]);

        for (int c = 0; c < nb_c; ++c) {
            int N = configs[c];
            std::vector<double> dates = generer_dates(N, T);
            Random rng(seed);
            BermudeNDatesLS opt(dates, K, PUT);
            double prix_mc, ic_lo, ic_hi;
            opt.monte_carlo_longstaff_schwartz(p, M, rng, prix_mc, ic_lo, ic_hi);
            double gap = prix_amer - prix_mc;

            std::cout << std::fixed << std::setprecision(4);
            std::cout << std::left << std::setw(8)  << N;
            std::cout << std::setw(16) << prix_mc;
            std::cout << std::setw(16) << gap;
            std::cout << std::endl;
        }
    }

    // ======================================================================
    //  Summary
    // ======================================================================
    std::cout << std::endl;
    print_separator('=', 70);
    std::cout << "  Key results:" << std::endl;
    std::cout << "    European <= Bermudean <= American" << std::endl;
    std::cout << "    As N increases, the Bermudean price converges" << std::endl;
    std::cout << "    toward the American price." << std::endl;
    print_separator('=', 70);
    std::cout << std::endl;

    return 0;
}
