// ============================================================================
//  BermudeNDatesLS.cpp - Longstaff-Schwartz N-Dates Bermudean Option Pricing
//
//  Reference: Longstaff, F.A. & Schwartz, E.S. (2001)
//  "Valuing American Options by Simulation: A Simple Least-Squares Approach"
//  The Review of Financial Studies, Vol. 14, No. 1, pp. 113-147
// ============================================================================

#include "BermudeNDatesLS.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

// ---------------------------------------------------------------------------
//  Constructor
// ---------------------------------------------------------------------------
BermudeNDatesLS::BermudeNDatesLS(const std::vector<double>& dates_exercice,
                                 double strike, OptionType t)
    : dates(dates_exercice), K(strike), type(t) {}

// ---------------------------------------------------------------------------
//  est_itm() - Check if the option is in-the-money for a given spot price
// ---------------------------------------------------------------------------
bool BermudeNDatesLS::est_itm(double S) const {
    return (type == PUT) ? (S < K) : (S > K);
}

// ---------------------------------------------------------------------------
//  payoff() - Compute the option payoff for a given spot price
// ---------------------------------------------------------------------------
double BermudeNDatesLS::payoff(double S) const {
    return (type == PUT) ? std::max(K - S, 0.0) : std::max(S - K, 0.0);
}

// ---------------------------------------------------------------------------
//  resoudre_moindres_carres() - Least-squares regression for a quadratic
//  polynomial Y = a*X^2 + b*X + c, using centered data for numerical stability
// ---------------------------------------------------------------------------
bool BermudeNDatesLS::resoudre_moindres_carres(const std::vector<double>& X,
                                               const std::vector<double>& Y,
                                               std::vector<double>& coeff) {
    size_t n = X.size();
    coeff.assign(3, 0.0);
    if (n < 3) return false;

    // Center X around its mean to avoid catastrophic cancellation in det
    double mean_x = 0.0;
    for (size_t i = 0; i < n; ++i) mean_x += X[i];
    mean_x /= static_cast<double>(n);

    double s_x0 = static_cast<double>(n);
    double s_xc1 = 0.0, s_xc2 = 0.0, s_xc3 = 0.0, s_xc4 = 0.0;
    double s_y  = 0.0, s_xcy = 0.0, s_xc2y = 0.0;

    for (size_t i = 0; i < n; ++i) {
        double xc  = X[i] - mean_x;
        double xc2 = xc * xc;
        s_xc1 += xc;   s_xc2 += xc2;
        s_xc3 += xc2 * xc;  s_xc4 += xc2 * xc2;
        s_y += Y[i];    s_xcy += xc * Y[i];  s_xc2y += xc2 * Y[i];
    }

    // Determinant of the normal equation matrix (centered)
    double detM = s_x0 * (s_xc2 * s_xc4 - s_xc3 * s_xc3)
                - s_xc1 * (s_xc1 * s_xc4 - s_xc3 * s_xc2)
                + s_xc2 * (s_xc1 * s_xc3 - s_xc2 * s_xc2);

    if (std::abs(detM) < 1e-12) return false;

    double detC = s_y   * (s_xc2 * s_xc4 - s_xc3 * s_xc3)
                - s_xc1 * (s_xcy * s_xc4 - s_xc3 * s_xc2y)
                + s_xc2 * (s_xcy * s_xc3 - s_xc2 * s_xc2y);

    double detB = s_x0  * (s_xcy * s_xc4 - s_xc2y * s_xc3)
                - s_y   * (s_xc1 * s_xc4 - s_xc2 * s_xc3)
                + s_xc2 * (s_xc1 * s_xc2y - s_xc2 * s_xcy);

    double detA = s_x0  * (s_xc2 * s_xc2y - s_xc3 * s_xcy)
                - s_xc1 * (s_xc1 * s_xc2y - s_xc3 * s_y)
                + s_y   * (s_xc1 * s_xc3 - s_xc2 * s_xc2);

    // Centered polynomial: Y = ac*(X-m)^2 + bc*(X-m) + cc
    double ac = detA / detM;
    double bc = detB / detM;
    double cc = detC / detM;

    // Convert back to standard form: Y = a*X^2 + b*X + c
    coeff[2] = ac;
    coeff[1] = bc - 2.0 * ac * mean_x;
    coeff[0] = ac * mean_x * mean_x - bc * mean_x + cc;

    return true;
}

// ---------------------------------------------------------------------------
//  evaluer_polynome() - Evaluate the quadratic polynomial at x
// ---------------------------------------------------------------------------
double BermudeNDatesLS::evaluer_polynome(double x, const std::vector<double>& coeff) {
    return coeff[2] * x * x + coeff[1] * x + coeff[0];
}

// ---------------------------------------------------------------------------
//  monte_carlo_longstaff_schwartz() - N-dates Bermudean pricing engine
//
//  Algorithm:
//    1. FORWARD phase: simulate M stock price paths across all exercise dates
//    2. BACKWARD phase: at each exercise date, regress continuation values
//       on a quadratic polynomial of spot prices (ITM paths only)
//    3. Decide early exercise vs. continuation at each date
//    4. Compute the discounted average payoff with 90% confidence interval
// ---------------------------------------------------------------------------
void BermudeNDatesLS::monte_carlo_longstaff_schwartz(const ParametresBS& p,
                                                     long long M, Random& rng,
                                                     double& prix_mc,
                                                     double& ic_inf,
                                                     double& ic_sup) {
    size_t N = dates.size();

    // Price grid: N rows (dates) x M columns (paths)
    std::vector<std::vector<double> > S(N, std::vector<double>(M));

    // Optimal cash-flows and their associated exercise dates
    std::vector<double> V(M);
    std::vector<double> exercise_date(M);

    // --- FORWARD PHASE: generate all stock prices from t=0 to T_N ---
    for (long long i = 0; i < M; ++i) {
        // First step: from t=0 to dates[0]
        double Z = rng.normale();
        S[0][i] = BlackScholes::simuler_ST(p, dates[0], Z);

        // Subsequent steps: from dates[t-1] to dates[t]
        for (size_t t = 1; t < N; ++t) {
            Z = rng.normale();
            ParametresBS p_prev = { S[t - 1][i], p.r, p.sigma, p.q };
            S[t][i] = BlackScholes::simuler_ST(p_prev, dates[t] - dates[t - 1], Z);
        }

        // Initialize at the final exercise date T_N
        V[i] = payoff(S[N - 1][i]);
        exercise_date[i] = dates[N - 1];
    }

    // --- BACKWARD PHASE: iterate from T_{N-1} back to T_1 ---
    for (long long t = static_cast<long long>(N) - 2; t >= 0; --t) {

        std::vector<double> X;
        std::vector<double> Y;
        std::vector<long long> itm_indices;

        double current_date = dates[t];

        // Collect in-the-money paths
        for (long long i = 0; i < M; ++i) {
            if (est_itm(S[t][i])) {
                X.push_back(S[t][i]);  // Spot price at current date

                // Discounted future optimal cash-flow
                double dt = exercise_date[i] - current_date;
                Y.push_back(V[i] * std::exp(-p.r * dt));
                itm_indices.push_back(i);
            }
        }

        // Perform least-squares regression
        std::vector<double> coeff;
        if (resoudre_moindres_carres(X, Y, coeff)) {
            // Apply the optimal exercise rule on ITM paths
            for (size_t k = 0; k < itm_indices.size(); ++k) {
                long long i = itm_indices[k];
                double immediate_payoff = payoff(S[t][i]);
                double estimated_continuation = evaluer_polynome(S[t][i], coeff);

                // Exercise immediately if payoff exceeds estimated continuation
                if (immediate_payoff > estimated_continuation) {
                    V[i] = immediate_payoff;
                    exercise_date[i] = current_date;
                }
            }
        }
    }

    // --- PRICING PHASE: discount all cash-flows to t=0 ---
    double sum_payoff    = 0.0;
    double sum_payoff_sq = 0.0;

    for (long long i = 0; i < M; ++i) {
        double gain_t0 = V[i] * std::exp(-p.r * exercise_date[i]);
        sum_payoff    += gain_t0;
        sum_payoff_sq += gain_t0 * gain_t0;
    }

    prix_mc = sum_payoff / M;

    double variance  = (sum_payoff_sq / M) - (prix_mc * prix_mc);
    double std_error = std::sqrt(variance / M);
    double z_90      = 1.645;

    ic_inf = prix_mc - z_90 * std_error;
    ic_sup = prix_mc + z_90 * std_error;
}
