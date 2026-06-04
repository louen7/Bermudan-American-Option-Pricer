#include "BermudeNDatesLS.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

// Constructeur
BermudeNDatesLS::BermudeNDatesLS(const std::vector<double>& dates_exercice, double strike, OptionType t)
    : dates(dates_exercice), K(strike), type(t) {}

// Helpers : condition ITM et payoff selon le type d'option
bool BermudeNDatesLS::est_itm(double S) const {
    return (type == PUT) ? (S < K) : (S > K);
}

double BermudeNDatesLS::payoff(double S) const {
    return (type == PUT) ? std::max(K - S, 0.0) : std::max(S - K, 0.0);
}

// FONCTION 1 : Moindres Carrés avec centrage des données pour stabilité numérique
bool BermudeNDatesLS::resoudre_moindres_carres(const std::vector<double>& X,
                                              const std::vector<double>& Y,
                                              std::vector<double>& coeff) {
    size_t n = X.size();
    coeff.assign(3, 0.0);
    if (n < 3) return false;

    // Centrage de X autour de sa moyenne pour eviter l'annulation catastrophique
    // dans le determinant quand les X sont loin de 0
    double mean_x = 0.0;
    for (size_t i = 0; i < n; ++i) mean_x += X[i];
    mean_x /= static_cast<double>(n);

    double s_x0 = static_cast<double>(n);
    double s_xc1 = 0.0, s_xc2 = 0.0, s_xc3 = 0.0, s_xc4 = 0.0;
    double s_y  = 0.0, s_xcy = 0.0, s_xc2y = 0.0;

    for (size_t i = 0; i < n; ++i) {
        double xc = X[i] - mean_x;
        double xc2 = xc * xc;
        s_xc1 += xc; s_xc2 += xc2; s_xc3 += xc2 * xc; s_xc4 += xc2 * xc2;
        s_y   += Y[i]; s_xcy += xc * Y[i]; s_xc2y += xc2 * Y[i];
    }

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

    // Coefficients du polynome centree : Y = ac*(X-m)² + bc*(X-m) + cc
    double ac = detA / detM;
    double bc = detB / detM;
    double cc = detC / detM;

    // Transformation inverse : Y = a*X² + b*X + c
    coeff[2] = ac;
    coeff[1] = bc - 2.0 * ac * mean_x;
    coeff[0] = ac * mean_x * mean_x - bc * mean_x + cc;

    return true;
}

// FONCTION 2 : Évaluation du polynôme
double BermudeNDatesLS::evaluer_polynome(double x, const std::vector<double>& coeff) {
    return coeff[2] * x * x + coeff[1] * x + coeff[0];
}

// FONCTION 3 : Le moteur général à N dates
void BermudeNDatesLS::monte_carlo_longstaff_schwartz(const ParametresBS& p, long long M, Random& rng,
                                                    double& prix_mc, double& ic_inf, double& ic_sup) {
    size_t N = dates.size(); // Nombre de dates d'exercice
    
    // Grille de prix : Matrice de taille N lignes x M colonnes
    std::vector<std::vector<double>> S(N, std::vector<double>(M));
    
    // Vecteurs pour stocker les cash-flows optimaux et leur date d'exercice associée
    std::vector<double> V(M);
    std::vector<double> date_exercice(M);

    // --- PHASE FORWARD : Génération de toute la grille de prix de 0 à T_N ---
    for (long long i = 0; i < M; ++i) {
        // Première étape : de t=0 à dates[0] (T1)
        double Z = rng.normale();
        S[0][i] = BlackScholes::simuler_ST(p, dates[0], Z);
        
        // Étapes suivantes : de dates[t-1] à dates[t]
        for (size_t t = 1; t < N; ++t) {
            Z = rng.normale();
            ParametresBS p_prec = { S[t-1][i], p.r, p.sigma, p.q };
            S[t][i] = BlackScholes::simuler_ST(p_prec, dates[t] - dates[t-1], Z);
        }

        // Initialisation à la date finale d'échéance T_N (indice N-1)
        V[i] = payoff(S[N-1][i]);
        date_exercice[i] = dates[N-1];
    }

    // --- PHASE BACKWARD : Grande boucle temporelle inversée (de T_{N-1} à T_1) ---
    // En indices informatiques : on part de l'avant-dernière ligne (N-2) jusqu'à la première (0)
    for (long long t = static_cast<long long>(N) - 2; t >= 0; --t) {
        
        std::vector<double> X;
        std::vector<double> Y;
        std::vector<long long> indices_itm;

        double t_actuel = dates[t];

        // 1. Filtrage des trajectoires In-The-Money (ITM) à la date courante t
        for (long long i = 0; i < M; ++i) {
            if (est_itm(S[t][i])) {
                X.push_back(S[t][i]); // X = Prix actuel de l'actif
                
                // Y = Gain futur optimal de cette trajectoire actualisé jusqu'à t_actuel
                double dt = date_exercice[i] - t_actuel;
                Y.push_back(V[i] * std::exp(-p.r * dt));
                indices_itm.push_back(i);
            }
        }

        // 2. Calcul de la régression pour la date courante t
        std::vector<double> coeff;
        if (resoudre_moindres_carres(X, Y, coeff)) {
            
            // 3. Application de la Règle d'or sur les trajectoires concernées
            for (size_t k = 0; k < indices_itm.size(); ++k) {
                long long i = indices_itm[k];
                double payoff_immediat = payoff(S[t][i]);
                double continuation_estimee = evaluer_polynome(S[t][i], coeff);

                // Si l'exercice immédiat rapporte plus que la continuation moyenne estimée
                if (payoff_immediat > continuation_estimee) {
                    V[i] = payoff_immediat;       // Remplacement par le gain immédiat
                    date_exercice[i] = t_actuel;   // Mise à jour de la date d'arrêt
                }
            }
        }
    }

    // --- PHASE CALCUL DU PRIX : Actualisation finale globale à t = 0 ---
    double somme_payoffs = 0.0;
    double somme_payoffs_carre = 0.0;

    for (long long i = 0; i < M; ++i) {
        double gain_t0 = V[i] * std::exp(-p.r * date_exercice[i]);
        somme_payoffs += gain_t0;
        somme_payoffs_carre += gain_t0 * gain_t0;
    }

    prix_mc = somme_payoffs / M;
    double variance = (somme_payoffs_carre / M) - (prix_mc * prix_mc);
    double ecart_type = std::sqrt(variance / M);

    double z_90 = 1.645; 
    ic_inf = prix_mc - z_90 * ecart_type;
    ic_sup = prix_mc + z_90 * ecart_type;
}