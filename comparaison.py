import numpy as np
import matplotlib.pyplot as plt
from math import log, sqrt, exp, erf
import subprocess

# ========================================================================
# Parametres avec dividendes pour rendre les calls interessants
# ========================================================================
r = 0.06
sigma = 0.25
K = 100.0
T = 5.0
q = 0.0        # Pas de dividende
N_berm = 20
N_arbre = 300
M = 300000
graine = 42

S0_range = np.linspace(50, 180, 27)

def norm_cdf(x):
    return 0.5 * (1.0 + erf(x / sqrt(2.0)))

def bs_call(S0):
    d1 = (log(S0 / K) + (r - q + 0.5 * sigma**2) * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    return S0 * exp(-q * T) * norm_cdf(d1) - K * exp(-r * T) * norm_cdf(d2)

def bs_put(S0):
    d1 = (log(S0 / K) + (r - q + 0.5 * sigma**2) * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    return K * exp(-r * T) * norm_cdf(-d2) - S0 * exp(-q * T) * norm_cdf(-d1)

prix_eur_call = np.array([bs_call(S0) for S0 in S0_range])
prix_eur_put = np.array([bs_put(S0) for S0 in S0_range])
intrinseque_call = np.array([max(S0 - K, 0.0) for S0 in S0_range])
intrinseque_put = np.array([max(K - S0, 0.0) for S0 in S0_range])

# ========================================================================
# Compilation
# ========================================================================
tmp_main = r"""
#include <iostream>
#include <vector>
#include "BlackScholes.hpp"
#include "Random.hpp"
#include "BermudeNDatesLS.hpp"
#include "PricerArbreBinomial.hpp"

int main() {
    double S0, r_val, sigma_val, K_val, T_max, q_val;
    int N_berm, N_arbre, type_int;
    long long M;
    unsigned int graine;
    std::cin >> S0 >> r_val >> sigma_val >> K_val >> T_max >> q_val >> N_berm >> N_arbre >> type_int >> M >> graine;

    ParametresBS p = { S0, r_val, sigma_val, q_val };
    OptionType otype = (type_int == 0) ? PUT : CALL;

    // Bermudien
    std::vector<double> dates(N_berm);
    for (int i = 0; i < N_berm; ++i) dates[i] = T_max * (i + 1) / N_berm;
    Random rng(graine);
    BermudeNDatesLS opt_b(dates, K_val, otype);
    double pb, ic_inf, ic_sup;
    opt_b.monte_carlo_longstaff_schwartz(p, M, rng, pb, ic_inf, ic_sup);

    // Americain
    PricerArbreBinomial opt_a(N_arbre, K_val, T_max);
    double pa = (type_int == 0) ? opt_a.prix_put_americain(p) : opt_a.prix_call_americain(p);

    std::cout << pb << " " << pa << std::endl;
    return 0;
}
"""

with open("/tmp/tmp_main_div.cpp", "w") as f:
    f.write(tmp_main)

result = subprocess.run(
    "c++ -std=c++11 -O2 -o /tmp/pricer_div /tmp/tmp_main_div.cpp "
    "BlackScholes.cpp Random.cpp Bermudes3Dates.cpp PricerArbreBinomial.cpp -I.",
    shell=True, capture_output=True, text=True,
    cwd="/Users/louenmarx/Documents/test_c/CONTROLE_C/bermudean_option_project")
if result.returncode != 0:
    print("Erreur:", result.stderr); exit(1)

print("Compilation OK.")

# ========================================================================
# Simulations PUT + CALL
# ========================================================================
prix_berm_put = []; prix_amer_put = []
prix_berm_call = []; prix_amer_call = []

total = len(S0_range) * 2
done = 0

for S0 in S0_range:
    for type_int, berm_list, amer_list in [(0, prix_berm_put, prix_amer_put),
                                            (1, prix_berm_call, prix_amer_call)]:
        inp = f"{S0} {r} {sigma} {K} {T} {q} {N_berm} {N_arbre} {type_int} {M} {graine}\n"
        res = subprocess.run(["/tmp/pricer_div"], input=inp, capture_output=True, text=True)
        parts = res.stdout.strip().split()
        berm_list.append(float(parts[0]))
        amer_list.append(float(parts[1]))
        done += 1

    bar = "#" * done; espace = " " * (total - done)
    print(f"\r[{bar}{espace}] S0={S0:.0f}", end="", flush=True)

print("\nSimulations terminees.")

prix_berm_put = np.array(prix_berm_put)
prix_amer_put = np.array(prix_amer_put)
prix_berm_call = np.array(prix_berm_call)
prix_amer_call = np.array(prix_amer_call)

# ========================================================================
# Graphique 2x2 : PUT (gauche) + CALL (droite)
# ========================================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

for col, (ax_prix, ax_prime, prix_eur, prix_berm, prix_amer, intrinseque, label, color) in enumerate([
    (axes[0][0], axes[1][0], prix_eur_put, prix_berm_put, prix_amer_put, intrinseque_put, 'PUT', 'blue'),
    (axes[0][1], axes[1][1], prix_eur_call, prix_berm_call, prix_amer_call, intrinseque_call, 'CALL', 'green'),
]):

    prime_berm = prix_berm - prix_eur
    prime_amer = prix_amer - prix_eur

    # Prix
    ax_prix.fill_between(S0_range, prix_eur, prix_amer, alpha=0.08, color='green')
    ax_prix.fill_between(S0_range, prix_eur, prix_berm, alpha=0.12, color=color)
    ax_prix.plot(S0_range, prix_amer, 'g-s', markersize=4, linewidth=2, label='Americain (CRR)')
    ax_prix.plot(S0_range, prix_berm, f'{color[0]}-o', markersize=3, linewidth=2, label='Bermudien (LS)')
    ax_prix.plot(S0_range, prix_eur, 'r--', linewidth=2, label='Europeen (BS)')
    ax_prix.plot(S0_range, intrinseque, 'k:', linewidth=1.2, alpha=0.4, label='Intrinseque')
    ax_prix.axvline(x=K, color='gray', linestyle=':', alpha=0.5)
    ax_prix.set_ylabel('Prix (EUR)', fontsize=12)
    ax_prix.set_title(f'{label}  -  K={K:.0f}  r={r*100:.0f}%  sigma={sigma*100:.0f}%  q={q*100:.0f}%  T={T:.0f}ans',
                      fontsize=12, fontweight='bold')
    ax_prix.legend(fontsize=9, loc='upper right' if col == 0 else 'upper left')
    ax_prix.grid(True, alpha=0.3)
    ax_prix.set_xlim(S0_range[0], S0_range[-1])

    # Primes barres cote a cote
    w = (S0_range[1] - S0_range[0]) * 0.35
    ax_prime.bar(S0_range - w, prime_berm, width=w, color=color, alpha=0.7,
                 edgecolor='navy' if col == 0 else 'darkgreen', linewidth=0.5, label='Prime Bermudien')
    ax_prime.bar(S0_range + w, prime_amer, width=w, color='limegreen', alpha=0.7,
                 edgecolor='darkgreen', linewidth=0.5, label='Prime Americain')
    ax_prime.axhline(y=0, color='black', linewidth=1)
    ax_prime.axvline(x=K, color='gray', linestyle=':', alpha=0.5)
    ax_prime.set_xlabel('S0 (EUR)', fontsize=12)
    ax_prime.set_ylabel('Prime (EUR)', fontsize=12)
    ax_prime.set_title(f"Prime exercice anticipe - {label}", fontsize=11)
    ax_prime.legend(fontsize=9)
    ax_prime.grid(True, alpha=0.3, axis='y')
    ax_prime.set_xlim(S0_range[0], S0_range[-1])

plt.tight_layout()
plt.savefig('comparaison_europeen_vs_bermudien.png', dpi=150, bbox_inches='tight')
print("Graphique sauvegarde dans comparaison_europeen_vs_bermudien.png")
