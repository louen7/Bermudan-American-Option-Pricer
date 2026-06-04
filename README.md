# Bermudean Option Pricing Project

Comparison of **European**, **Bermudean**, and **American** put option pricing using three numerical methods. The project validates the theoretical price ordering:

```
European <= Bermudean <= American
```

and demonstrates the convergence of the Bermudean price toward the American price as the number of exercise dates increases.

---

## Pricing Methods

| Option Type | Method | Description |
|---|---|---|
| European | Black-Scholes | Closed-form analytical solution |
| Bermudean | Longstaff-Schwartz | Monte Carlo with least-squares regression for continuation value estimation |
| American | Cox-Ross-Rubinstein | Binomial tree with backward induction and early exercise |

---

## Project Structure

```
.
├── Makefile                    # Build configuration
├── main.cpp                    # Main program: comparison tables and convergence study
├── BlackScholes.hpp / .cpp     # Black-Scholes model (analytical pricing, MC simulation)
├── Random.hpp / .cpp           # Random number generation (Mersenne Twister, Box-Muller)
├── BermudeNDatesLS.hpp / .cpp  # Longstaff-Schwartz method for N-dates Bermudean options
├── Bermudes2Dates.hpp / .cpp   # 2-dates Bermudean option (analytical boundary + antithetic MC)
├── Bermudes3Dates.hpp          # 3-dates Bermudean option (Longstaff-Schwartz)
└── PricerArbreBinomial.hpp/.cpp # CRR binomial tree for American options
```

---

## Build & Run

```bash
make            # Compile the project
./projet_bermude   # Run the program
make clean      # Remove object files and executable
```

**Requirements:** C++11 compatible compiler (`c++`, `g++`, or `clang++`).

---

## Default Parameters

| Parameter | Value | Description |
|---|---|---|
| `K` | 100.0 | Strike price |
| `r` | 15% | Risk-free interest rate |
| `σ` | 15% | Volatility |
| `T` | 10 years | Maturity |
| `M` | 300 000 | Number of Monte Carlo paths |
| `N_berm` | 20 | Exercise dates (Bermudean) |
| `N_tree` | 300 | Time steps (binomial tree) |
| `seed` | 42 | Random seed (reproducibility) |

---

## Program Output

### Part 1 — Price Comparison Table

For 15 spot prices (S0 = 30 to 150), the program displays:

- **European** price (Black-Scholes closed-form)
- **Bermudean** price (Longstaff-Schwartz MC) with 90% confidence interval
- **American** price (CRR binomial tree)
- **Bermudean premium** = Bermudean − European
- **American premium** = American − European
- **Gap A − B** = American − Bermudean

### Part 2 — Convergence Study

For S0 = 60, the Bermudean price is computed with increasing numbers of exercise dates (N = 1, 2, 3, 5, 10, 15, 20, 30, 50, 75, 100). The gap with the American price shrinks as N grows, confirming convergence.

---

## Key Algorithms

### Longstaff-Schwartz Method

1. **Forward phase** — Simulate M stock price paths across all exercise dates
2. **Backward phase** — At each exercise date:
   - Filter in-the-money (ITM) paths
   - Estimate continuation values via quadratic least-squares regression
   - Decide: exercise immediately or continue
3. **Pricing** — Discount optimal cash-flows to t = 0, compute 90% CI

Data centering is used in the regression to ensure numerical stability.

### CRR Binomial Tree

1. Build forward stock price tree with parameters u = exp(σ√Δt), d = 1/u
2. Initialize option values at maturity
3. Backward induction: max(exercise, continuation) at each node

### 2-Dates Analytical Boundary

The critical exercise boundary K_bar at T1 is found by bisection on:

```
f(S) = Payoff(S, T1) − ContinuationValue(S, T1)
```

Antithetic variance reduction is available for the Monte Carlo simulation.

---
