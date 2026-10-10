# Kinetic Parameter Estimation and Uncertainty Quantification for a Batch Reaction A → B

Estimation of kinetic parameters (rate constant, activation energy, reaction order) from noisy concentration–time data using nonlinear least squares in MATLAB, with uncertainty quantification and a Monte Carlo noise-sensitivity study.

**Author:** Lyna LAZZOUZI, École Nationale Polytechnique d'Alger (Chemical Engineering)
**Context:** Project 01 of a personal Chemical Engineering Research & Simulation portfolio (MATLAB → Python).

---

## 1. Problem Statement

Sizing or operating a chemical reactor requires reliable kinetic parameters, but they cannot be measured directly: they must be inferred by fitting a model to concentration measurements, which are always noisy.

**Research question:** How accurately can the rate constant *k* and the Arrhenius parameters (*k₀*, *Eₐ*) be recovered from noisy concentration data, and how does the measurement noise level affect the quality of the estimates?

Synthetic data with known "true" parameters are used on purpose: this makes it possible to check that the method recovers the right values and to quantify the estimation error.

## 2. Objectives

- Simulate a first-order batch reaction A → B and generate noisy synthetic data
- Estimate *k* by nonlinear least squares (`ode45` + `lsqnonlin`)
- Quantify the uncertainty (95 % confidence intervals) and validate the fit (residuals, RMSE, R²)
- Study the influence of the noise level (Monte Carlo, 200 runs per level)
- Estimate *k₀* and *Eₐ* from experiments at several temperatures (Arrhenius law)
- Extend to a free reaction order *n* and discuss parameter identifiability

## 3. Mathematical Model

Closed (batch) reactor, isothermal, perfectly mixed, constant volume, single reaction A → B.

Mass balance:

$$\frac{dC_A}{dt} = -k\,C_A^{\,n}, \qquad C_A(0) = C_{A0}$$

For first order ($n = 1$) the analytical solution is $C_A(t) = C_{A0}\,e^{-kt}$, i.e. $\ln C_A = \ln C_{A0} - kt$.

Arrhenius law:

$$k(T) = k_0\,e^{-E_a/(RT)}$$

For the global fit, the law is reparametrized around a reference temperature $T_{ref}$ to reduce the strong correlation between $k_0$ and $E_a$:

$$k(T) = k_{ref}\exp\left[-\frac{E_a}{R}\left(\frac{1}{T}-\frac{1}{T_{ref}}\right)\right], \qquad k_0 = k_{ref}\,e^{E_a/(RT_{ref})}$$

## 4. Methodology

| Element | Choice |
|---|---|
| ODE solver | `ode45`, RelTol = 1e-8, AbsTol = 1e-10 |
| Estimation | `lsqnonlin` (Optimization Toolbox), lower bound k ≥ 0 |
| Confidence intervals | `nlparci` (95 %, from the Jacobian) |
| Noise model | Multiplicative Gaussian noise (relative standard deviation 1 % to 20 %) |
| Monte Carlo | 200 repetitions for each of 5 noise levels (1000 fits) |
| Arrhenius, method A | Fit *k* at each temperature, then linear regression of ln *k* vs 1/*T* |
| Arrhenius, method B | Global fit of $k_{ref}$ and $E_a$ on all data (residuals in log form) |
| Free order | Fit of [*k*, *n*] on three initial concentrations (1, 2, 4 mol/L) |
| Reproducibility | Fixed random seeds (`rng(1)` to `rng(4)`) |

Base case: $k = 0.35\ \mathrm{min^{-1}}$, $C_{A0} = 2\ \mathrm{mol/L}$, 16 time points from 0 to 15 min, 3 % relative noise.
Arrhenius case: $k_0 = 10^7\ \mathrm{min^{-1}}$, $E_a = 50\ \mathrm{kJ/mol}$, $T$ = 300, 310, 320, 330, 340 K, $T_{ref}$ = 320 K, 15 points per temperature.

## 5. Results

### 5.1 Estimation of k (3 % noise, first order)

| Quantity | Value |
|---|---|
| True *k* | 0.3500 min⁻¹ |
| Estimated *k* (nonlinear fit) | 0.3527 min⁻¹ |
| 95 % confidence interval | [0.3445 ; 0.3610] |
| Relative error | 0.78 % |
| RMSE | 0.0179 mol/L |
| R² | 0.9990 |
| *k* from linearized regression (ln *C*ₐ vs *t*) | 0.3491 min⁻¹ |

![Kinetic fit and residuals](figures/fig01_fit_k.png)

### 5.2 Influence of the noise level (Monte Carlo, 200 runs per level)

| Noise σ | Bias on *k* | Std. deviation | CV (%) |
|---|---|---|---|
| 1 % | 0.00025 | 0.0016 | 0.45 |
| 3 % | 0.00061 | 0.0045 | 1.27 |
| 5 % | 0.00165 | 0.0078 | 2.24 |
| 10 % | 0.00190 | 0.0136 | 3.88 |
| 20 % | 0.00538 | 0.0334 | 9.54 |

![Dispersion of k versus noise level](figures/fig02_noise_boxplot.png)

### 5.3 Arrhenius parameters

| Method | *Eₐ* (kJ/mol) | *k₀* (min⁻¹) |
|---|---|---|
| True values | 50.00 | 1.000 × 10⁷ |
| A: two-step regression | 50.15 | 1.058 × 10⁷ |
| B: global fit | 49.94 | 9.797 × 10⁶ |

The 95 % confidence interval of *Eₐ* (method B) is [49.8 ; 50.1] kJ/mol and contains the true value.

![Arrhenius plot](figures/fig03_arrhenius.png)

### 5.4 Free reaction order (true k = 0.2, n = 2)

| Data used | *k* (95 % CI) | *n* (95 % CI) | CI width on *n* |
|---|---|---|---|
| 3 initial concentrations | 0.199 [0.197 ; 0.202] | 1.991 [1.959 ; 2.023] | 0.064 |
| 1 initial concentration (2 mol/L) | 0.200 [0.196 ; 0.204] | 1.999 [1.910 ; 2.089] | 0.179 |

![Free reaction order fit](figures/fig04_order_n.png)

## 6. Validation

- **Solver vs analytical solution:** maximum difference between `ode45` and $C_{A0}e^{-kt}$ is 7.5 × 10⁻¹⁰ mol/L.
- **Estimates vs known truth:** *k*, *Eₐ*, *k₀*, and *n* are recovered within a few percent, and the true values lie inside the 95 % confidence intervals in this run.
- **Convergence:** `exitflag = 1` for the base fit.
- **Residuals:** no strong systematic trend. Their absolute amplitude is larger at early times, as expected with relative noise (the error is proportional to the concentration).

## 7. Discussion

- **Noise and uncertainty.** The standard deviation of *k* grows with the noise level, and the coefficient of variation stays at roughly 0.4 to 0.5 times the relative noise (for example, 20 % noise gives a CV of about 9.5 %). The bias remains small: it is between about 14 % and 21 % of the standard deviation at every noise level, so the estimation error is dominated by variance rather than by bias.
- **Linearization vs nonlinear fit.** Both give values close to the truth (0.3491 and 0.3527 versus 0.3500). They are not identical because they minimize different criteria (errors on *C*ₐ versus on ln *C*ₐ). With relative noise, the log regression is reasonable here; the nonlinear fit remains more general because it applies to any simulable model.
- **Arrhenius.** The global fit (method B) gives *Eₐ* and *k₀* closer to the truth in this run and directly provides confidence intervals. The reparametrization around *T*ᵣₑf reduces the correlation between the two parameters. This comparison comes from a single noise realization.
- **Reaction order.** *k* and *n* are recovered accurately even with one initial concentration, but the confidence interval on *n* is about 2.8 times wider (0.179 versus 0.064). Using several initial concentrations therefore improves identifiability substantially.

## 8. Limitations

- Synthetic data only: real measurements include systematic errors (calibration, sampling) that are not represented.
- Ideal closed, isothermal, perfectly mixed reactor with a single reaction; no transport limitations.
- Gaussian, independent noise with a known relative standard deviation.
- Confidence intervals rely on a local linearization of the model (Jacobian approximation).
- The Arrhenius comparison uses a single noise realization, and the effect of the temperature range on *Eₐ* was not quantified.

## 9. Future Work

- Apply the method to real experimental data
- Study the effect of the temperature range and number of points on *Eₐ* (Monte Carlo)
- Handle series and parallel reactions (several species)
- Use a Bayesian approach for more reliable intervals
- Reuse these building blocks in the next projects: CSTR/PFR comparison, then the exothermic CSTR

## How to Run

**Requirements:** MATLAB R2020a or later (for `boxchart`), Optimization Toolbox, Statistics and Machine Learning Toolbox.

1. Clone or download this repository and set the project folder as the MATLAB current folder.
2. Open `Projet01_Report_EN.mlx` and click **Run All** (the Monte Carlo section performs 1000 fits).
3. Figures are saved in `figures/`, tables in `results/`, generated data in `data/`.

## Repository Structure

```
Projet01_Cinetique/
├── README.md
├── Projet01_Report_EN.mlx      Live Script: full study (text, code, results)
├── report/                     PDF export of the Live Script
├── src/
│   ├── model_first_order.m     First-order ODE model
│   ├── model_order_n.m         Order-n ODE model
│   ├── global_residual.m       Residuals of the global Arrhenius fit
│   └── resid_n.m               Residuals of the free-order fit
├── data/                       Generated synthetic data
├── figures/                    Result figures
└── results/                    Monte Carlo results (CSV)
```

## References

- H. S. Fogler, *Elements of Chemical Reaction Engineering*, Pearson.
- O. Levenspiel, *Chemical Reaction Engineering*, Wiley.
- MATLAB documentation: `ode45`, `lsqnonlin`, `nlparci`, `boxchart`.
