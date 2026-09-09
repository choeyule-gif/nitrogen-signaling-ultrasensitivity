# Two separable sources of ultrasensitivity in bifunctional modification cycles

Code and data for

> Choi Y. *Two separable sources of ultrasensitivity in bifunctional modification cycles, and how a bacterial nitrogen cascade uses them.* PLOS Computational Biology (submitted, 2026).

The manuscript studies the effector-switched bifunctional network (one enzyme carries both the
modifying and the demodifying activity; a small effector bound to the enzyme selects between
them), proves that the enzyme concentration cancels from its steady state, that the sharpness of
the response factorises into a ladder factor and a switch factor, and applies the result to the
GlnD–PII and ATase–GS cycles of *E. coli* nitrogen assimilation.  Everything numerical in the
paper — the verifications of the theorems, the reanalysis of the published titration, the cascade
accounting and the eight figures — is produced by the scripts in this repository, and
`tests/test_reproduce.py` checks the output against the values quoted in the text.

**Archived release:** Zenodo DOI `10.5281/zenodo.<to be assigned>` (see `CITATION.cff`).

## Contents

| path | what it holds |
|---|---|
| `esbm/` | the model as a small Python package: compartmental matrix and kernel steady state of the raw network (`model.py`), closed forms proved in the paper, the figure style, and the loaders for `data/` |
| `data/` | inputs only: the twelve digitised points of the published titration and the published constants (`data/README.md` gives sources and the digitisation uncertainty model) |
| `scripts/` | one script per result; each prints what it reproduces and writes `results/<name>.json` |
| `results/` | generated numbers (`*.json`), logs, and the bootstrap / profile / MCMC arrays; `results/expected/` is the set the manuscript was written from |
| `figures/` | generated `Fig1.pdf` … `Fig8.pdf` |
| `S1_Data/` | generated Supporting Information: the values behind every figure panel (`Fig*.csv`) and the rate constants, totals and roots of every multistationary example (`S1_counterexamples.json`) |
| `tests/` | comparison of `results/*.json` with the manuscript |
| `extras/` | the Re *F*(iω) figure of the companion stability manuscript, which the Results cite |

## Requirements

Python ≥ 3.10 with `numpy`, `scipy`, `sympy`, `mpmath` and `matplotlib` (`requirements.txt` pins the
versions the manuscript was produced with: CPython 3.14.6 on macOS 26 / Apple silicon, numpy 2.5.2,
scipy 1.18.1, sympy 1.14.0, mpmath 1.3.0, matplotlib 3.11.1).  No compiled code, no GPU, no
network access.

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Reproducing the paper

```bash
python run_all.py          # everything, in dependency order (~40 min on a laptop)
python run_all.py --quick  # skips the four long searches (~8 min)
python tests/test_reproduce.py      # or: python -m pytest -q tests
```

`run_all.py` runs the scripts below in order and writes a log per script to `results/`.  Every
script can also be run on its own from any directory (`python scripts/fig2_sites.py`); the few
that depend on another's output say so in their docstring (`fig8_budget.py` and
`sensitivity_major1.py` need `cascade_budget.py`; `fig5_integer.py` needs `fit_titration.py`
and `mcmc_h.py`; `fig6_uniqueness.py` needs `counterexamples.py`).  All random draws are seeded;
the seeds are in the scripts.

## What each script reproduces

Equation and figure numbers refer to the manuscript.  The right-hand column is the statement in
the text that the script's printed output and `results/<name>.json` are checked against.

### Theory: verifications of the proved statements

| script | manuscript statement |
|---|---|
| `verify_symbolic.py` | computer-algebra checks with `n` symbolic: the chain survives three simultaneous dead ends (eq. nogo); the Möbius relation and the quartic with its extreme coefficients (eqs quartic, quarticco); the floored exponent-2 Hill identity (eq. hillid); the two-form chain with variable *l*₁/*l*₂ (eq. ratio) and the linear/quadratic form of *l*₂ under condition (i); eq. omcov on 4,000 random ladders; the sign counts of det *J* (18/0 at *n* = 1, 44/7 at *n* = 2); the active bridging complex at *n* = 1 |
| `verify_adair.py` | eq. adair against the kernel steady state over 1,500 parameter sets; the two-form network over a thousandfold range of enzyme; *n*_H = *h* exactly for ε = *A e l^h*; ν → 1.999 at *n* = 2 with *K*₁ = 10⁻³, *K*₂ = 10³ |
| `verify_nu_bound.py` | eq. nubound in 50–80-digit arithmetic: no exceedance of *n* over 28 decades of *c_i*, the bound approached to 10⁻²⁰, deviation of ν from 1 below 6 × 10⁻⁴³ for a binomial ladder and non-constant under any perturbation |
| `verify_deadend_60digit.py` | Fig 1b in 60-digit arithmetic: the chain holds to 1.3 × 10⁻²⁹ with a dead end swept over eleven decades; the ratios drift by 5.4 × 10⁻⁵ at *T_E*/*T_L* = 1.1 × 10⁻⁴ and 5.4 × 10⁻¹⁰ at 1.1 × 10⁻⁹ (eq. sandwich) |
| `sequestration_sweeps.py` | eq. zeroorder: 0.988 → 0.976 on doubling the enzyme; the Bhat/Ghat/*A* sweep; the target sweep *K_m*/60 … 100 *K_m* giving [0.981, 0.9996] |
| `deficiency_acr.py` | δ = 8 − 3 − 4 = 1 at *n* = 1, class deficiencies 0, no Shinar–Feinberg pair, δ = *n* in general; ACR of *s*₁ in the complex-mediated design; the buffering example (*e* 0.1128 → 0.0326, *s*₀ 0.2669 → 4.7647) |
| `counterexamples.py` | Fig 6b and S1 Appendix: three positive steady states with condition (i) broken by a spread of 2.89 in ρ_i at *n* = 4, and with condition (ii) broken at *n* = 4, 5, 6 (elasticity of ω up to 4.98 at *n* = 5); all roots refined in 60-digit arithmetic; writes `S1_Data/S1_counterexamples.json` |
| `mono_adversarial.py` | the adversarial search for non-monotone *T_L*(*l*): 21 pairs (*n*, *h*), 60 Nelder–Mead starts each, every optimum positive in 50-digit arithmetic between 9 × 10⁻¹⁵ and 10⁻¹⁰ |

### Application: reanalysis of the published data

| script | manuscript statement |
|---|---|
| `fit_titration.py` | Table tab:ic (SSE, RMSE, AIC, BIC for *h* = 1…4 and free), Shapiro–Wilk *W* = 0.964, profile interval [1.86, 2.20], standardised residuals 6.9 / 1.9 / 4.7, the floorless model (2Δ𝓛 = 31.9, *h* = 1.42), the refractory-fraction model (*q* = 0.15, *h* = 1.96, 2Δ𝓛 = 1.30), end-point drops 2.06 / 2.03 / 2.04, *U*_min = 0.45 ± 0.04, η = 1.63 and 1.71 |
| `mcmc_h.py` | posterior of *h*: acceptance 0.245, R̂ = 1.0002, ESS 1.6 × 10⁴, median 2.025, 95 % interval [1.81, 2.28], *P*(*h* < 1.5) = 3 × 10⁻⁵ |
| `cascade_budget.py` | R4/R5: readout factors 2.44 / 1.87 / 1.76 and composites 4.92 / 3.78 / 3.55 (linear reading), 2.32 / 1.67 / 1.57 and 4.67 / 3.37 / 3.17 (complementary); supremum 2.4399 and matched factor 1.25935, both independent of *h*; predicted half-point 1.08 mM; the digitisation bootstrap (536 of 600 usable) giving residuals 0.92 [0.83, 1.03], 1.70 [1.57, 1.84], 1.47 [1.35, 1.58] |
| `sensitivity_major1.py` | error on the three published cascade coefficients: the peak survives with probability 0.98 / 0.85 / 0.76 at 5 / 10 / 15 % relative error; 6.46 would have to fall 13.9 % to 5.57 |
| `nu2_major2.py` | letting ν₂ float: residuals 0.74 / 1.38 / 1.20 at ν₂ = 1.34 and 0.51 / 0.96 / 0.83 at ν₂ = 2.2; the peak survives at every ν₂ |
| `mixture_bias.py` | R1: 15 % of competent trimers (27 % of modifiable subunits) di- or trivalent; at ν(*n* = 3) = 1.40 the mixture raises the monovalent coefficient by 2.8 % |
| `conversion_factor.py` | Methods: the conversion from a reported coefficient to the local ν on the three ladder families (1.02–1.04), on the full two-coordinate locus, and its independence of *h* |
| `swing_profile.py` | R3: θ₀ = 0.989, θ_∞ = 0.151, *w*₀/K̂ = 0.011, swing ≈ 5 × 10², 95 % interval [175, ∞), 2Δ𝓛 = 54 at swing = *R*_rev, inhibition ≥ 26-fold |
| `direct_input.py` | the first candidate of 'What is left': additivity of the two log-slopes and the ten-to-ninety gain of a direct glutamine input (a re-implementation; see the note below) |
| `ratios_table.py` | Table tab:ratios: *T_E*/*T_L* = 9 × 10⁻⁵ (first cycle), 2 × 10⁻⁴ at *h* = 2 for the refit, 0.003 / 0.02 / 0.2 (second cycle) |

### Figures

| script | figure |
|---|---|
| `fig1_invariance.py` | Fig 1: enzyme swept 5,000-fold (agreement with *K_i*/*l* to 5 × 10⁻¹⁴); dead end over eleven decades (drift 5 × 10⁻¹⁰) |
| `fig2_sites.py` | Fig 2: θ(*l*) and the local coefficient (= 1 to 4 × 10⁻⁸) for *n* = 1, 2, 3, 12; the total-fraction coefficient against *T_E*/*T_S* |
| `fig3_plateau.py` | Fig 3: the surface of the local coefficient over both totals, flat to 2.5 × 10⁻³ (0.99745–1) on the interior |
| `fig4_plane.py` | Fig 4: the (ν, η) plane with the two ceilings and the measured point |
| `fig5_integer.py` | Fig 5: the exact Hill functions, the fits at fixed *h*, the standardised residuals, the profile and the posterior |
| `fig6_uniqueness.py` | Fig 6: the residual 𝓡(*l*) with one zero (conditions satisfied) and three zeros (condition (i) broken) |
| `fig7_tests.py` | Fig 7: the titrations predicted at two enzyme and two target levels under the three accounts (1.000 in every condition; 1.08 against 1.26 for the bridging complex) |
| `fig8_budget.py` | Fig 8: the budget, the residual with its bootstrap intervals, and the readout factor over 24 decades of K̂₂ |

## Notes on reproducibility

* Every number in the Results, Methods and figure captions is produced by a script here and
  compared in `tests/test_reproduce.py` at the rounding used in the text.  The long searches
  (`mono_adversarial.py`, `counterexamples.py`, `mcmc_h.py`, the profile in `fit_titration.py`)
  are seeded and deterministic on one machine; across numpy/scipy versions their qualitative
  conclusions are what the tests check.
* `direct_input.py` is a re-implementation: the original interactive calculation behind the three
  additivity numbers quoted in the text (3.0009 + 0.4913 = 3.4922) was not preserved.  The script
  states its construction, reproduces the exactness of the additivity and the gain range
  0.95–1.2 at the measured constants, and reports what its construction gives elsewhere.
* Steady states are always the kernel of the compartmental matrix (an exact linear solve), never
  a multistart Newton iteration; in this network Newton iteration produces spurious roots from
  cancellations of order *u^n*.  In double precision the kernel solve is limited by the condition
  number of the matrix (of order 10⁶ for rate constants spanning three decades), so
  `verify_adair.py` reports three numbers: the wide random family in double precision, a
  well-conditioned family in double precision, and the wide family with the kernel solved in
  40-digit arithmetic, which is the figure that tests the algebra rather than the arithmetic.
* `conversion_factor.py` evaluates the conversion on the three ladder families named in the
  Methods and, separately, along the full two-coordinate locus of three-site ladders that
  reproduce the observed ratio; the locus range is wider than the three-family envelope and is
  reported as such in `results/conversion_factor.json`.
* The twelve titration points were read from the printed figure; the uncertainty model
  (±0.06 in ordinate, ±6 % in abscissa) is stated in `data/README.md` and used by every bootstrap.

## License and citation

Code: MIT (`LICENSE`).  Please cite the manuscript and the Zenodo record (`CITATION.cff`).
