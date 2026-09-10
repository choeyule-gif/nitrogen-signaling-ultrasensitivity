# Reproducibility record

Every quantitative statement of the manuscript, the script that produces it, the value the
script gives, and whether the text agrees. The table below was written against the longer first
version of the manuscript, whose numbers the revised version carries over unchanged; where the
revised version states a quantity differently, `results/revision_numbers.json` holds the value it
quotes and `tests/test_reproduce.py` checks both.  Values were regenerated from a clean checkout on
2026-09-09 (CPython 3.14.6, numpy 2.5.2, scipy 1.18.1, sympy 1.14.0, mpmath 1.3.0).
"agrees" means agreement at the rounding used in the text; the last section lists every
statement that does **not** agree and what the code supports instead.

## Theory

| statement in the text | script | reproduced value | |
|---|---|---|---|
| Fig 1a: enzyme swept over a 5,000-fold range, ratios agree with *K_i*/*l* to 5 × 10⁻¹⁴ | `fig1_invariance.py` | 5,278-fold; 5.0 × 10⁻¹⁴; drift 1.5 × 10⁻¹⁴ | agrees |
| Fig 1b: dead end over eleven decades at *T_E*/*T_L* = 1.1 × 10⁻⁹, ratios drift by 5 × 10⁻¹⁰ | `fig1_invariance.py` | 5.4 × 10⁻¹⁰ | agrees |
| eq. nogo in 60-digit arithmetic: chain to 1.3 × 10⁻²⁹; ratios drift 5.4 × 10⁻⁵ at 1.1 × 10⁻⁴ and 5.4 × 10⁻¹⁰ at 1.1 × 10⁻⁹ | `verify_deadend_60digit.py` | 1.25 × 10⁻²⁹; 5.37 × 10⁻⁵; 5.37 × 10⁻¹⁰ | agrees |
| eq. nogo symbolically with three simultaneous dead ends | `verify_symbolic.py` A | no dead-end constant in *s*₁*l*/*s*₀ | agrees |
| Fig 2b: local coefficient equal to 1 to 4 × 10⁻⁸ | `fig2_sites.py` | 4.4 × 10⁻⁸ | agrees |
| eq. adair verified to 10⁻¹⁵ over 1,500 parameter sets | `verify_adair.py` | see discrepancy 1 | **differs** |
| eq. nubound: no ν above *n* over 28 decades; bound approached to 10⁻²⁰; binomial deviation below 6 × 10⁻⁴³; perturbation makes ν non-constant | `verify_nu_bound.py` | 0 exceedances; gap 1.0 × 10⁻²⁰ at ε = 10⁻²⁰; 5.3 × 10⁻⁴³; non-constant | agrees |
| *n* = 2, *K*₁ = 10⁻³, *K*₂ = 10³: ν rises to 1.999 | `verify_adair.py` | 1.999000 | agrees |
| eq. integer: *n*_H = *h* exactly for *h* = 1, 2, 3 | `verify_adair.py` | *h* to 10 decimals | agrees |
| two-form network: chain variable *l*₁/*l*₂, enzyme cancels (*n* = 4, symbolic) | `verify_symbolic.py` D | identity holds | agrees |
| two-form network numerically: ratios move by 8 × 10⁻¹⁶ over a thousandfold range of enzyme | `verify_adair.py` | see discrepancy 1 | **differs** |
| eq. zeroorder: 0.988 at *T_E*/*T_S* = 0.1, 0.976 on doubling | `sequestration_sweeps.py` | 0.98779; 0.97615 | agrees |
| sweep of Bhat, Ghat over two decades and *A* over three: 0.92–0.999, change on doubling 0.001–0.08 | `sequestration_sweeps.py` | see discrepancy 5 | **partly** |
| target sweep *K_m*/60 … 100 *K_m*: 0.981–0.9996 | `sequestration_sweeps.py` | 0.9813–0.9996 | agrees |
| Fig 3: interior flat between 0.99745 and 1, width 2.5 × 10⁻³ | `fig3_plateau.py` | 0.997451–1.000000; 2.5 × 10⁻³ | agrees |
| eq. defic: δ = 8 − 3 − 4 = 1; class deficiencies 0; no Shinar–Feinberg pair; |𝒞| = 4*n*+4, ℓ = 3, rank 3*n*+1, δ = *n* | `deficiency_acr.py` | all as stated, *n* = 1…6 | agrees |
| ACR in the complex-mediated design; buffering *e* 0.1128 → 0.0326, *s*₀ 0.2669 → 4.7647 in the effector-switched one | `deficiency_acr.py` | *s*₁ constant to 10⁻¹⁰ over six pairs of totals; 0.112825 → 0.0326392; 0.266936 → 4.76472 | agrees |
| eq. mobius, eq. quartic, eq. quarticco (degree 4 for symbolic *n*; extreme coefficients) | `verify_symbolic.py` B | degree 4; both coefficients; factorised form equals the numerator up to (1+ρ)⁻² | agrees |
| eq. hillid: exact floored Hill function of exponent 2 | `verify_symbolic.py` C | residual 0 | agrees |
| eq. omcov | `verify_symbolic.py` F | 0 violations in 4,000 ladders | agrees |
| det *J*: 18 positive, 0 negative at *n* = 1; 44 positive, 7 negative at *n* = 2 | `verify_symbolic.py` G | 18/0; 44/7 | agrees |
| bridging complex at *n* = 1: ratio free of *e*, carries *s*₀, reduces to the chain | `verify_symbolic.py` H | all three | agrees |
| two-form network under (i): *l*₂ = linear/quadratic in *e* | `verify_symbolic.py` E | degrees (1, 2) | agrees |
| Fig 6b: condition (i) broken by a spread of 2.89 at *n* = 4, three positive steady states, residuals below 10⁻⁶⁰ | `counterexamples.py` | *n* = 4, spread 2.890, three roots, residuals ≤ 10⁻⁶⁰ | agrees |
| S1 Appendix: condition (ii) broken at *n* = 4, 5, 6; elasticity of ω reaches 4.98 | `counterexamples.py` | *n* = 6, 5, 4; max elasticity 5.95, 4.98, 2.98 | agrees |
| adversarial search: 21 pairs, 60 starts, cap 2,500; double precision reaches order −10⁻⁹; all optima positive in 50 digits, between 9 × 10⁻¹⁵ and 10⁻¹⁰ | `mono_adversarial.py` | min −1.1 × 10⁻⁸ in double precision; 8.5 × 10⁻¹⁵ to 1.1 × 10⁻¹⁰ verified | agrees (the "order −10⁻⁹" is −1.1 × 10⁻⁸ at worst) |

## Application

| statement in the text | script | reproduced value | |
|---|---|---|---|
| Table tab:ic (SSE 0.3658/0.0149/0.1193/0.2919/0.0149; ΔAIC 38.38/0/24.93/35.67/1.95) | `fit_titration.py` | identical to four decimals | agrees |
| free-*h* optimum 2.02 (*S*₀.₅ 0.56 mM); Shapiro–Wilk *W* = 0.964, *p* = 0.83 | `fit_titration.py` | 2.0179; 0.5595; 0.9635; 0.832 | agrees |
| profile: minimum 2.02, interval [1.86, 2.20], 38.4 at *h* = 1, 25.0 at *h* = 3 | `fit_titration.py` | 2.02; [1.86, 2.20]; 38.4; 25.0 | agrees |
| standardised residuals 6.9 (h=1), 4.7 (h=3), inside 1.9 (h=2), σ̂ = 0.0407 | `fit_titration.py` | 6.92; 4.70; 1.83; 0.0407 | agrees |
| floorless model: 2Δ𝓛 = 31.9, SSE fourteenfold, *h* = 1.42 | `fit_titration.py` | 31.9; 14.2; 1.42 | agrees |
| refractory model: *q* = 0.15, *h* = 1.96, 2Δ𝓛 = 1.30, marginally better on AIC | `fit_titration.py` | 0.148; 1.96; 1.30; ΔAIC −0.70 | agrees |
| end-point drops 2.06 / 2.03 / 2.04; *U*_min = 0.45 ± 0.04 | `fit_titration.py` | 2.065 / 2.030 / 2.045; 0.454 ± 0.039 | agrees |
| η = 1.63 at θ = ½ and 1.71 at *l* = *S* for the *h* = 2 fit | `fit_titration.py` | 1.635; 1.715 | agrees |
| MCMC: acceptance 0.245, R̂ 1.0002, ESS 1.6 × 10⁴, median 2.025, [1.81, 2.28], *P*(*h* < 1.5) = 3 × 10⁻⁵, none below 1 | `mcmc_h.py` | 0.245; 1.00016; 15,768; 2.0250; [1.809, 2.282]; 3.4 × 10⁻⁵; 0 | agrees |
| bootstrap of the digitisation (4,000 replicates): 2.01 [1.71, 2.41] | `fit_titration.py` (Gaussian, `cascade_budget.py` uniform) | the 4,000-replicate Gaussian bootstrap is in `sim/` of the working folder, not re-run here; the uniform 600-replicate bootstrap gives 2.03 [1.87, 2.25] | see discrepancy 6 |
| R1: 15 % of competent trimers, 27 % of modifiable subunits; bias 2.8 % at ν(*n* = 3) = 1.40 | `mixture_bias.py` | 15.0 %; 26.5 %; 1.0276 | agrees |
| conversion 1.026–1.037 (*n*_H^fit) and 1.036–1.059 (*n*_H⁹⁰), quoted 1.04 ± 0.02; exactly independent of *h* | `conversion_factor.py` | see discrepancy 3 | **partly** |
| R3: θ₀ = 0.989, θ_∞ = 0.151, *w*₀/K̂ = 0.011, swing ≈ 5 × 10², [175, ∞), 2Δ𝓛 = 54 at *R*_rev, ≥ 26-fold | `swing_profile.py` | 0.989; 0.151; 0.0109; 516 (ML 501 on the grid); lower bound 174.5, open above; 54.3; 25.7 | agrees |
| R4: readout factors 2.44/1.87/1.76, composites 4.92/3.78/3.55 (linear); 2.32/1.67/1.57, 4.67/3.37/3.17 (complementary) | `cascade_budget.py` | 2.4355/1.8730/1.7602, **4.9145**/3.7794/3.5518; 2.3165/1.6678/1.5711, 4.6744/3.3654/3.1703 | agrees except 4.92 (discrepancy 2) |
| matched factor 1.26 and predicted half-point 1.08 mM; supremum 2.44, *h*-independent 2.43990, matched 1.25935 | `cascade_budget.py` | 1.25935; 1.082 mM; 2.43989 for *h* = 1…4; 1.25935 | agrees (last digit of 2.43990: computed 2.439885) |
| factor runs from 1.15 to 2.44 over 24 decades | `cascade_budget.py` | 1.146 (complementary) / 1.186 (linear) to 2.4399 | agrees |
| R5: residuals 0.91/1.71/1.47; bootstrap 536 usable: 0.92 [0.83, 1.03], 1.70 [1.57, 1.84], 1.47 [1.35, 1.58] | `cascade_budget.py` | 0.912/1.709/1.472; 536; 0.924 [0.829, 1.025], 1.702 [1.565, 1.835], 1.465 [1.353, 1.576] | agrees |
| Major 1: middle value would have to fall 13.9 % to 5.57; peak survives with 0.98/0.85/0.76 at 5/10/15 % | `sensitivity_major1.py` | 13.9 %; 5.57; 0.982/0.854/0.759 | agrees |
| Major 2: ν₂ = 1.34 → 0.74/1.38/1.20; ν₂ = 2.2 → 0.51/0.96/0.83; shape invariant | `nu2_major2.py` | 1.341 → 0.735/1.380/1.196; 2.234 → 0.509/0.955/0.830; peak at every ν₂ | agrees |
| direct glutamine input: additivity 3.0009 + 0.4913 = 3.4922; gain 0.95–1.2 at *R* = 3.7; close to one and usually below it over four decades | `direct_input.py` | see discrepancy 4 | **partly** |
| Table tab:ratios: 9 × 10⁻⁵; 3.6 × 10⁻³; 0.003/0.02/0.2; free effector to 2 × 10⁻⁴ at *h* = 2 | `ratios_table.py` | 9.1 × 10⁻⁵; 0.0028/0.020/0.20; 1.8 × 10⁻⁴ (= *h* × *T_E*/*T_L*) | agrees (3.6 × 10⁻³ is transcribed from the paper's concentrations, not computed) |
| Fig 7: 1.000 in every condition; bridging complex 1.08 against 1.26 | `fig7_tests.py` | 1.000 ×4; 1.080 / 1.259 | agrees |
| Fig 4: ν η = 1.11 × 1.82 = 2.02; η ceilings for the Mg²⁺ and Mn²⁺ affinity ratios | `fig4_plane.py` | 2.020; 1.033 (ratio 1.14), 1.367 (ratio 4.67) | agrees |

## Statements that the code does not reproduce as written

1. **Precision of the numerical checks of eq. adair and of the two-form network.**  The text says
   "to 10⁻¹⁵ over 1,500 parameter sets" and "the ratios moving by 8 × 10⁻¹⁶ over a thousandfold
   range of total enzyme".  With the kernel of the compartmental matrix solved in double precision
   over rate constants spanning three decades, the worst relative deviation is 6 × 10⁻⁹ (1,500
   sets, dynamic range < 10⁸) and 8 × 10⁻¹¹ on a well-conditioned family (rate constants over one
   decade); the two-form ratios agree with eq. ratio to 4 × 10⁻¹⁰ and drift by 5 × 10⁻¹⁰.  These are
   condition-number limits of the linear solve, not deviations of the identities: with the same
   1,500 sets (no conditioning filter) and the kernel solved in 40-digit arithmetic the worst
   deviation is 4 × 10⁻²⁴, and the two-form ratios agree with eq. ratio to 5 × 10⁻³⁵ and drift by
   8 × 10⁻³⁵ over the thousandfold range.  **Suggested wording:** "to 10⁻²⁴ over 1,500 parameter
   sets with the kernel solved in 40-digit arithmetic (10⁻⁹ in double precision, limited by the
   conditioning of the linear solve)" and "moving by 10⁻³⁴ over a thousandfold range of enzyme".
2. **Composite 4.92.**  The computed composite at 36 µM PII is 4.9145, which rounds to 4.91; the
   residual 0.91 quoted with it is unaffected.
3. **Conversion factor.**  On the three ladder families named in the Methods the two estimators
   give the ranges reported in `results/conversion_factor.json` (`fit_range`, `range_range`), inside
   1.04 ± 0.02, and the conversion is independent of *h* to 10⁻¹² once the fitting window is scaled
   with *h*.  Along the full two-coordinate locus of three-site ladders reproducing the observed
   ratio 1.1106, however, the conversion runs from 1.024 to 1.125 (*n*_H^fit) and from 1.036 to
   1.225 (*n*_H⁹⁰); the larger values belong to ladders whose doubly modified state is strongly
   depleted (*c*₂ ≪ 3), for which a Hill function fits the trimer curve poorly (RMSE up to 0.02 on
   the θ scale, against 0.0015 for the ladders at the low end).  The Methods
   already say the conversion "varies by order 10⁻²" over admissible ladders; the envelope
   1.04 ± 0.02 holds for the three named families and not for the whole locus.
4. **Direct glutamine input.**  The three additivity numbers (3.0009, 0.4913, 3.4922) came from an
   interactive calculation whose script was not kept; `direct_input.py` re-implements the check
   with a stated construction (single-site activation of the forward activity, φ = (1 + *R G*/*K_G*)/(1 + *G*/*K_G*),
   K̂₂ pinned by the PII path alone) and reproduces additivity exactly (4.2168 + 0.3039 = 4.5207 at
   the composite's half-point) but not those three values.  At the measured *R* = 3.7 the gain is
   0.951–1.172 across *K_G*, as the text says.  Over four decades of *K_G* and *R* the
   re-implementation gives gains from 0.66 to 1.57 with only 8 % of the grid below one, so
   "usually below it" is not reproduced; the conclusion the text draws (closing 1.5–1.7 needs three
   or four cooperative sites) is unaffected, since a single site never reaches 1.6.
5. **Sweep of Bhat, Ghat and *A*.**  With the design that produced the text's change-on-doubling
   range (Bhat = Ghat ∈ {0.2, 2, 20}, *A* ∈ {1, 29}, *T_S* = 0.5, *T_E* = 0.05 → 0.10) the change on
   doubling is 0.0008–0.079, as stated, but the coefficient runs from 0.84 (after doubling, *n* = 1,
   Bhat = 20) to 0.999, not from 0.92; and *A* spans 1.5 decades in that design, not three.  With
   *A* over three decades the coefficient runs from 0.82 to 0.999 and the change from 0.001 to 0.09.
6. **The 4,000-replicate digitisation bootstrap (2.01 [1.71, 2.41]).**  This interval was produced
   by `sim/digit_err.py` of the working folder (Gaussian errors, seed 7, 4,000 replicates); it is not
   part of this repository's pipeline, which carries the 600-replicate uniform bootstrap of the
   cascade accounting (536 usable; *h* = 2.03 [1.87, 2.25]).  `fit_titration.py` could be extended
   with the Gaussian bootstrap if the interval is retained in the text.
