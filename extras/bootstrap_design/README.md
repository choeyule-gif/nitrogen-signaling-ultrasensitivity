# Expanded raw-assay calibration benchmark

All numerical work for this layer was run on a dedicated Linux computation server. The numerical environment is recorded in each
results directory. The original 100-trial analysis remains archived separately.

## Design and interpretation

Two fitting strategies are compared using matching random seeds and hence
matching synthetic raw observations: shared background/gain (five parameters),
and separate background/gain for each glutamine condition (nine parameters).
Each trial contains eight independently generated complete titration sets,
each with three glutamine conditions, eleven protein doses and two readings.

Each fitting strategy uses 3,000 raw-simulation/refit trials at each of 27 null
settings, followed by independent evaluation with 2,000 trials at each of 210
settings. These give 4,008,000 nonlinear fits per strategy, 8,016,000 in total.
The raw datasets are paired between strategies; these are not 1,002,000
independently generated datasets. Calibration and evaluation phases have
different seeds. The trial files retain means, standard errors, failure and
boundary indicators, residual error and fitting effort.

The null reference is alpha=0.17. The 27 null settings cross three KG values,
three stock/enzyme error scales and three background/gain noise scales. The
cutoff is the maximum of per-setting parametric null-bootstrap order statistics.
Each simulated null trial includes generation of the raw measurements and
nonlinear re-estimation of its halfpoints. A binomial tolerance construction
accounts for calibration Monte Carlo uncertainty, with Bonferroni adjustment
across the 27 null settings. This is fixed-null grid calibration, not a
data-dependent plug-in bootstrap. Its guarantee applies to the finite grid;
it does not prove a supremum over the intervening continuous parameter region.

The independent audit has four components:

- 108 settings: three KG values × three concentration-error scales × three
  signal-error scales × four regulatory coefficients.
- 33 settings: three KG values × eleven continuous-path mixture strengths.
- 42 settings: three KG values × seven dose biases × with/without independent
  concentration correction.
- 27 settings: three KG values × nine differential background or gain biases.

All errors are prospective scenarios. They are not estimates of experimental
uncertainty. Separate background/gain fits correct an affine observation
distortion under the modeled signal relationship; they do not correct hidden
affinity mixtures or identify the native molecular mechanism.

## Verified results

The calibrated cutoffs are 2.06252095 (shared fit) and 2.10355999
(condition-specific fit). Independent null rejection rates are at most 0.0415
and 0.0385 over their respective 27 grid settings. A maximum estimated rate
below 0.05 is not an empirical proof of uniform frequentist coverage.

At KG=15.6 mM, a differential background of -0.02 gives rejection 0.924 with
the shared fit and 0.027 with the condition-specific fit. A log-gain shift of
-0.05 gives 0.992 versus 0.0295. The same raw-noise realizations are used for
each comparison. With the condition-specific fit, a mixture strength of 0.2
still gives rejection 0.567; at strength 0.5 it reaches 1.0. Rejecting the
single-class reference in that case does not identify a change in alpha.

No optimizer failures occurred in the saved evaluation runs. Active parameter
bounds occurred in 17 and 22 fits, respectively; their trials remain in the
denominator and do not reject. They are not discarded as if successful.
The verification script checks both indicators separately, all 840,000 saved
evaluation decisions, and complete trial replay from seeds. Fifty alternative
fits use a free-ligand root, a different starting point and a different
optimization method. Maximum score disagreement is 1.38e-6.

## Reproduce

From the repository root, using the installed project environment:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python extras/bootstrap_design/code/bootstrap_binding.py --trials 2000 --calibration-trials 3000 --workers 16
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python extras/bootstrap_design/code/bootstrap_binding.py --trials 2000 --calibration-trials 3000 --workers 16 --condition-affine --output results_condition_affine
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python extras/bootstrap_design/code/verify_bootstrap.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python extras/bootstrap_design/code/plot_bootstrap.py
```

Complete cell files are resumable. Use a new output directory to perform a
fresh Monte Carlo regeneration; existing compatible cells are reused.

## Figure caption for manuscript integration

**Calibrated raw-assay simulations separate correctable signal errors from
mechanistic ambiguity.** (A) Independently evaluated null rejection at 27
settings per fit strategy. Each point represents 2,000 trials; nine error
settings occur at each KG. Both cutoffs were calibrated from separate raw
null simulations including nonlinear fitting. (B) Condition-specific-fit
rejection versus alpha; lines show means and shading the range across nine
error settings at each KG. This range is a nuisance-scenario range, not a
confidence band. (C) A condition-dependent signal background can cause strong
rejection under the shared-fit observation model; fitting condition-specific
background/gain removes the tested affine distortion. (D) A continuously
increasing hidden affinity mixture still changes the single-class inference.
Panels C–D fix KG=15.6 mM and show 95% Wilson Monte Carlo intervals. The
dashed reference is 0.05. Every trial contains eight complete titrations, not
eight dose points. Signal and concentration uncertainties are stipulated
scenarios; this figure does not validate a native assay or unrestricted model.
