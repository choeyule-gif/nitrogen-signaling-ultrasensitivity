# Assumption boundaries and nuisance-aware experiment design

Run `python code/analyze.py` with NumPy, SciPy and Matplotlib. The portable
`regulatory_core.py` carries the previous GlnD reference model; its `main()` is
not invoked. No original experimental data are added. The same twelve digitized
coordinates determine the reference coefficients. All wider parameter and error
ranges below are explicit scenarios, not fitted physiological confidence intervals.

## Exact statements

1. With ligand exchange on free enzyme, three enzyme states, six catalytic
   complexes and two target states give an 11-concentration mass-action system.
   Complex balances recover the free-state ratio for arbitrary positive finite
   ligand-binding speed and arbitrary substrate saturation. The model uses one
   target modification unit per complex, not an entire resolved PII trimer.
2. With common free modified fraction and common totals, the difference of two
   total modified fractions is at most `min(1, enzyme_total/target_total)`.
   Symmetric forward/reverse effective affinities make total and free fractions
   identical even at saturation. The bound does not bound derivatives.
3. Ligand exchange on target-bound enzyme preserves the construction when its
   generator balances the catalytic complex weights. Binding-cycle consistency
   alone does not ensure this: state-dependent dissociation/catalysis ratios can
   generate regulatory currents and break the identity.
4. Target-free singly liganded GlnD at nominal K has independently profiled ranges
   [4/9,1/2] for lambda=1 and [1/6,1/5] for lambda=4, if each true K lies between
   K_nominal/2 and 2*K_nominal. Their gap is 11/45. For a general alternative,
   scalar ranges overlap when lambda<=1.25; more repeats alone cannot resolve it.
5. With PII-UMP absent, normalizing a GlnE PII-binding halfpoint by its value at
   zero glutamine removes KP. At half occupancy, total PII equals the free
   halfpoint plus half the GlnE total. Correct that depletion before taking ratios.
   Titrations at glutamine=0,1,10 mM identify alpha1 and KG ideally, except for
   alpha1=1; noisy data require fitting the original binding model.

## Numerical model and parameter provenance

All GlnD coefficients and reference K are inherited from the same fit, with a
shared conversion-time normalization. In the asymmetric saturation illustration,
`hplus=(.2,2,10)`, `hminus=(5,.3,.1)` in inverse target-concentration units and
`target_total=1`; these are chosen scenarios. `kcat=a/h`, `koff=kcat`,
`kon=2*h*kcat` realize the stated effective specificity. Ratios E/T range from
1e-4 to 1. The random verification independently varies h entries from .01 to100,
T from .1 to10, and ligand speed from .001 to1000.

The substrate-bound exchange extension changes `koff/kcat` to (.1,1,10) for UT
and (10,1,.1) for UR, retaining h and a/b. Then `kon=(1+rho)*h*kcat`.
Complex ligand transitions use equilibrium association weights
`h*(1+rho)/rho`, enforcing noncatalytic binding-cycle ratios. Five G/K values
(.3,1,3,10,30), 26 exchange speeds (zero and 1e-4 to100), and both lambdas are
solved positively. Maximum free-fraction gaps near eta=.1 are .0349 at10K and
.0703 at30K. These are mechanism sensitivity examples, not biological estimates.

The GlnE topology and nominal parameters come from the preceding extension
(Jiang, Mayo & Ninfa 2007, DOI 10.1021/bi0620510). Here the two disjoint test
classes are alpha1 in [.136,.2125] versus [.544,.850], with independent KG in
[7.8,31.2] mM; KP cancels. These are specified effect-size classes, not measured
parameter ranges. Broader classes can overlap and invalidate the guarantee.

## Precision and replication

`noise_design.csv` labels the observation and error units. GlnD occupancy SD is
an absolute fraction SD; GlnE halfpoint SD is a log SD. Every halfpoint estimate
requires a complete independent titration. The one-sided normal test uses known
pilot SD, size .05 PER assay, and target power .80. The least-favourable null
boundary and nearest alternative define the guarantee over the declared classes.
At SD=.20 this gives five occupancy estimates or eight halfpoint estimates at
each of two glutamine levels. It does not mean eight individual wells. Shared
calibration factors may cancel; unknown dose-specific bias does not. An assay
with estimated variance, correlation, biological variability, or clipped fraction
errors requires a new power analysis. Independent replication and synthetic
noise draws are not new biological data.

Kinetic contrasts are separately marked CONDITIONAL. Paired capacities cancel a
shared speed, but an unknown baseline K and catalytic shape can restore overlap
when refitted independently. The capacity grid optimum (.020K,.855K) is therefore
not used for the nuisance-profiled guarantees in main Fig10. Reciprocal doses
z and 1/z are uninformative for the symmetry multiplier.

## Files and verification

- `total_observation_boundary.csv`: total-readout differences and rigorous bound.
- `finite_binding_dynamics.csv`, `rapid_binding_error.csv`: finite regulatory lag.
- `substrate_bound_ligand_exchange.csv`: topology extension and full residuals.
- `noise_design.csv`: effect-size-specific sample counts and power.
- `continuous_alternative_resolution.csv`: unresolved overlapping classes are Inf.
- `GlnE_depletion_correction.csv`: exact total-to-free halfpoint correction.
- `summary.json`: scope, parameters, errors and 100,000-draw Monte Carlo checks.
- `figures/`: vector PDF and 600-dpi PNG, four panels each for main Fig9 and Fig10.

Checks include 300 random mass-action steady states, positivity and conservation,
260 extended-model steady states, independent BDF convergence for both models in
each topology, exact binding inversion, and analytic versus synthetic test power.
The authoring/package helpers are workspace-only and excluded from the portable
bundle. The extension is included in v1.1.0.
