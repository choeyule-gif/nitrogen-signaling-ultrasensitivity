# GlnD catalytic states and the six-state GlnE regulator

Run `python code/analyze.py` with NumPy, SciPy and Matplotlib. The authoring and
packaging scripts are workspace helpers and are not required for reproduction.

This extension develops two mechanistic directions without postulating a
permanently refractory PII population. It uses the same twelve supplied digitized
PII coordinates. It adds no new experimental observations.

## GlnD

A three-parameter reference has two identical independent glutamine sites and
state-dependent catalytic activities. A binding/catalytic transformation changes
K1 to lambda*K1, K2 to K2/lambda, and both middle-state specificities by lambda.
The entire steady PII distribution is invariant, while GlnD occupancy and reduced
relaxation rates differ. Lambda=1,2,4 all preserve monotone UT/UR regulation.

The inherited UR-specificity fold is an explicit cross-assay scenario, not a
microscopic parameter measurement. The states are not assigned to individual ACT
domains. The exact free-state identity admits the earlier enzyme-complex
realization; the time courses additionally require rapid regulatory binding and
unsaturated conversion. Model time is not minutes.

## GlnE

The topology and activity assignments are transcribed from Figure 12 of
Jiang, Mayo & Ninfa (2007), DOI 10.1021/bi0620510:

    Z = 1 + g + p + u + p*g/alpha1 + p*u/alpha2
    AT = kP * (p + beta*g + beta_prime*p*g/alpha1) / Z
    AR = kR * u / Z

The steady drive depends on beta_prime/alpha1 and is independent of alpha2.
Paired scaling of alpha1 and beta_prime preserves the drive, yet changes the
PII-binding halfpoint and catalytic speed. Changing alpha2 supplies another
common-rate transformation. The direct PII-independent AT term gives a drive
affine in inverse free-target concentration at fixed glutamine and target state.

Parameter provenance is explicit in results/summary.json. KG, KP and alpha1 are
a source-reported fitting example, KU is a selected apparent-scale scenario,
and beta, beta_prime and the reference alpha2 are illustrative choices. No joint
fit of all historical GlnE assays or cascade coefficients is claimed. Binding
halfpoints predicted here are distinct from blindly treating every apparent
activation constant as a microscopic dissociation constant.

## Reproducible outputs

- `results/glnd_state_families.csv`: identical target means, enzyme occupancies,
  and monotone catalytic capacities for three finite models.
- `results/glnd_occupancy_design.csv`: mean and single-state ligand occupancy at
  K/3, K, and 3K. Mean ligand occupancy alone is uninformative at the central K.
- `results/glne_cascade.csv`: identical steady GS curves and separate AT/AR rates.
- `results/glne_PII_binding_isotherms.csv`: discriminating binding titrations.
- `results/glne_purified_predictions.csv`: occupancy and catalytic-speed changes
  along a glutamine titration with controlled free PII/PII-UMP inputs.
- `results/concentration_signature.csv`: affine GS-odds scenario; the intercept
  subtraction is only for display, not new experimental data.
- `results/local_sensitivity.csv`: analytic route decomposition checked by finite
  differences. The direct partial elasticity bound is not a Hill-coefficient bound.
- `results/step_responses.csv`: joint dynamical illustrations in reference units.
- `results/summary.json`: coefficients, assumptions and verification errors.
- `figures/`: vector PDFs and 600-dpi PNGs for the two main figures.

Verification includes 100 independent PII stationary linear solves, 250 six-state
receptor solves, 36 coupled mechanism/target-total combinations, positivity and
monotonicity checks, and analytic-versus-numerical derivative checks.

This extension is included in v1.1.0.
