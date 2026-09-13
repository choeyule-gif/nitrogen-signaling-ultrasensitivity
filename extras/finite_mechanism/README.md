# Finite biochemical realizations and constraint tests

This extension addresses two issues: whether a population nonidentifiability
requires an arbitrary switch function, and whether biochemical measurements
actually constrain the alternative mechanisms. It contains no new experiments.

Run `python run.py` with NumPy, SciPy, and Matplotlib installed. The manuscript
authoring script is workspace-specific and is not needed for reproduction.

## Models and observation assumptions

Three GlnD ligand states with positive sequential binding weights and
state-dependent catalytic specificities produce a rational mean. A positive
coefficient transformation creates an exactly mean-matched refractory population.
The main numerical construction has two identical independent ligand-binding
sites, monotone UT/UR activities, three fitted parameters, and one UR specificity
fold transported from older enzyme assays as an explicit scenario.

The effective rational coefficient ratio is NOT the microscopic binding ratio.
The main fitted model has microscopic K1/K2=1/4, while a different realization
that transports two kinetic summaries requires strong binding cooperativity.
The primary model does not fit all the older kinetic measurements. A stricter
cross-assay transfer deliberately reports a poor fit rather than discarding it.

A separate-site GlnE saturation scenario can break gain-only equivalence at fixed
affinities. Rescaling the unknown unmodified-PII affinity restores the receptor
occupancies with one shared transformation over all inputs and target totals.
Reverse recognition is assumed proportional to modification count. The four
occupancy classes represent eight underlying receptor microstates. This is not
a measured state-recognition law or a full closed-total cascade reconstruction.

## Inputs and outputs

- `data/jiangninfa2011_fig2B_digitised.csv`: same 12 author-supplied digitized
  coordinates used by the earlier analysis; not raw replicate measurements.
- `data/primary_source_constraints.csv`: selected, checked source-table entries,
  assay contexts, and limits on their use. Apparent constants are not microscopic Kd.
- `results/summary.json`: fitted parameters, success flags, identity errors,
  independent-binding and cooperative scenarios, and failed strict transfer.
- `results/finite_pair_curves.csv`: complete generated state/rate predictions.
- `results/cooperativity_profile.csv`: effective coefficient-cap scenarios.
- `results/transported_kinetics_profile.csv`: two-summary transfers at fixed ratios.
- `results/strict_cross_assay_transfer.csv`: failed stricter-transfer predictions.
- `results/separate_site_saturation_scenarios.csv`: 663 scenario combinations.
- `results/network_verification.json`: every species balance for 10 explicitly
  constructed finite-network steady states and 1,000 random coefficient checks.
- `figures/Fig7_finite_mechanisms.pdf`: vector figure, plus a 600-dpi PNG companion.

## Interpretation

Finite realizability is established for the stated free-input architecture.
Physiological validity is not established. The fit's baseline UT/UR ratio and
individual activity half-ranges differ from older enzyme assays. Apparent
binding-activation constants are context-dependent and are used only as scenario
envelopes. No experimental likelihood, acceptance probability, or detection power
is inferred from these calculations.

The explicit enzyme-complex verification checks free-input steady states, whose
computed conserved totals can differ across inputs. It is not a fit to a titration
at common closed-system totals. New analyses are local revisions and are not part
of the earlier v1.0.0 Zenodo archive.
