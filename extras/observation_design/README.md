# Observation-conditioned cascade bounds and end-to-end binding design

This layer addresses the reviewer revision 4 methodological requirements. It adds a concentration decision based on an observed PII count, a raw binding-assay benchmark, and a GlnD kinetic-constraint profile. All new experimental settings are prospective scenarios. No new biological measurements, native occupancy validation or experimentally estimated confidence regions are claimed.

## Reproduce

From the repository root:

```sh
python extras/observation_design/code/observation_certificate.py
python extras/observation_design/code/verify_observation_certificate.py
python extras/observation_design/code/constrained_glnd.py
python extras/observation_design/code/glnd_constraint_profile.py
python extras/observation_design/code/binding_observation.py --trials 100
python extras/observation_design/code/binding_local_effects.py
python extras/observation_design/code/verify_revision.py
python extras/observation_design/code/plot_revision.py
```

The integration command `python tests/run_research_suite.py --regenerate` includes these numerical analyses; default mode verifies the saved 8,700 trial summaries, freshly checks binding mass balance and refits noiseless curves, and solves a new full-reaction pair. The larger Monte Carlo regeneration is separate from a default check. Figure rendering is explicit because historical scripts retain old figure numbers.

## Files and assumptions

- `OBSERVATION_CERTIFICATE.md`: exact quantity bounded, proof, source-coefficient decision example, unit sensitivity and measurement gate.
- `SOURCE_AND_ASSAY_AUDIT.md`: primary-source assay conditions, fit versus measured parameters, and practical observable requirements.
- `data/jiangninfa2011_fig2B_digitised.csv`: identical copy of the repository's twelve author-supplied digitized coordinates, DOI 10.1021/bi201410x; not replicate-level observations.
- `model_source.py`: parameterized version of the existing finite-cascade realization. Defaults are checked against the previous model. Source-rate tests hold GlnE sigma at one and vary the GlnD equivalence parameter.
- `binding_observation_trials.csv`: 7,500 trial records; 100 trials at each of 75 declared configurations. Eight full experiments, each containing all three glutamine conditions and 66 technical observations, per trial.
- `binding_local_effects_trials.csv`: 1,200 extra trial records at twelve small-effect corrected-dose settings; independent seed.
- `binding_observation_example.json`: one raw synthetic titration set used in Fig4A.
- `certificate_full_reaction_checks.json`: 120 paired conserved-reaction realizations and their observed-count/generic bounds.
- `glnd_constraint_profile.csv`: fixed-binding-scale stress profile, not a confidence interval. The upper-bound optimum is not a validated affinity.
- `bound_decision.csv`: retained exploratory fixed-GS-capacity sweep using the earlier illustrative rate ratios. It is not the source-anchored main decision example.

The 2007 primary supplement used for the source-fit ratios is DOI 10.1021/bi0620510.s001, Fig S10, PDF page 16. Its catalytic coefficients are 7.5, 175 and 1315. The main-text K_P unit is used and the inconsistent literal supplementary unit is explicitly tested. The s002 supplement is a different file. Copyrighted primary PDFs are not redistributed here; only source-attributed values needed for the analyses are retained.

The bound is on pairwise differences within the specified free-state-equivalent mechanism class. It is not an error relative to biological truth. Source-fit coefficients and the PII count gate require independent justification in the proposed assay. A large upper bound does not establish detectability. The Monte Carlo intervals quantify simulation sampling uncertainty, not native-assay calibration. Rejection under a hidden affinity mixture does not identify a specific change in the single-class regulatory coefficient.
