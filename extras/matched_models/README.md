# Matched-model research revision (13 September 2026)

This addition investigates what remains identifiable after different mechanisms
explain the same existing PII mean titration. No new experimental observations
are generated. It supplements, and does not silently replace, the version 1.0.0
baseline analysis archived at https://doi.org/10.5281/zenodo.22731769.

Run `python3 run.py` from this folder (or by absolute path). Requirements:
Python 3.10+, numpy, scipy, matplotlib. The recorded execution used Python 3.14,
numpy 2.5.2 and scipy 1.18.1. All paths for computational reproduction are relative
to this folder. The two `revise_*` scripts are workspace-specific manuscript
builders and are excluded from the portable analysis distribution.

## Inputs and provenance

- `data/jiangninfa2011_fig2B_digitised.csv`: the same twelve author-supplied,
  approximate digitized coordinates analyzed in the baseline manuscript.
- `data/cascade_comparison.csv`: only its first three columns are observational
  inputs here (PII concentration, reported midpoint, reported coefficient).
  Remaining columns are old baseline calculations, not additional observations.
- Original experimental source: Jiang P, Ninfa AJ (2011), Biochemistry,
  https://doi.org/10.1021/bi201410x. The source-data limitations and original
  broken-axis endpoint issue remain documented in the manuscript.

## Analyses

1. Three four-parameter finite-switch fits at fixed ladder scenarios, a 61-point
   ladder profile, and leave-one-coordinate-out prediction checks.
2. Exact mean matching for positive ladders by inverse construction, and a
   distinct fully modified refractory population with the same mean.
3. An exact gain-rescaling equivalence that preserves a downstream response
   based on the unmodified-state/mean-uridylylation ratio.
4. Sixteen fixed-upstream cascade comparisons; four additional shared-upstream
   fits; one asymmetric GS ladder comparison. Added parameters are shared across
   all three assay conditions, with three midpoint-calibration nuisance gains.
5. Three hundred shared coordinate perturbations; separate, deterministic
   +/-10% summary-tolerance scenarios; candidate-specific state-measurement design.
6. Independent root, derivative, conservation, and stationary-generator checks.

## Interpretation boundaries

Finite-switch fits and inverse-constructed mean-matched models are different
analyses. The latter do not identify a finite ligand-binding implementation and
are not assigned an AIC. The transfer-power model is phenomenological. Effective
ladder constants may encode state-specific recognition, not only intrinsic site
interactions. The refractory population is hypothetical.

Cascade data consist of only three midpoint/coefficient pairs, not replicated
curves. No synthetic points reconstructed from those summaries are treated as
measurements. No experimental likelihood or model posterior is claimed.
Boundary fits and very large calibrated gains are reported rather than treated
as measured kinetic constants. Agreement of two optimizers is not proof of a
global optimum. Coordinate percentiles and summary tolerances are not assay
confidence intervals. State-design predictions currently use an effective free
target interpretation; total-state observations may require binding corrections.

The robust findings are conditional mathematical invariants and observational
equivalences. A universal low-PII amplification deficit or 26-fold transferase
inhibition requirement is **not** claimed. In particular, a fully modified
refractory population does not itself prove weak regulation in accessible target.

## Files

`results/summary.json` records parameters and model comparisons; CSVs give every
curve, fit profile, held-out prediction, perturbation, and design score.
`results/verification.json` and `results/final_checks.json` record scientific
checks. `figures/` contains generated vector PDF, 600-dpi TIFF, and preview PNG
files for the revised main figures. The original baseline and original design
figures are retained as supplementary figures in the revised manuscript.
