# Protein-conserving cascade extension

This extension is included in v1.1.0 and accompanies the revised manuscript. The v1.0.0 archive covers the earlier baseline only.

## Reproduce the numerical research

Python 3.10+ with NumPy, SciPy, and Matplotlib is sufficient. From this directory:

```sh
python3 code/analyze.py
python3 code/topology.py
python3 code/total_input.py
python3 code/verify_algebra.py
python3 code/plot.py
```

The scripts use paths relative to their own files. Results are written to `results/`; The source plot for final Figure 10 is written to `figures/`. The manuscript-building scripts are editorial workspace utilities that require the previous manuscript directory; they are not required to reproduce the numerical results. The delivered LaTeX files can be compiled independently with pdfLaTeX twice using the included figure files and style files.

## Scientific scope

- The base model has 120 species and 300 directed reactions, including every count-indexed free PII and GS modification state and explicit catalytic complexes. Four protein totals are conserved exactly.
- Catalytic complexes retain their regulatory state until release. Regulatory ligands exchange only on enzymes free of catalytic target. This is a sufficient architecture, not a complete physiological nitrogen model.
- Glutamine is a buffered free input in the principal calculations. Nucleotide donors/products and metabolic regulation are not dynamic species.
- A separate steady-state extension conserves regulatory glutamine binding. It has no independent glutamine-dynamic ODE validation.
- Partly modified PII recognition, catalytic association coefficients, regulatory rate scales, GS interaction strengths, and pool sizes are illustrative model choices. No new experimental fit was performed.
- Fractions describe UMP or adenylyl modification, not measured GS activity. The user's original Figure 1 and all prior figures are preserved.
- The GlnD fitted scales were converted from mM to micromolar. The borrowed specificity numbers set arbitrary rate units; simulation time is not calibrated biological time.

## Results and limits

`summary.json` records 8,200 paired buffered-input comparisons, 60 randomized steady realizations, and three full mass-action integrations. The finite-state conservation identities and the cascade fraction bounds hold in these checks. `topology_summary.json` records the productive EPU boundary and the separate compensation bound that survives that extension. `total_input_summary.json` records the distinction between fixed free and fixed total glutamine. `algebra_verification.json` checks the hidden-conformer example and arbitrary positive twelve-state-ladder bounds (13 count states).

The protein-pool inequality is an upper bound on differences, not a prediction that the largest difference will occur. Its inversion provides a sufficient experimental condition for controlling interpretation error. It gives neither a slope bound nor a guaranteed discrimination budget.

The productive EPU construction is an alternative reaction model, not evidence that native GlnE has that activity. Splitting a regulatory state into active and inactive conformations demonstrates why target plus ligand occupancy cannot identify an unrestricted hidden state space. The proposed binding designs therefore resolve the stated model families, not every possible biochemical mechanism.

## Source provenance

The GlnD fit and GlnE parameter distinctions are inherited from the accompanying regulatory-state study. The GlnE class topology is motivated by Jiang, Mayo and Ninfa (2007), DOI 10.1021/bi0620510. The new state-resolved refinement is explicitly a modeling choice. Prior scaling/identifiability theory is acknowledged through Castro and de Boer (2020), DOI 10.1371/journal.pcbi.1008248, and Villaverde and Massonis (2021), DOI 10.1371/journal.pcbi.1009032. No priority claim is made for scaling symmetries, general non-identifiability, or the variance identity.
