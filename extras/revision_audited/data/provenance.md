# Data provenance

The canonical twelve-point input, `jiangninfa2011_fig2B_digitised.csv`, was supplied by the author in the nitrogen-papers repository. Its accompanying README attributes it to by-eye digitization of an enlarged Fig. 2B of Jiang and Ninfa (2011), with 0.5 µM PII and 0.05 µM UTase/UR. It is not raw experimental replicate data.

`reconstructed_titration.csv` retains an independent vector-marker extraction from page 19, Fig. 5b of the supplied manuscript PDF. It reproduces the same plotted coordinates with small numerical extraction differences. Use the author-supplied rounded CSV for analysis. SI Eq. S21 records the extraction calibration.

`author_revision_numbers.json` is the supplied numerical summary. `analysis_summary.json` is the independent Python reproduction from the supplied rounded coordinates. Their conclusions and displayed values agree. `published_constants.json` preserves author-supplied literature constants; reported cascade summary values are treated as fixed, not replicated measurements.

The 2,000 perturbation draws assume a log-input SD of 0.06 and ordinate SD of 0.06 UMP/trimer; their percentile intervals are coordinate sensitivity ranges, not experimental confidence intervals. The nominal likelihood-profile cutoff is a conditional descriptive criterion.

`figdata/` and `figures/matlab/` preserve the supplied MATLAB panel data and graphics snapshot. The five assembled figures retain vector graphics and add letters only. In Fig. 2A the blue curve is a comparison contour at product 2.018, not a measured local sensitivity; the historical CSV filename `Fig2A_measured_hyperbola.csv` does not change that interpretation. Seven supplementary figures are now included. S1 and S3 preserve the supplied vector panels after checks; S2 and S4–S7 were rebuilt in Python from the supplied models, data, and explicitly documented corrections. The original MATLAB panels remain archived for comparison.

The full Jiang–Ninfa 2011 primary paper was checked locally. Its coefficients use 10–90% of the full response. The Fig. 2B upstream trace was measured within the bicyclic system, at PII 0.5 µM, UTase/UR 0.05 µM, and ATase 0.1 µM. The three cascade conditions vary both PII (0.5, 5, 36 µM) and UTase/UR (0.05, 0.5, 0.8 µM). A shared upstream fit across these conditions is a modeling assumption.

The supplied CSV first coordinate, x=0.02 mM, represents a plotting proxy for the zero-input point before the source's broken logarithmic axis. The canonical CSV is preserved for reproducibility. `audit/endpoint_checks.json` reports refits with zero input, zero input and full ordinate, and the first point removed. These are sensitivity analyses, not new experimental measurements.

`audit/*_author.json` files preserve author inputs or source-script outputs; `audit/omission_checks.json`, `replacement_condition_i_checks.json`, `exact_root_count.json`, `exact_physical_roots.json`, `ref_checks.json`, `endpoint_checks.json`, `locus_checks.json`, and `figure_checks.json` record the additional audit. The new condition-(i) counterexample is a mathematical construction with large kinetic heterogeneity, not a fitted physiological parameter set. `audit/new_condition_i.json` is an intermediate search result superseded by `replacement_condition_i.json`.

See `../audit_report_ko.md` for the complete figure crosswalk and limits of the audit. Historical symbolic certificates and all original ODE bistability claims were not independently re-certified in this revision.
