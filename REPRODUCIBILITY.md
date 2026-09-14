# Reproducibility of the working revision

## Verified clean regeneration

The v1.2 numerical source snapshot completed all 16 regeneration and check stages in 745 seconds, including original-workbook extraction and count-state bounds. See `validation/clean_regeneration_v1.2.json`, `clean_environment_v1.2.txt` and `clean_source_v1.2.json`. The final rendering and a baseline-export docstring changed after this numerical snapshot; those differences are listed explicitly. This statement covers the earlier v1.2 numerical layer; later calibration and observation-design additions are verified separately.


A fresh virtual environment and source copy with all generated results removed regenerated every pre-existing numerical layer. `validation/clean_regeneration_pre_paralog.json` records all 12 commands, durations and zero exit codes; `validation/clean_environment.txt` records dependencies. The baseline passed 22 checks with no skips. Matched-model verification used 123 independent stationary solves, 300 coordinate draws and cascade comparisons. Three extension test groups check regulatory equivalences, conservation/positivity and independent integration, and ladder/design bounds.

This evidence covers numerical regeneration, not visual reproduction of final publication layouts, global optimization, experimental uncertainty coverage or native mechanism validation. The subsequent paired-paralog addition is checked separately and by the current integrated entry point. The saved clean report intentionally identifies its earlier scope.

## Commands

```sh
pip install -e .
python tests/run_research_suite.py --regenerate
```

Without the flag, saved baseline values are checked and extension identities are solved afresh. `results/research_suite.json` is successful only when every scheduled command completed and passed. `run_all.py` alone reproduces the baseline; it does not run the later research layers.

To re-extract the unchanged independent workbook:

```sh
pip install openpyxl
python extras/paired_paralogs/code/extract_data.py
python extras/paired_paralogs/code/state_bounds.py
python extras/paired_paralogs/code/paired_reporter.py
python extras/paired_paralogs/code/reporter_intervals.py
```

The source DOI is 10.6084/m9.figshare.4880003. Cell-level provenance, mean/SE/n and the original files are retained. No replicate-level covariance is supplied. Peptide ratios must not be mistaken for intact-trimer distributions. Reporter interval checks cover 10,000 admissible synthetic configurations; synthetic examples are not observed sequestration. Dynamic fitting scripts are exploratory and are not evidence that the native network has been validated.

## Manuscript and packaging

The working main manuscript has five main figures and six supporting figures; the mapping is in README.md. S2 Appendix covers data and sensitivity; extended legacy analyses are preserved in `docs/extended_analyses/v1.2.0/`. A synchronized working submission ZIP is built separately from the released archive. Do not use the old release package as the revised manuscript. Compile the current LaTeX sources twice with pdfLaTeX; upload main figure files separately according to the journal's submission workflow. Archive v1.1.0 has DOI 10.5281/zenodo.22735182 and does not contain the new paired-paralog extension.

## Revision 4 / v1.3.0 observation-design verification

The new analyses were generated from the supplied revision3 S1 Data baseline; matching files had no content differences from the baseline repository. The recorded 75 main and twelve local-effect configurations comprise 8,700 Monte Carlo trials, 69,600 complete experiments and 4,593,600 synthetic technical readings. One experiment contains three curves, eleven doses and two readings per dose. The noise magnitudes are declared scenarios.

`extras/observation_design/results/revision_verification.json` independently recounts all trial decisions, checks the null envelope on 3,000 parameter points, checks binding occupancy against 500 independently solved ligand mass balances, verifies exact curve recovery and 500 covariance identities, compares the parameterized reaction model with baseline defaults, and solves a fresh source-rate reaction pair. The separate certificate report checks 120 paired finite realizations, 5,000 positive ladders and 5,000 free/bound mixtures. The previous clean-environment reports are retained with their original scope, not relabeled as a clean regeneration of the new Monte Carlo study.

Regenerate and check this layer using the commands in `extras/observation_design/README.md`. All dependencies, source hashes and the final integrated validation status accompany the current revision. These are computational checks of declared scenarios; matched biological calibration and direct occupancy validation remain outstanding.


## Verification patch v1.3.1

The **downloaded public v1.3.0 submission ZIP**, not the local working outputs, was unpacked into a new source copy. The nine generated-output directories listed in `validation/public_source_audit_v1.3.1.json` were removed. A new isolated virtual environment installed that package with `pip install -e .` and completed all **24 scheduled regeneration/check stages in 810 seconds**, with zero failures. The original 22 baseline checks had zero skips. This includes both observation benchmarks: 8,700 trials, 69,600 complete experiments and 4,593,600 synthetic technical readings. The nine central observation-design result files were **byte-identical** to the public release copies; see `regeneration_comparison_v1.3.0.json`. Legacy numerical searches are checked by their stated tolerances rather than a byte-identity requirement. Optional exploratory scripts outside the scheduled suite were not all rerun.

The recorded source is v1.3.0, commit `39774c720ef713fc20a007748e74ace131cc80af`. The v1.3.1 patch changes two rounded sufficient thresholds, metadata, documentation and validation code; it does not change the numerical models or figures. `tests/test_observation_independent.py` was additionally run against the clean regenerated outputs. It audits 65-digit guarantees, 15 count-gate linear programs, nine limiting ladders, two independent source-rate BDF integrations, an alternative raw-data refit and all 8,700 trial decisions. The current full regeneration command therefore schedules **25 stages**; the historic 24-stage report is preserved with its actual source version.

All ten documents in the public portable LaTeX ZIP compiled twice without undefined references or overfull boxes. Reference audit covers the titles and DOIs of all 36 bibliography entries; the biochemical audit also directly checked the 1998 catalytic table, 2007 main text/S10 coefficients and the 2011 assay conditions. DOI resolution does not validate every cited scientific inference. The 2007 K_P unit discrepancy is present in the original sources and remains explicitly included in sensitivity analysis.

Environment, public asset hashes, source identity and machine-readable results are in `validation/`. Direct integrations verify two selected starting conditions, not global convergence. All assay-design guarantees remain conditional on the manuscript's topology, free-input, observation-gate and coefficient assumptions.
