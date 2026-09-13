# Reproducibility of the working revision

## Verified clean regeneration

The v1.2 numerical source snapshot completed all 16 regeneration and check stages in 745 seconds, including original-workbook extraction and count-state bounds. See `validation/clean_regeneration_v1.2.json`, `clean_environment_v1.2.txt` and `clean_source_v1.2.json`. The final rendering and a baseline-export docstring changed after this numerical snapshot; those differences are listed explicitly. All numerical analysis source files match the verified snapshot.


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
