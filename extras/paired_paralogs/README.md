# Paired GlnB/GlnK measurements and conditional sequestration inference

This extension adds independent published measurements, count-state upper bounds on binding, and a prospective calibrated assay construction. It does not establish the native GlnD mechanism or infer measured AmtB occupancy from the available time series.

## Source and attribution

Gosztolai et al. (2017), *GlnK facilitates the dynamic regulation of bacterial nitrogen assimilation*, Biophysical Journal, doi:10.1016/j.bpj.2017.04.012. The original workbook and model are from doi:10.6084/m9.figshare.4880003, distributed under CC BY 4.0. Their source metadata, checksums and unchanged files are in `sources/`. This third-party license is distinct from the repository code license. `results/observations.csv` retains source cells, means, SE and n. Modified-peptide measurements do not resolve intact trimer states.

## Reproduction

Run from the repository root with numpy, scipy, matplotlib and openpyxl installed:

```sh
python extras/paired_paralogs/code/extract_data.py
python extras/paired_paralogs/code/state_bounds.py
python extras/paired_paralogs/code/paired_reporter.py
python extras/paired_paralogs/code/reporter_intervals.py
```

The first program re-extracts the original workbook. The remaining programs check sharp moment bounds, solve 2,000 independent count-state generators, and verify interval coverage for 10,000 admissible synthetic configurations. These are mathematical checks, not experimental coverage guarantees.

`PAIRED_REPORTER_RESULT.md` states the assay assumptions and calibration requirement. `reporter_intervals.py` propagates interval-valued control and assay measurements, and explicitly permits bounded error in the free-pool relation. The synthetic design examples are not biological estimates.

## Exploratory dynamics

`native_trimers.py` and `profile_paralog.py` fit a positive count-state model with buffered AmtB, imposed protein totals and assumed unmodified synthesis. `verify_native.py` compares independent BDF and LSODA integrations at saved parameters. The fits do not establish predictive biological validation: rate estimates reach a search boundary and mutant comparisons remain imperfect. `dynamic_pilot.py` is a diagnostic linear-interpolation port of the older model, not a faithful reproduction claim. Its scores are not used as biological evidence. The prospective stationary assay does not assume the published time points are equilibrated.

Main Figure 5 is rendered by `python extras/paired_paralogs/code/plot_paired_reporter.py` after data extraction and state-bound calculation. It uses 7.5-inch width, Arial labels and 600-dpi raster export. Panel D propagates the observation-model alternatives into binding ceilings; it is not an occupancy estimate.
