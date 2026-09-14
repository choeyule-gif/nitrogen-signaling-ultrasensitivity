# Observation-conditioned bounds define experimental decision criteria for bifunctional enzyme systems

Code, published measurement summaries, numerical results and manuscript sources for Yule Choi's study. **Source version: v1.4.0.**

The contribution is an observation-conditioned enclosure for a declared biochemical model class. It derives sufficient conditions for decisions despite microscopic ambiguity, extends beyond enzyme-concentration cancellation to AceK–IDH, and checks how raw binding measurements support or confound discrimination. The nitrogen and AceK concentration gates are prospective requirements, not already verified native-assay decisions. Source-summary calibration, independent-activity checks and synthetic performance have distinct evidential roles.

The version-specific source is [GitHub release v1.4.0](https://github.com/choeyule-gif/nitrogen-signaling-ultrasensitivity/releases/tag/v1.4.0). The [Zenodo concept DOI 10.5281/zenodo.22731768](https://doi.org/10.5281/zenodo.22731768) indexes the versioned archive; select **1.4.0** for this manuscript. The earlier v1.3.4 DOI 10.5281/zenodo.22741777 does not contain these new analyses. The final submission manifest records the independently verified archive identifier and file hashes.

## Reproduce

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
python tests/run_research_suite.py
```

The current default suite has **13 stages**. It checks saved baseline outputs, independently recounts the stored observation trials, and freshly solves extension identities, reaction models, general graph enclosures and withheld-activity checks. It does not rerun all eight million nonlinear fits. `--regenerate` additionally recreates the earlier numerical layers; the new large Monte Carlo generator has separate commands in `extras/bootstrap_design/README.md`. Existing compatible Monte Carlo cells are resumable; use fresh output directories for clean regeneration. See `REPRODUCIBILITY.md` for the exact scope of each verification record.

Build the complete S1 Data archive with `python scripts/package_research_data.py`. The older `scripts/make_s1_data.py` produces only a legacy baseline export.

## Research layout

| Location | Content |
|---|---|
| `esbm/`, `scripts/`, `data/`, `results/`, `figdata/` | Baseline theory, digitized titration and plot inputs |
| `extras/matched_models/` | Matched mean responses, distinct ladders and downstream calibration |
| `extras/regulatory_states/`, `extras/finite_mechanism/` | Finite GlnD/GlnE mechanisms and regulatory equivalences |
| `extras/closed_cascade/` | Conserved cascades, topology and free/total-input boundaries |
| `extras/observation_design/` | Count-conditioned bounds, GlnD profiles and earlier 8,700-trial benchmark |
| `extras/bootstrap_design/` | Expanded raw-assay calibration: 27 × 3,000 calibration and 210 × 2,000 evaluation trials per fit strategy |
| `extras/bifunctional_generalization/` | General graph/LP algorithm, full AceK dimer network, six source summaries and three withheld ATPase checks |
| `extras/robust_design/` | Earlier nuisance-aware binding and occupancy designs |
| `extras/paired_paralogs/` | Published GlnB/GlnK data and prospective calibrated reporter extension |
| `extras/revision_audited/`, `matlab/` | Legacy figure calculations and rendering inputs |
| `manuscript/` | Current LaTeX, figures, appendices and cover letter |
| `validation/`, `tests/` | Version-scoped evidence and executable checks |

## Figure and supplement mapping

The five main figures are: (1) the author's unchanged pathway diagram; (2) GlnD equivalence; (3) the conditional nitrogen concentration criterion; (4) expanded raw-assay calibration; (5) generalization to AceK. The previous paired-paralog Fig. 5 is now **S7 Fig**, with its interpretation preserved in S3 Appendix E7. **S8 Fig** retains raw-titration and GlnD-profile diagnostics; **S9 Fig** shows every audited nuisance cell for the null and two alternative effects. S1–S6 Figs retain their earlier artwork.

S1 Appendix contains A1–A15 (core mechanisms and proofs), S2 Appendix D1–D7 (provenance and earlier observation analyses), and S3 Appendix E1–E7 (general algorithm, AceK evidence, expanded statistics and preserved extension). The old 100-trial benchmark remains explicitly marked as exploratory and superseded for the main statistical claim. Earlier figures and text remain in `docs/extended_analyses/`; historical filenames do not imply current figure numbering.

Original code is MIT licensed. The unchanged Gosztolai et al. workbook and model retain CC BY 4.0 attribution. Measurement-table provenance and cross-assay limits are recorded with each dataset. No primary literature PDFs are redistributed.
