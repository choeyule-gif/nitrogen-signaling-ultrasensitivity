# Mechanistic equivalence and measurement design in models of bacterial nitrogen signaling

Code, numerical outputs and manuscript sources for Yule Choi's nitrogen-signaling study. This is a conditional analysis of mechanistic equivalence and measurement design; the calculations do not establish a unique native GlnD–PII–GlnE–GS mechanism.

**Current revision: v1.3.3.** This final submission-preparation patch sharpens and compresses the initial-submission cover letter. It retains the v1.3.2 GlnD comparison, assay-range and Fig5C corrections and the v1.3.1 verification fixes. The v1.3.0 extension adds the observation-conditioned concentration certificate, end-to-end binding measurements and biochemical constraint profiles, and synchronizes the main manuscript, appendices and Figures 3–4. The earlier [v1.1.0 archive](https://doi.org/10.5281/zenodo.22735182) and v1.2.0 release do not contain these extensions. A version DOI for this revision is not asserted until its archival record is verified.

## Reproduce

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
python tests/run_research_suite.py --regenerate
```

The command regenerates the baseline, matched-model, finite-mechanism, regulatory-state, robust-design and closed-cascade numerical layers and runs their checks, followed by the paired-reporter and observation-design layers. Run without `--regenerate` to check saved baseline results while freshly solving extension identities. See `REPRODUCIBILITY.md` for scope and independent-data extraction. A passing calculation verifies its stated model, not physiological truth or global parameter optimality.

Build the full current S1 Data archive after regeneration with `python scripts/package_research_data.py`. The older `scripts/make_s1_data.py` called by the baseline runner produces a baseline export folder only; it is not the submission archive.

## Research layout

| Location | Content |
|---|---|
| `esbm/`, `scripts/`, `data/`, `results/`, `figdata/` | Baseline theory, digitized titration, calculations and plot inputs |
| `extras/matched_models/` | Matched titrations, different ladders and downstream calibration |
| `extras/regulatory_states/` | Explicit GlnD/GlnE regulatory-state equivalences |
| `extras/finite_mechanism/` | Finite reaction kinetics and assumption boundaries |
| `extras/closed_cascade/` | Protein-conserving cascades, topology and free/total input bounds |
| `extras/observation_design/` | Count-conditioned bounds, raw binding simulation, concentration controls and constraint profiles |
| `extras/robust_design/` | Nuisance-aware occupancy and binding designs |
| `extras/paired_paralogs/` | Independent published GlnB/GlnK data and calibrated reporter design |
| `extras/revision_audited/`, `matlab/` | Figure calculations and rendering inputs; some filenames retain older numbering |
| `manuscript/` | Working manuscript, appendix and supporting figures |
| `validation/`, `tests/` | Clean-environment evidence and executable checks |

## Current figure mapping

The compressed working manuscript has five main figures and six supporting figures. The original pathway diagram is unchanged.

| Main figure | Subject | Previous main figure |
|---|---|---|
| 1 | Author's pathway diagram | 1 |
| 2 | GlnD regulatory equivalence and complementary observables | 7 |
| 3 | Count-conditioned concentration certificate and topology limits | Redesigned from 10 |
| 4 | Raw binding data, discrimination, calibration and GlnD constraint profile | Replaces analytic-budget main panels |
| 5 | Independent paired-paralog data and conditional assay | New |

Current S1–S6 Figs correspond to the v1.2.0 S11, S12, S13, S14, S16 and S17 figures. S1 Appendix (A1–A15) contains core mechanisms, proofs and prior art; S2 Appendix (D1–D7) contains data and sensitivity analyses. The unchanged former supplements and all 17 older figures are preserved in `docs/extended_analyses/v1.2.0/`. Old script/input filenames are not current submission figure numbers. The generalized calibration extension and observation-design layer are included in v1.3.0 and S1 Data. Prior reviewer-revision text is also retained under `docs/extended_analyses/reviewer_revision3/`; current claims and numbering are in `manuscript/`.

Code is MIT licensed. The Gosztolai et al. source workbook and model retain their CC BY 4.0 license and attribution, documented in `extras/paired_paralogs/README.md`. Other published measurements retain their source provenance in `data/README.md`.
