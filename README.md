# Mechanistic equivalence and measurement design in models of bacterial nitrogen signaling

Code, numerical outputs and manuscript sources for Yule Choi's nitrogen-signaling study. This is a conditional analysis of mechanistic equivalence and measurement design; the calculations do not establish a unique native GlnD–PII–GlnE–GS mechanism.

**Archived release:** [v1.1.0](https://doi.org/10.5281/zenodo.22735182). **Current release [v1.2.0](https://github.com/choeyule-gif/nitrogen-signaling-ultrasensitivity/releases/tag/v1.2.0):** adds independent paired-paralog data, assay interval calculations and integrated validation. These additions are not contained in the v1.1.0 archive. The v1.2.0 Zenodo version DOI is pending verification. The older release ZIP is not the current manuscript.

## Reproduce

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
python tests/run_research_suite.py --regenerate
```

The command regenerates the baseline, matched-model, finite-mechanism, regulatory-state, robust-design and closed-cascade numerical layers and runs their checks, followed by the paired-reporter checks. Run without `--regenerate` to check saved baseline results while freshly solving extension identities. See `REPRODUCIBILITY.md` for scope and independent-data extraction. A passing calculation verifies its stated model, not physiological truth or global parameter optimality.

Build the full current S1 Data archive after regeneration with `python scripts/package_research_data.py`. The older `scripts/make_s1_data.py` called by the baseline runner produces a baseline export folder only; it is not the submission archive.

## Research layout

| Location | Content |
|---|---|
| `esbm/`, `scripts/`, `data/`, `results/`, `figdata/` | Baseline theory, digitized titration, calculations and plot inputs |
| `extras/matched_models/` | Matched titrations, different ladders and downstream calibration |
| `extras/regulatory_states/` | Explicit GlnD/GlnE regulatory-state equivalences |
| `extras/finite_mechanism/` | Finite reaction kinetics and assumption boundaries |
| `extras/closed_cascade/` | Protein-conserving cascades, topology and free/total input bounds |
| `extras/robust_design/` | Nuisance-aware occupancy and binding designs |
| `extras/paired_paralogs/` | Independent published GlnB/GlnK data and calibrated reporter design |
| `extras/revision_audited/`, `matlab/` | Figure calculations and rendering inputs; some filenames retain older numbering |
| `manuscript/` | Working manuscript, appendix and supporting figures |
| `validation/`, `tests/` | Clean-environment evidence and executable checks |

## Current figure mapping

The compressed working manuscript has five main figures and 17 supporting figures. The original pathway diagram is unchanged.

| Main figure | Subject | Previous main figure |
|---|---|---|
| 1 | Author's pathway diagram | 1 |
| 2 | GlnD regulatory equivalence and complementary observables | 7 |
| 3 | Protein-conserving cascade and topology bounds | 10 |
| 4 | Nuisance-aware measurement design | 11 |
| 5 | Independent paired-paralog data and conditional assay | New |

Previous main figures 2, 3, 4, 5, 6, 8 and 9 are now S11, S12, S13, S14, S15, S16 and S17 Figs. Original S1–S10 Figs remain. All extended pre-compression result and method arguments are retained in S2 Appendix; S1 Appendix contains S1–S29, including prior-art comparison and the paired-reporter derivation. Old script names are not current submission figure numbers.

Code is MIT licensed. The Gosztolai et al. source workbook and model retain their CC BY 4.0 license and attribution, documented in `extras/paired_paralogs/README.md`. Other published measurements retain their source provenance in `data/README.md`.
