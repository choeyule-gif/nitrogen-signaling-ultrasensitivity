# data/

Inputs only. Nothing in this directory is produced by the code.

| file | content | source |
|---|---|---|
| `jiangninfa2011_fig2B_digitised.csv` | twelve points (glutamine, mM; mean UMP per PII trimer, 0–3) read from Fig. 2B of Jiang & Ninfa, *Biochemistry* 50:10929 (2011) — the UTase/UR–PII part of the reconstituted bicyclic system, PII 0.5 µM, UTase/UR 0.05 µM; the authors report n_H ≈ 1.95 and S₀.₅ = 0.55 mM for this curve | digitised by eye from the enlarged printed figure |
| `published_constants.json` | Hill coefficients, half-points, concentrations and affinities quoted from the cited papers | Ventura et al. 2010; Jiang & Ninfa 2011; Jiang et al. 1998a,b |

**Digitisation uncertainty model** (used by every bootstrap in `scripts/`): each ordinate is
perturbed by ±0.06 (12 % of the 0.5 tick spacing) and each abscissa by ±6 %, independently,
uniform in the bootstrap of `cascade_budget.py` and Gaussian in the bootstrap of
`fit_titration.py`; both are stated in the manuscript's Methods.
