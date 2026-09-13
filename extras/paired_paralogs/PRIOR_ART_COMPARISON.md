# Boundary between established results and current contributions

This is a research comparison, not a claim that a literature search proves priority.

| Result | Established foundation | Present contribution / required revision |
|---|---|---|
| Adjacent stationary ratios and adjacent-state contrast | Gunawardena (2007), Eqs. 3–5, DOI 10.1529/biophysj.107.110866 | Cite directly. The contrast in current Eq. 14 is an application of the invariant, not a new theorem. Model-dependent constancy can test the ladder, but does not uniquely establish a native mechanism. |
| Elimination of catalytic complexes and rational stationary forms | Thomson & Gunawardena (2009), DOI 10.1016/j.jtbi.2009.09.003; Feliu & Wiuf (2013), DOI 10.1007/s00285-012-0510-4 | Treat as the mathematical method. State explicitly which reaction restrictions permit reduction. |
| Scaling non-identifiability | Castro & de Boer (2020), DOI 10.1371/journal.pcbi.1008248; Villaverde & Massonis (2021), DOI 10.1371/journal.pcbi.1009032 | The GlnD/GlnE transformations are sufficient constructive examples, not an exhaustive identifiability test or invention of scaling symmetry. |
| Free/total observation bounds in a closed cascade | Built on conservation and finite-support response bounds | Candidate central result: combine upstream sequestration capacity, downstream drive elasticity and GS ladder support into a finite error envelope; distinguish topology-dependent identities from the bound that survives a productive EPU route. Further prior-art search is needed before a general novelty claim. |
| Complementary binding or kinetic designs | Experimental design can determine identifiability; Silk et al. (2014), DOI 10.1371/journal.pcbi.1003650 | Specific assay predictions for stated biochemical families must survive nuisance variation and independent data. A generic recommendation to measure another observable is not the new contribution. |
| Oligomer-state bounds from peptide-level means | Finite-support moment inequalities / linear programming | New biological application under investigation: quantify what existing GlnB/GlnK peptide measurements can imply about intact-trimer availability without assuming independent sites. The inequality itself is elementary and must not be advertised as a new general theorem. |

## Independent-data study being tested

Gosztolai et al. (2017), DOI 10.1016/j.bpj.2017.04.012, data/code DOI 10.6084/m9.figshare.4880003. The released workbook supplies means, SE and n, not replicate vectors. Measurements concern modified peptides and total abundance. Two observation/calibration models (modified/total and modified/(modified+unmodified)) must not be silently interchanged.

The original glutamine factors 1/(1+g) and g/(1+g) admit a two-ligand realization with weights (1,2g,g²) and catalytic coefficients (1,1/2,0) and (0,1/2,1). Scaling the singly liganded weight and its coefficients oppositely preserves the stationary rate ratio but changes their common speed. This gives a direct way to test stationary-equivalent families against the published time-course setting. The initial numerical port is a pilot: interpolation, changing protein pools, nonnegative accessible pools and input uncertainty must pass before any biological inference.
