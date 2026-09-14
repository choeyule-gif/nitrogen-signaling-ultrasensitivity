# Verification of revision 4 and v1.3.1 corrections

The principal numerical conclusions reproduce. The verification patch changes two printed sufficient thresholds and citation/package metadata; it does not change the numerical model, fitted results or figure artwork.

## Public-package regeneration

The public v1.3.0 release was downloaded and all three asset hashes checked. The outer upload manifest (16 files) and inner data manifest (451 files) matched. Nine generated-output directories were deleted from an extracted source copy before installation into a fresh virtual environment. The 24 scheduled regeneration and verification stages all completed successfully in 809.76 seconds. The baseline had 22 passing checks and no skips. The nine central observation-design files, including the 8,700 trial records, were byte-identical to the public copies.

This is a regeneration of source commit `39774c720ef713fc20a007748e74ace131cc80af` (v1.3.0). The current v1.3.1 suite adds one independent audit stage; it has 25 stages with regeneration and eight in default check mode. The additional check passed both on the regenerated source/results and on the final working tree. Optional exploratory scripts outside the scheduled suite were not all rerun. Model calculations are unchanged by the patch.

## Mathematical and numerical checks

- Independent 65-digit arithmetic gives the nominal total-GS bound **0.04117568260412924**, consistent with the reported 4.12%.
- The exact coefficient cutoff for a 0.05 tolerance is **0.1417991000508063 micromolar**. The former main-text cutoff 0.142 was rounded upward and gives a bound **0.05001375278416886**, slightly above the guarantee. The printed sufficient threshold is now **0.1417 micromolar**, explicitly rounded down; it gives **0.04999321559770619**.
- The generic twelve-site pool threshold for delta=0.01 and epsilon_S=0.001 is **0.002995585257912924**. The former 0.002996 was also rounded upward. It is now **0.002995**. The independent-site value 0.03536 was already conservative and is unchanged.
- Fifteen independent linear programs recover the free zero-count lower bound. Nine positive-ladder endpoint constructions recover the tanh constant to better than 1e-45. These establish sharpness of the component inequalities, not of every full-cascade realization.
- Two source-rate mass-action BDF integrations start from free, unmodified target pools and unliganded enzymes. At G=2 mM, P_T=5, S_T=0.5, D_T=0.5 and E_T=0.02 micromolar, gain=0.001 and lambda=1 and 4, both approach the constructed stationary solution. Maximum species error is below 1.3e-11 micromolar, protein-balance error below 3.3e-13, and final vector-field residual below 5e-16. The time scale is arbitrary; two trajectories do not establish global convergence or physiological equilibration time.
- A separate free-ligand root calculation and Levenberg-Marquardt fit recover the archived raw binding example, differing in halfpoints by at most 8.5e-9 micromolar. Independent tail-probability calculations reproduce all 8,700 trial decisions.

## Primary sources and document integrity

All 36 bibliography DOIs resolve to matching titles in Crossref metadata. This check does not validate every scientific inference attributed to those sources. The directly reviewed primary pages confirm the 1998 Table 2 catalytic values, the 1998/2011 MgCl2 assay differences, and the 2007 GlnE source-fit coefficients.

The 2007 main article reports K_P=0.18 micromolar, while its S10 supplementary legend prints 0.18 millimolar. This inconsistency is in the original source. The manuscript already discloses it and retains the literal-unit sensitivity calculation (bound approximately 0.300). Source coefficients are fitted examples, not independently validated microscopic parameters for the proposed assay.

All ten documents in the public portable LaTeX archive compiled twice without missing files, undefined references or overfull boxes. The final main manuscript remains 18 pages, S1 Appendix 20 pages and S2 Appendix 13 pages. Figure artwork is preserved. Updated manuscript pages were rendered and inspected.

## Metadata corrections

The CFF schema identifier was incorrectly changed to 1.3.0 during versioning. It is restored to `cff-version: 1.2.0`, distinct from software `version: 1.3.1`, and validates against the [official CFF schema](https://raw.githubusercontent.com/citation-file-format/citation-file-format/main/schema.json). Package metadata previously reported 1.2.0 despite the v1.3.0 release; the current package, citation and archive metadata consistently report software version 1.3.1. Dependency bounds are no longer described as proof that every allowed version was tested; the verified environments are recorded explicitly.

## Remaining scientific scope

The computed bound is a pairwise guarantee within the declared reaction topology, shared free input, protein conservation and shared positive ladder. The PII count gate and coefficient envelope still require justification in the actual assay. Synthetic noise scenarios and nuisance ranges are not biological confidence regions. Hidden binding heterogeneity can defeat mechanistic attribution despite a good single-class curve fit. No new biological experiments, native-mechanism identification, global stability theorem or journal acceptance guarantee follows from this audit.

Machine-readable evidence is in the accompanying validation files and `extras/observation_design/results/independent_verification.json`. The previous released archive is preserved; the author manages Zenodo and no unverified new DOI is asserted.
