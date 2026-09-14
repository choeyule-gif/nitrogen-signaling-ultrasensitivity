# Response to pre-submission review: revision 4, v1.3.2

## A. GlnD fit comparison

The review identified a valid error in the interpretation of the RMSE comparison. The earlier fit was described as a constant-capacity model, and the improvement was attributed to varying state capacities. In fact, its capacities are proportional to (1,1,d2), not all equal, and the newer fit releases several restrictions simultaneously.

The main Results and S2 Appendix D4 now identify these differences explicitly:

| Constraint | Earlier RMSE 0.539934 model | Newer RMSE 0.040935 model |
|---|---|---|
| Unliganded UT/UR specificity ratio | 38.9012346 | Retained |
| UR endpoint specificity fold | 6.7524842 | Retained |
| Effective coefficient ratio C | Fixed at 100 | Not fixed |
| Intermediate-state constraint | q1=q0 | Released; monotone state activities retained |
| UR activity half-range | Imposed at 0.080 mM | Not imposed; fitted value 0.256624292315 mM |

The UR half-range means the glutamine concentration at which the free-enzyme UR capacity R(G)=sum(b_j w_j)/sum(w_j) reaches (b0+b2)/2. It is distinct from the titration midpoint, the singly liganded occupancy maximum, and a microscopic binding constant.

The independent verification reconstructs the older model from its saved parameters and obtains half-ranges 0.080000000000 and 0.256624292315 mM using both a scalar root and the equivalent quadratic. It checks both endpoint specificity ratios and the nonuniform older capacities. The newer GlnD fit was rerun: all fitted parameters, predictions, RMSE and equivalence-family values are unchanged; scope metadata is corrected.

The revised explanation attributes the comparison to joint relaxation of fixed-C, intermediate-state and transported half-range constraints. It makes no single-relaxation attribution and no claim that all previous kinetic constraints remain satisfied. The boundary affinity and profile retain their conditional interpretation. Fig4D's caption now specifies the two retained endpoint ratios and the released UR half-range.

## B. GlnE range and Monte Carlo interpretation

The main Results and S2 Appendix D5 now distinguish the earlier analytic K_G range [7.8,31.2] mM, a factor-two sensitivity range around 15.6 mM, from the newer source-anchored [7.8,15.6] mM observation benchmark. The latter spans the source's direct-scale comparison and joint-fit value. Neither range is a measured confidence region, and the newer Monte Carlo results are not claimed to cover the wider range.

The main text explicitly states that 100 trials do not establish uniform 5% false-positive control. S2 retains the Wilson intervals, including approximately [0.028,0.125] for 6/100 rejections, and limits conclusions to the tested configurations. No new power or error-control result is claimed.

## C. Pairwise certificate

The 0.04117568 bound remains a guarantee on differences between members of the declared mechanism class. The qualification distinguishing this from error relative to biological truth is retained. Figure 3 artwork and its numerical results are unchanged.

## D. Fig5C reference

The caption now points directly to "the control and assay ranges in S1 Appendix A13." Figure 5 remains a prospective application and no new biological validation is claimed.

## E. Verification scope

This patch reran the GlnD fit and the default eight-stage integrated checks. These include fresh reaction integrations and binding refits, and recomputation of all 8,700 saved trial decisions. It did not regenerate those 8,700 trials from raw observations again. The full 24-stage clean regeneration remains attributed to the previously audited public v1.3.0 source. The suite still contains 25 stages when regeneration is requested.

No external biochemical source was newly audited in this patch. The earlier direct primary-source checks and DOI/title metadata checks retain their stated scope. Figures are preserved. Source, manuscripts, metadata and submission files are synchronized as v1.3.2; earlier releases and verified archives are preserved.
