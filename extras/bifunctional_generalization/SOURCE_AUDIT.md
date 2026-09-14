# Independent AceK–IDH sources

These are transcriptions of published summary tables. They are not raw replicate
measurements, and reconstructed curves must not be labeled experimental points.

## Miller et al. (1996), JBC 271:19124–19128

DOI: https://doi.org/10.1074/jbc.271.32.19124

Table I provides initial-rate parameters for wild-type AceK, AceK3 and AceK4,
each in kinase and phosphatase assays. The concentration units are micromolar;
the text conversion on ResearchGate incorrectly renders the Greek mu as `m`.
Dexter & Gunawardena (2013), Table 2, independently confirms the WT protein
Michaelis constants as 2.3e-7 M and 4e-7 M.

The original rate law (Eq. 1) is

    v = Vmax A B / (Ka Kmb + Kma B + Kmb A + A B).

Table I also lists apparent dissociation constants and rounded interaction
ratios; these rounded quantities are not algebraically exact identities. Use
the printed Ka, Kma and Kmb in the rate law and retain Kb as an independently
reported rounded value. The statement that standard errors are below 20%
does not supply a covariance matrix or a simultaneous 95% confidence region.
Do not turn a ±20% sensitivity box into an experimental confidence set.

Accessible full article was read at:
https://www.researchgate.net/publication/14498934_Isocitrate_dehydrogenase_kinasephosphatase_Kinetic_characteristics_of_the_wild-type_and_two_mutant_proteins

The same paper reports ATPase maxima of 590, 710 and 650 nmol/mg/min and ATP halfpoints of 230, 770 and 2900 uM for WT, AceK3 and AceK4. The latter halfpoints are withheld from the kinase construction. ATPase maxima set an additional branch rate; they are not themselves held-out predictions. See `miller1996_atpase_holdout.csv` and the independently checked minimal embedding.

## LaPorte, Thorsness & Koshland (1985), JBC 260:10563–10568

DOI: https://doi.org/10.1016/S0021-9258(19)85122-0

Table II measures total IDH activity after dephosphorylation and activity during
growth on acetate, in units/mg. The table's fraction is rounded, so multiplying
the rounded fraction by total activity need not reproduce the independently
printed activity. DEK2010 without a plasmid and with pBR322 have duplicate
reported values; do not count them as independent experimental replicates.
The pCK505 overexpression pair compares total activity 2.3 with 34.5 units/mg,
while the reported active activity is 0.6 units/mg in both conditions.

Units/mg cannot be silently converted into micromolar IDH dimers. Conversion
requires specific activity, intracellular volume and normalization metadata.
These in vivo observations are also not matched to Miller's in vitro assays.

Accessible full article was read at:
https://www.researchgate.net/publication/19265162_Compensatory_phosphorylation_of_isocitrate_dehydrogenase_A_mechanism_for_adaptation_to_the_intracellular_environment

## Dexter & Gunawardena (2013), JBC 288:5770–5778

DOI: https://doi.org/10.1074/jbc.M112.339226
Open primary article: https://pmc.ncbi.nlm.nih.gov/articles/PMC3581427/

Table 1 specifies the ten-species conserved mass-action network (free AceK
plus nine IDH species). One IDH dimer can bind two AceK molecules. Table 2
assigns on-rates 1e8 M^-1 s^-1 (assumed, not measured), off-rates 23 and 40 s^-1,
and catalytic rates 0.098 and 0.065 s^-1. The latter derive from Miller's WT
measurements. In micromolar units the on-rates are 100 uM^-1 s^-1.

The paper's nearly constant-Ioo illustration additionally fixes Iop. That is
an observation/biological assumption, not a consequence of imposing just the
two mass conservation laws. A new analysis must retain this distinction.
The biochemical data support the model's topology and coefficients at
different levels; neither the old paper nor the new computations establish
an unrestricted native mechanism.
