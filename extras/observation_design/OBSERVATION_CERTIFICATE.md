# An observed PII count can make a cascade bound informative

## Quantity certified

The target is the maximum difference in total GS modified-site fraction between the specified free-PII-equivalent mechanisms. This is not a bound on error relative to the unknown biological truth, nor a guarantee of discrimination when the bound is large. All protein totals below count PII trimers, GS dodecamers, and GlnD/GlnE monomers.

Let epsilon=(D_T+2 E_T)/P_T<1 and e=E_T/S_T. Suppose a concurrently checked upper bound on the total PII mean UMP count is M. No independent PII-site assumption is required. If pi is the common free count distribution, nonnegative bound-pool modification gives m_free <= M/(1-epsilon). Since m_free >= 1-pi_0,

    pi_0 >= p_min = max(0, 1-M/(1-epsilon)).

For the source GlnE topology at fixed glutamine, the GS odds drive has form k(a+b/T), with a=(pi_0/K_P)(1+beta_prime*g/alpha_1), b=beta*g, and g=G/K_G. The common positive factor k contains the AR gain and recognition normalization. Defining ell=b/a, the free-PII pool T lies in [L,U]=[P_T(1-epsilon),P_T]. Therefore

    d = log[(1+ell/L)/(1+ell/U)],
    |Y_total,1-Y_total,2| <= min(1, e+tanh(nu*d/4)).

Here nu=12 covers every positive twelve-site GS ladder; nu=1 is sufficient for independent GS sites. The gain k cancels from d. Substitution of p_min into ell provides an upper bound because ell decreases with pi_0. If p_min=0, the generic log contrast -log(1-epsilon) is recovered. The mixture step and logit/variance inequality are the established ingredients; their composition turns a concurrently measured total PII count into a tighter cascade certificate.

The source rates and observation gate are not an unconditional native model fit. A gate must be checked in the assay in which the certificate is used. In particular, a rounded historical coordinate and an arbitrarily assigned error range are not a measured confidence interval.

## Complete prospective decision example

Fix P_T=5 uM, S_T=0.5 uM, D_T=0.5 uM and glutamine 2 mM. The protein scales are motivated by the reconstituted system; the proposed GlnE total is E_T=0.02 uM. Use the gate M<=0.67 UMP per trimer, motivated by the historical digitized point 0.65 at 2 mM but requiring a new concurrent measurement/calibration rather than transferring that point between assays.

The 2007 GlnE source-fit example gives K_G=15.6 mM, K_P=0.18 uM (main-text units), alpha_1=0.17 and catalytic rates 7.5,175,1315 for EP,EG,EGP. Thus beta=175/7.5 and beta_prime=1315/7.5. These are fitted coefficients, not independently known microscopic constants or a confidence set.

At E_T=0.02 uM, epsilon=0.108, e=0.04 and p_min=0.248879. The generic twelve-site bound is 0.370035; the observation-conditioned bound is 0.041176. The latter certifies the prescribed 0.05 tolerance. At the original E_T=0.10 uM scale, the bound remains above 0.20, so increasing PII alone cannot remove the GS sequestration allowance. Reducing GlnE fivefold is sufficient in this declared example without also reducing GlnD. It is not claimed to be the unique or globally optimal concentration choice.

The certificate can be stated directly as an experimentally testable coefficient requirement: ell<=0.141799 uM is sufficient at the proposed totals and tolerance. The nominal source-fit ell upper bound is 0.016240 uM, leaving an 8.73-fold coefficient margin. An alternative K_P=1 uM scale gives a bound of 0.046431; that apparent kinetic scale is a sensitivity scenario, not a measured binding constant. Literal use of the supplementary legend's conflicting 0.18 mM K_P gives 0.299534 and fails the certificate. This explicit unit sensitivity is retained.

Consequently, the decision is conditional but operational: verify the count gate and a coefficient upper bound, then choose the GlnE concentration. If the gate, units or coefficient bound cannot be justified, use the generic bound or obtain the missing binding measurement. A small computed error in one fitted realization is insufficient for this certificate.

## Checks

`observation_certificate.py` tests 5,000 arbitrary positive GS ladders over unknown catalytic gains and 5,000 independent free/bound PII mixtures. `verify_observation_certificate.py` tests 120 paired full 120-species conserved-reaction realizations with source-rate ratios, multiple gains, GlnE totals and K_P alternatives. No bound violations were found; maximum stationary vector-field residual was 1.23e-15. GlnE sigma is fixed at one in those source-constrained realizations, while GlnD lambda differs. These calculations verify the stated mechanism class; they do not establish the complete native cascade or equilibration times after lowering enzyme concentration.
