# Calibrated response curves and sharp sequestration intervals

This extension replaces the binomial odds identity with a calibrated free-pool response curve. It permits bounded reporter sequestration and partially modified bound target. It supplies an assay framework, not an estimate of native occupancy from the published time series.

## 1. Calibration beyond independent sites

Let a free target have a positive ladder distribution p_i(u) proportional to a_i u^i, i=0,...,n. Its mean modified-site fraction is m(u)=sum(i p_i)/n. Differentiation gives d m / d log u = Var(i)/n > 0 for finite positive u. The variance identity is established binding-polynomial mathematics, not a new theorem.

For two ladders with the same scalar drive, define C=m_K composed with inverse(m_B). Then q_K,free=C(q_B,free). This relation need not be a ratio of odds: positive interacting ladders are allowed. A no-binding titration can measure C directly, without separately identifying the two ladders or catalytic drive. Calibration must cover the assay range; extrapolation is not guaranteed. A full curve requires more controls than the one-point calibration sufficient for binomial ladders with a fixed relative specificity.

The observation construction also applies to an empirical, independently validated curve outside this ladder class. In that case the curve's transferability, monotonicity and error envelope require experimental support. Multiple independent inputs, condition-dependent recognition or hysteresis need not collapse to one curve. Repeated controls across ATP, 2-OG and protein concentration test that collapse; they do not follow from the algebra. A shared dynamic protocol can be calibrated empirically, but a stationary calibration cannot be transferred to a transient assay without a justified discrepancy bound.

## 2. A reporter need not be completely unsequestered

Let its measured total fraction lie in [B_L,B_U]. Suppose an independent experiment bounds the fraction of reporter protein sequestered by epsilon<1, with unknown modification in the bound pool. Mixture conservation gives the sharp free-reporter interval

    B_free in [max(0,(B_L-epsilon)/(1-epsilon)), min(1,B_U/(1-epsilon))].

For epsilon=1 the interval is [0,1]. Sharpness follows by selecting bound-pool modification at its extreme and an allowed bound fraction. Thus nonzero reporter sequestration is allowed, but the method cannot create its bound from an uncalibrated reporter measurement alone.

For an increasing calibrated curve C with a jointly valid absolute transfer error delta, the free K interval is

    F_L=max(0,C(B_free,L)-delta),
    F_U=min(1,C(B_free,U)+delta).

Calibration uncertainty can be included through lower and upper curve envelopes evaluated over the full allowed reporter interval. These envelopes must cover both calibration error and biological transfer error. Without monotonicity, use the minimum and maximum over the reporter interval rather than just its endpoints.

## 3. Imperfect modification selectivity

Let free K modification f lie in [F_L,F_U], bound K modification b in [H_L,H_U], and total K modification q in [Q_L,Q_U]. For a bound protein fraction s,

    q=(1-s)f+s*b.

The sharp identified set is the intersection of [0,1] with

    (H_L-F_L)*s <= Q_U-F_L,
    (F_U-H_U)*s <= F_U-Q_L.

Proof: at fixed s the continuous attainable mixture is exactly [(1-s)F_L+s H_L, (1-s)F_U+s H_U]. This interval intersects the observed interval exactly when the two displayed inequalities hold. Every point in the resulting set has a feasible mixture; hence the bound is sharp for the stated rectangular uncertainty set. Dependencies between uncertainties can tighten it; the rectangular set is then a conservative relaxation, not necessarily sharp for the richer model.

No site independence, homotrimer composition or stationarity is required for this mixture identity. Those conditions can still be needed to justify the supplied free-response calibration. Here s is a protein-subunit fraction. Converting it to an oligomer count requires composition information.

If the free and bound responses are identical, the measurement cannot distinguish s. If the two inequalities have no feasible s, the joint calibration, selectivity and observation assumptions are incompatible. That inconsistency does not identify which assumption failed.

## 4. What becomes experimentally distinguishable

With f in [0.70,0.80] and q in [0.35,0.39], assuming b=0 yields s in [0.443,0.563]. Permitting b up to 0.15 widens this to [0.443,0.692]. With no bound-pool selectivity information, b in [0,1], the upper endpoint becomes 1 but the positive lower endpoint persists. These are synthetic numbers, not estimates from Gosztolai et al.

This separates three experimental tasks: a no-binding control constrains the free response; reporter partition measurements bound contamination of the calibration coordinate; direct analysis of bound K constrains selectivity. A direct bound-fraction measurement is then a validation endpoint, rather than an input substituted for the desired result. Controls should be selected before evaluating a held-out binding condition.

For fixed f>b and exact measurements, s=(f-q)/(f-b). A small separation f-b amplifies measurement and calibration error, so an assay should seek independently verified free/bound contrast, not simply a steep total titration. This is a consequence of mixture inversion, not a claim that a particular physiological input already provides that contrast.

## 5. Independent computational checks and limits

`code/robust_calibration.py` checks 3,000 interval configurations against independent linear programs in bound fraction and free/bound modified masses. It also checks curve reconstruction in 1,000 positive ladder pairs of different sizes and interacting coefficients. `results/robust_calibration_verification.json` records discrepancies and synthetic examples. The integrated research suite includes this extension.

These tests validate the calculations. They do not validate native stationarity, empirical calibration transfer, assembly composition, bound-pool selectivity, or actual AmtB occupancy. The published peptide time series lacks the matched no-binding calibration and direct occupancy endpoint needed for such a validation. The contribution is the explicit distinction between experimentally supplied bounds and consequences that can be inferred sharply from them; the underlying mixture mathematics is elementary.
