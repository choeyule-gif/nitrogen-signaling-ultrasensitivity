# A graph-to-observation workflow for bifunctional regulators

This is a working mathematical extension. Empirical validation and comparison
against existing prediction-design methods remain to be completed before the
extension can be presented as a finished result in the submitted manuscript.

## 1. Define the class before claiming enzyme cancellation

Fix the free external ligands and free converter-enzyme concentration e. Assume
each complex contains exactly one target molecule (which may itself be an
oligomer), and all target-state transitions are first order at fixed e. Let
Q(e) be the irreducible column-conservative transition generator on these states.
Binding a second converter to a target is allowed. Aggregating several target
molecules into one complex, changing free ligands, synthesis/degradation, and
multiple stationary classes require an enlarged model or a separate argument.

By the established matrix-tree theorem, a stationary distribution is

    p_i(e) = tau_i(e) / sum_j tau_j(e),

where tau_i is the sum of positive directed spanning-tree weights rooted at i.
If edge rates are polynomials in e with nonnegative coefficients, tau_i also
has nonnegative coefficients. This result and graph elimination are prior
theory; novelty is not claimed for them.

For a designated subset F of unbound target states, free-state enzyme
cancellation holds if and only if the polynomial vector tau_F(e) has rank one
as a coefficient matrix. Equivalently, every ratio tau_i(e)/tau_j(e), i,j in F,
is constant. A common factor in the tree polynomials does not change this
criterion: all polynomials must be scalar multiples of the same polynomial.
Thus bifunctionality alone does not establish enzyme-independent free states.

Proof: equality of normalized free distributions implies tau_i(e)=c_i H(e)
on a positive open interval. Polynomial identity extends the equality to every
e. The converse follows immediately on normalization.

For a fixed e, positive column scaling Q_lambda=Q D_lambda changes stationary
weights to D_lambda^-1 tau. Scaling only bound-state columns leaves normalized
free states identical, while changing dwell times and bound fractions. This
is a constructive family, not proof that measured initial-rate assays permit
arbitrary scaling. Such assays must constrain D_lambda when included in the
observation set.

## 2. Enclose positive stationary polynomials on an enzyme interval

Let tau_F(e)=C(1,e,...,e^n)^T and restrict e to [L,U]. Express each polynomial
in the degree-n Bernstein basis on that interval:

    tau_i(L+(U-L)t) = sum_k b_ik binom(n,k) t^k(1-t)^(n-k).

For positive power-basis coefficients and L>=0, every b_ik is nonnegative.
Normalize each nonzero column to v_k=b_:k/sum_i b_ik. Then

    pi_free(e) belongs to conv{v_0,...,v_n}.

Proof: normalize the preceding polynomial sum; the remaining column weights
are nonnegative and sum to one. This is an outer enclosure, not a claim that
every convex combination is kinetically realizable. Subdivision of [L,U]
and union of the resulting polytopes can tighten it. Rank-one cancellation
produces a one-point enclosure. Noncancelling networks are handled by the same
algorithm rather than being forced into a cancellation approximation.

For an observable h in [0,1], its free-state range is bounded by
max_k h.v_k - min_k h.v_k. A distribution-wide bound is max_kl TV(v_k,v_l).
This enclosure uses all supplied kinetic coefficients; it can be much tighter
than a bound using polynomial degree alone. It complements an exponential-family
log-oscillation bound rather than replacing that established result.

## 3. Condition on measurements with conservation-aware linear programs

Let V contain the enclosure vertices. Let lambda>=0 weight the free pool and
z>=0 represent the bound pool, categorized by the same modification count.
If the target-bound fraction is at most epsilon, impose

    sum(lambda)+sum(z)=1, sum(z)<=epsilon.

The total count distribution is V lambda+z. For matched observations with
matrix A and bounds [a,b], add a<=A(V lambda+z)<=b. Minimize and maximize the
desired linear reporter h.(V lambda+z) over this polytope. Every compatible
stationary mechanism lies inside the resulting interval. Its width bounds
pairwise disagreement within that class; it is not error relative to native
biology. If the polytope is empty, the observation interval and model class
are incompatible. An empty set is not a successful certificate.

The observation intervals must come from measured uncertainty or be labeled
prospective thresholds. Source-table standard errors, digitization errors,
and arbitrary parameter sensitivity boxes must not be exchanged silently.

## 4. Independent topology test: AceK–IDH

Use the nine target states in Dexter & Gunawardena (2013), Table 1, plus free
AceK through conservation. Both singly and doubly AceK-bound IDH dimers are
included. Eliminating the three doubly bound states gives, with B_j denoting
singly bound dimers,

    B_0/B_1 = (k6 + e*b)/(k3 + e*a),
    B_2/B_1 = (k7 + e*c)/(k10 + e*d),

where a=k17*k20/(k18+k20), b=k11*k14/(k12+k13+k14),
c=k11*k13/(k12+k13+k14), d=k15*k19/(k16+k19).

The unbound states satisfy

    e*k1*U0 = k2*B0+k6*B1,
    e*k4*U1 = k3*B0+k5*B1+k10*B2,
    e*k8*U2 = k7*B1+k9*B2.

Clearing the two positive linear denominators yields three positive
quadratic polynomials. Their coefficient rank tests cancellation, and their
Bernstein enclosure supplies the observation-conditioned bound. This extends
the workflow beyond the original single-converter-binding ladder while
retaining the actual second-enzyme complexes in the published model.

Published kinetic values and the independently reported in vivo activity
table are transcribed separately with their assay and unit limitations.
The numerical example must not treat an in vivo units/mg value as an in vitro
micromolar concentration, or a fixed monophosphorylated IDH level as a derived
conservation law.

## Comparison that still requires quantitative completion

Prediction profile likelihood supplies likelihood-based prediction intervals;
model-discrimination OED optimizes prospective differences/information. This
workflow instead constructs a conservative class enclosure before selecting
a measurement. It is not automatically superior: the enclosure may be wider,
and profile methods may use likelihood information that this polytope omits.
A fair comparison must use the same kinetic class, observation model, data,
noise assumptions, and design budget, and report both informativeness and
computational cost. A source-comparison table alone is not that comparison.
