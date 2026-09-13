from pathlib import Path
import shutil,re
BASE=Path(__file__).resolve().parents[3];ROOT=Path(__file__).resolve().parents[1]
OLD=BASE/'output/pdf/plos_regulatory_states';NEW=BASE/'output/pdf/plos_robust_design';NEW.mkdir(exist_ok=True)
for pat in ['*.sty','*.bst']:
 for p in OLD.glob(pat):shutil.copy2(p,NEW/p.name)
shutil.copytree(OLD/'figures',NEW/'figures',dirs_exist_ok=True)
for n,name in [(9,'Fig9_assumption_boundary'),(10,'Fig10_discriminating_design')]:
 for ext in ['pdf','png']:shutil.copy2(ROOT/f'figures/{name}.{ext}',NEW/f'figures/Fig{n}.{ext}')
RESULT=r'''
\subsection*{Finite reaction kinetics delimit equivalence in free and total observations}

The GlnD equivalence does not require instantaneous ligand equilibration at steady state. We replace averaged catalytic capacities by three free enzyme states, six enzyme--target complexes, and two target modification states, with finite mass-action association, dissociation, and catalysis. When regulatory ligand exchange occurs on free enzyme, complex balances recover the same catalytic specificities and free-target ratio at any positive binding-rate scale (S22). Explicit saturation and enzyme conservation therefore retain the free-state identity. The target here is one modification unit; this calculation does not assign the same sequestration stoichiometry to a whole PII trimer.

Total modification can nevertheless differ because complexes are counted in the assay (Fig 9A,B). For one target unit per enzyme complex, two mechanisms with the same free modified fraction and the same totals obey
\[
 |Y_{\mathrm{tot}}^{(1)}-Y_{\mathrm{tot}}^{(2)}|
 \leq \min(1,E_{\mathrm{tot}}/T_{\mathrm{tot}}).
\]
The bound follows from conservation, not from a low-occupancy expansion. If forward and reverse effective substrate affinities are equal within each enzyme state, total and free fractions coincide even at saturation. Otherwise the enzyme/target ratio bounds how much the observation rule alone can separate the models. This is a fraction bound, not a bound on titration slopes or Hill coefficients. It also limits experimental information: if $E_T/T_T=0.01$ and a single total-fraction estimate has known Gaussian SD 0.05, even a separation attaining the bound requires at least 155 independent estimates for 80\% power at one-sided size 0.05. The calculation is an optimistic lower bound under specified errors, not a measured assay budget; it motivates a complementary observable when sequestration effects are small.

Finite ligand-binding rates change transients even when the endpoints remain identical (Fig 9C). Extending ligand exchange to substrate-bound enzyme probes a different restriction. Exchange that balances the catalytic-complex weights retains the constructed steady state. However, thermodynamically consistent binding cycles with state-dependent dissociation/catalysis ratios can produce regulatory currents and break the identity. In a disclosed scenario, the free-fraction difference reaches 0.035 at $G=10K$ and 0.070 at $G=30K$ across the scanned exchange rates (Fig 9D; S22). These are sensitivity scenarios, not estimated physiological effects. They show why persistence under slow binding and failure under altered catalytic coupling are different questions.

'''
DESIGN=r'''
\subsection*{Complementary binding measurements remain informative after nuisance variation}

A useful experiment must distinguish mechanisms after allowing unknown quantities to vary within each candidate, rather than only separate two fixed prediction curves \cite{r70}. A single-condition relaxation curve can be made identical by refitting an unknown common speed. Likewise, a single GlnE binding isotherm can be matched by rescaling an unknown PII affinity. We therefore seek observations that eliminate these nuisance scales or retain separated prediction ranges (S23).

For GlnD, measure the singly liganded enzyme fraction without target at the nominal binding concentration. Allow the true $K$ to range independently from half to twice that concentration in each mechanism. The $\lambda=1$ model then predicts a fraction between 0.444 and 0.500, while $\lambda=4$ predicts 0.167--0.200 (Fig 10A). Their minimum separation is 0.244, independent of a shared catalytic speed. In contrast, mean ligand occupancy at the central symmetry point is uninformative. The target-free binding assay also avoids the substrate-sequestration and bound-state catalytic currents examined in Fig 9.

For GlnE, perform PII-binding titrations with PII-UMP absent, at zero and nonzero buffered glutamine. The normalized halfpoint is
\[
 r(G)=\frac{K_{P,\mathrm{eff}}(G)}{K_{P,\mathrm{eff}}(0)}
      =\frac{1+G/K_G}{1+G/(\alpha_1K_G)}.
\]
The unknown PII affinity cancels, and removing PII-UMP removes its interaction parameter. If enzyme binding appreciably depletes PII, the nominal total-PII halfpoint is $K_{P,\mathrm{eff}}+E_{\mathrm{tot}}/2$; subtracting the measured half-enzyme concentration restores the free halfpoint. Titrations at 0, 1, and 10 mM glutamine can additionally identify both $K_G$ and $\alpha_1$ in the ideal six-state binding model, except at the no-synergy degeneracy $\alpha_1=1$ (S23). This supplies a continuous parameter-estimation design as well as a comparison of specified mechanisms.

We quantify prospective precision requirements using independently varying nuisance ranges, not experimental confidence intervals. At 0 and 10 mM glutamine, the classes $\alpha_1\in[0.136,0.2125]$ and $[0.544,0.850]$, each with $K_G\in[7.8,31.2]$ mM, have a minimum log-halfpoint-ratio gap of 0.256 (Fig 10B). With known independent Gaussian measurement errors, a one-sided 5\% test has at least 80\% power with five GlnD occupancy estimates at fraction SD 0.20, or eight independent GlnE halfpoint estimates per condition at log SD 0.20 (Fig 10C). A halfpoint estimate requires a complete binding titration; these counts are not numbers of individual concentration points. Pilot measurements must establish precision and control systematic bias before these calculations are used as an experimental budget.

The guarantee is effect-size dependent. With factor-two uncertainty in $K$, single-dose GlnD occupancy ranges overlap for $\lambda\leq1.25$; additional replicates cannot remove that ambiguity (Fig 10D). Broader overlapping GlnE parameter classes similarly require additional conditions or independent constraints. Paired GlnD capacity ratios cancel a common speed but remain conditional on calibration of the baseline input-response shape. Thus the design specifies when a new observable is required, when concentration normalization suffices, and when greater precision can genuinely resolve a mechanism.

'''
CAP9=r'''\textbf{Reaction and observation assumptions define different boundaries.} (A) Explicit enzyme--target complexes preserve the free modified fraction but change total fractions; enzyme/target total is 0.3. (B) Total-fraction separation between $\lambda=1,4$ across enzyme/target ratios for one disclosed asymmetric-affinity scenario. Conservation bounds this separation by $\min(1,E_{\rm tot}/T_{\rm tot})$ in the stated monovalent model. (C) Finite regulatory binding changes the transient without changing its steady endpoint. Binding rates and time are dimensionless reference units. (D) Ligand exchange in substrate-bound complexes with state-dependent dissociation/catalysis ratios breaks free-state equivalence. This extension preserves thermodynamic binding-cycle ratios but permits catalytic regulatory currents. Its parameters and positive steady-state checks are specified in S22. All curves are calculations, not new observations.'''
CAP10=r'''\textbf{Nuisance-aware designs identify informative measurements and required precision.} (A) Target-free singly liganded GlnD fractions; shading allows $K$ to vary independently by a factor of two in each candidate. At nominal $G=K$, prediction ranges are separated by at least 0.244. (B) GlnE PII-binding halfpoints normalized to zero glutamine, with PII-UMP absent. Shading independently varies $K_G$ from 7.8 to 31.2 mM and $\alpha_1$ within the two stated effect-size classes; it is not a confidence band. (C) Independent assay estimates required for at least 80\% power at one-sided size 0.05. GlnD SD is in fraction units; GlnE SD is in log-halfpoint units. GlnE counts apply to complete titration estimates at each of two glutamine conditions. (D) Smaller GlnD alternatives require more precise measurements; overlapping nuisance ranges cannot be separated by replication alone. The Gaussian error model assumes known precision and no uncorrected condition-dependent bias. S23 gives the test and depletion correction.'''
METHOD=r'''
\subsection*{Assumption boundaries and prospective discrimination}

The finite reaction model tracks 11 concentrations and conserves enzyme and monovalent target totals. Effective affinities $h_j^\pm$ and specificities $a_j,b_j$ define positive association, dissociation, and catalytic rates. We verify constructed steady states in 300 random parameter cases, including enzyme/target ratios $10^{-4}$--1 and effective affinities $10^{-2}$--$10^2$ in reciprocal target-concentration units. Two independent integrations from free enzyme and unmodified target check convergence. A second extension permits substrate-bound ligand exchange with binding-cycle ratios derived from association/dissociation constants; positive stationary solves cover five glutamine levels and 26 exchange rates for both mechanisms, with independent integrations at a selected condition. These scenarios distinguish exact identities from topology-dependent numerical predictions.

For experimental design, nuisance ranges vary independently under the two hypotheses. GlnD occupancy bounds are analytic; GlnE ratio bounds follow the monotonic binding formula over specified parameter boxes. A normal test at the least-favourable null boundary controls the one-sided false-positive rate at 0.05. Sample sizes use known measurement SD and target power 0.80. Synthetic draws (100,000 per null and alternative for each selected assay) verify the probability calculations; they are not biological replicates. We separately report conditional GlnD capacity comparisons, for which calibration of the baseline shape is required. Full reaction rates, grids, observation definitions, design thresholds, and code are in S1 Data.

'''
ABSTRACT=r'''Sharp biochemical responses can arise from different mechanisms, limiting what aggregate titrations reveal about regulation. We analyze distributive modification cycles controlled by an effector-switched bifunctional enzyme. The shared enzyme pool cancels from adjacent free-target ratios, and local sensitivity separates into modification-ladder and effector-switch contributions. Finite GlnD binding-and-catalysis models preserve a published PII titration and every free modification-state fraction while predicting different enzyme occupancies. The published six-state GlnE topology permits a corresponding compensation between binding and catalytic synergy. Explicit reaction models show that slow binding and enzyme saturation need not remove steady-state equivalence, whereas altered catalytic coupling can do so. Conservation bounds the difference between total observations of free-state-equivalent mechanisms. We then construct complementary binding designs that account for unknown affinity and rate scales, ligand depletion, and independently varying nuisance parameters. Target-free GlnD occupancy and normalized GlnE binding titrations separate specified mechanisms; analytical noise calculations quantify the required precision and identify alternatives that replication alone cannot resolve. These results connect mechanistic ambiguity to experimentally testable conditions, while distinguishing exact identities and prospective design scenarios from physiological parameter validation.'''
AUTHOR=r'''Different enzyme mechanisms can generate the same sharp biochemical response. Measuring that response more precisely may therefore leave its mechanism unresolved. Using bacterial nitrogen signaling, we construct enzyme models that produce the same protein modification states and ask when this ambiguity survives more detailed reaction kinetics. Slow binding can change response times without changing the endpoint; counting enzyme-bound protein changes the observation; other catalytic couplings can remove the equivalence. We use these distinctions to design complementary measurements. Measuring GlnD's own ligand-bound states avoids ambiguities in its target response. Comparing GlnE binding curves across glutamine conditions removes an unknown affinity scale. The calculations specify how parameter uncertainty and measurement precision affect discrimination, and when adding another type of measurement is more useful than collecting more replicates. The proposed designs are conditional predictions that can be tested, not new experimental observations.'''
for name,embed in [('manuscript_review_with_figures',True),('manuscript_submission',False)]:
 t=(OLD/(name+'.tex')).read_text()
 t=re.sub(r'(\\section\*\{Abstract\}\s*).*?(?=\\section\*\{Author summary\})',lambda m:m.group(1)+ABSTRACT+'\n\n',t,flags=re.S)
 t=re.sub(r'(\\section\*\{Author summary\}\s*).*?(?=\\clearpage)',lambda m:m.group(1)+AUTHOR+'\n\n',t,flags=re.S)
 def fig(n,cap):
  art=rf'\includegraphics[width=\textwidth]{{figures/Fig{n}.pdf}}' if embed else ''
  return '\n'+r'\begin{figure}[!ht]\centering'+art+'\n'+r'\caption{'+cap+'}\n'+r'\end{figure}'+'\n'
 block=RESULT+fig(9,CAP9)+DESIGN+fig(10,CAP10)+'\n\\clearpage\n'
 t=t.replace(r'\section*{Discussion}',block+r'\section*{Discussion}')
 t=t.replace(r'\section*{Data and code availability}',METHOD+r'\section*{Data and code availability}')
 t=t.replace('The finite regulatory-state constructions in Figs 7 and 8 give explicit kinetic and occupancy predictions.', 'The finite regulatory-state constructions in Figs 7 and 8 give explicit kinetic and occupancy predictions. Figures 9 and 10 further distinguish an exact steady-state identity from a transient approximation, bound the total-observation gap, and translate separated nuisance ranges into prospective measurement requirements.')
 t=t.replace('A positive steady-state distribution does not establish uniqueness or dynamic stability.', 'The selected full-reaction integrations establish convergence for the tested conditions, not a general uniqueness or stability theorem.')
 t=t.replace('These tests validate the implemented calculations; they are distinct from full mass-action simulations, which are not needed to evaluate the closed-form identities presented here.', 'These tests validate the closed-form calculations; the additional full mass-action models and independent integrations are described below.')
 t=t.replace('The new matched-model and regulatory-state analyses,', 'The new matched-model, regulatory-state, and assumption-boundary/design analyses,')
 (NEW/(name+'.tex')).write_text(t)

SUPP=r'''
\clearpage
\section*{S22 Finite reaction realization, total-readout bound, and catalytic coupling}

Let $X_0,X_1$ be two free states of one modification unit. Three free enzyme states exchange glutamine with finite rates. For $j=0,1,2$ introduce
\[
 E_j+X_0\mathop{\rightleftharpoons}^{k_{\rm on,j}^+}_{k_{\rm off,j}^+} C_j^+
 \mathop{\longrightarrow}^{k_{\rm cat,j}^+} E_j+X_1,
\qquad
 E_j+X_1\mathop{\rightleftharpoons}^{k_{\rm on,j}^-}_{k_{\rm off,j}^-} C_j^-
 \mathop{\longrightarrow}^{k_{\rm cat,j}^-} E_j+X_0.
\]
All inputs are buffered free glutamine; enzyme and target totals are conserved. Set $h_j^\pm=k_{\rm on,j}^\pm/(k_{\rm off,j}^\pm+k_{\rm cat,j}^\pm)$, $a_j=h_j^+k_{\rm cat,j}^+$, and $b_j=h_j^-k_{\rm cat,j}^-$. With ligand exchange restricted to free enzyme, the steady complex equations give $C_j^+=h_j^+E_jX_0$ and $C_j^-=h_j^-E_jX_1$. Each complex's net association minus release-and-catalysis is zero, leaving the free-enzyme equations $Q_EE=0$. Any positive common multiplier of the ligand-binding generator changes its relaxation but not its stationary weights $w=(1,G/K_1,G^2/(K_1K_2))$. Consequently
\[
 \frac{X_1}{X_0}=\frac{\sum_j a_jw_j}{\sum_j b_jw_j}
\]
is exact for this full reaction system, without a rapid-equilibrium or unsaturated-target approximation. The earlier $\lambda$ transformation preserves it. This statement concerns the stationary equations, not uniqueness for every possible topology or a universal first-order time course.

Let $q=X_1/(X_0+X_1)$, $x=X_0+X_1$, and let $o_j=w_j/\sum w$. Define $H_+=\sum o_jh_j^+$, $H_-=\sum o_jh_j^-$, and $H=(1-q)H_++qH_-$. The free enzyme pool is $E_f=E_T/(1+Hx)$ and target conservation gives
\[
 Hx^2+[1+H(E_T-T_T)]x-T_T=0.
\]
Its positive root yields all 11 concentrations. The measured total modified fraction is $Y_T=(qx+E_fqxH_-)/T_T$. If $H_+=H_-$, then $Y_T=q$ exactly, including when most enzyme is complexed. Equality $h_j^+=h_j^-$ state by state is sufficient.

For general affinities write $B=\sum(C_j^++C_j^-)$ and $M_B=\sum C_j^-$. Since $0\leq M_B\leq B\leq\min(E_T,T_T)$,
\[
Y_T=q+\frac{M_B-qB}{T_T},\qquad
q-q\epsilon\leq Y_T\leq q+(1-q)\epsilon,
\quad \epsilon=\min(1,E_T/T_T).
\]
Two models with the same $q,E_T,T_T$ lie in the same interval of width $\epsilon$, proving the pairwise bound in the main text. If one enzyme can bind more than one counted modification unit, the upper bound on $B$ must be adjusted to that stoichiometry. In particular, a monovalent unit calculation must not be relabelled a fully resolved trimer sequestration model. No derivative bound follows from the concentration bound alone.

The numerical realization sets $k_{\rm cat}^+=a/h^+$, $k_{\rm cat}^-=b/h^-$, $k_{\rm off}=k_{\rm cat}$, and $k_{\rm on}=2h k_{\rm cat}$. A shared scale makes the reference summed capacity at $G=K$ equal to one. The illustrative asymmetric affinities are $h^+=(0.2,2,10)$ and $h^-=(5,0.3,0.1)$ in reciprocal target-concentration units. Target total is one; enzyme/target ratios range from $10^{-4}$ to one. Random verification uses 300 positive cases with all affinity coordinates independently spanning $10^{-2}$--$10^2$, binding-rate multipliers $10^{-3}$--$10^3$, and target totals $10^{-1}$--$10$. Full balances, conserved totals, the pairwise bound, and the symmetric-affinity identity are checked. Two BDF integrations from initially free enzyme and unmodified target converge to the independently constructed steady states.

For the finite-binding transient in Fig 9C, enzyme fractions obey $\dot o=Q_Eo$ and the dilute target obeys $\dot q=(o\cdot a)(1-q)-(o\cdot b)q$. The glutamine step is $10K\to K$, initialized at the preceding steady state. Binding multipliers span 0.01--100. These equations are an unsaturated mean-field conversion model with explicit regulatory relaxation, distinct from the 11-concentration saturation model. The units are reference times, not fitted seconds or minutes.

\textbf{Exchange in substrate-bound complexes.} Adding ligand exchange to $C^+$ and $C^-$ retains the constructed steady state if their generators annihilate the vectors $w_jh_j^+$ and $w_jh_j^-$. This is a sufficient kinetic compatibility condition and can hold for nonzero exchange rates. Thermodynamic compatibility of binding alone is weaker. Write $\rho_j=k_{\rm off,j}/k_{\rm cat,j}$; the equilibrium association constant is $k_{\rm on}/k_{\rm off}=h_j(1+\rho_j)/\rho_j$, which need not be proportional to the catalytic-complex weight $h_j$ across states.

To test this distinction without arbitrary inconsistent binding cycles, retain the same $h,a,b$ but choose $\rho^+=(0.1,1,10)$ and $\rho^-=(10,1,0.1)$. Set $k_{\rm on}=(1+\rho)h k_{\rm cat}$. For each complex chain, take reverse ligand-exchange rate $\eta$ and forward rate equal to $\eta$ times the free-chain association ratio times the adjacent ratio of $h(1+\rho)/\rho$. All noncatalytic binding cycles then satisfy the equilibrium ratio constraints. Catalysis can nevertheless sustain regulatory currents; the buffered nucleotide-driven conversion is not assumed at thermodynamic equilibrium.

Positive log-concentration stationary solves cover $G/K=0.3,1,3,10,30$, $\eta=0$ and 25 logarithmically spaced values from $10^{-4}$ to $10^2$, $E_T/T_T=0.1$, and both $\lambda$ values. The full ODE residual and both conservation laws are checked; independent BDF integrations at $G=10K,\eta=1$ verify two stationary solutions. In this scenario the maximum scanned free-fraction gaps are about 0.035 at $10K$ and 0.070 at $30K$, occurring near $\eta=0.1$. These values establish a conditional counterexample to unrestricted topology invariance; they are not estimates of biological rates or evidence that this coupling occurs in GlnD.

\clearpage
\section*{S23 Nuisance elimination, effect-size classes, and assay precision}

\textbf{GlnD occupancy.} In a target-free equilibrium binding assay,
\[
 o_1(z,\lambda)=\frac{2z/\lambda}{1+2z/\lambda+z^2},\qquad z=G/K.
\]
At nominal $G=K_0$, allow $K/K_0\in[1/2,2]$ separately under each hypothesis. For $\lambda=1$, the fraction lies in $[4/9,1/2]$; for $\lambda=4$, it lies in $[1/6,1/5]$. The minimum gap is $4/9-1/5=11/45$. This interval separation profiles the nuisance ranges independently, rather than imposing the same unknown $K$ on the two candidates. The catalytic scale and substrate affinities do not enter this target-free assay. The nominal concentration is approximately 0.0520 mM from the reference fit; the factor-two range is a design scenario, not an inferred confidence interval.

For a general alternative $\lambda\geq1$, its largest single-state occupancy over the same range is $1/(1+\lambda)$. The worst-case gap from the reference is $4/9-1/(1+\lambda)$, positive only when $\lambda>5/4$. When ranges overlap, replication of this scalar measurement alone cannot give uniform discrimination over the nuisance classes. Multiple glutamine levels or a tighter independently supported $K$ range are then required.

\textbf{GlnE binding and depletion.} Set PII-UMP to zero and keep glutamine buffered. The three PII-bound states reduce to the unmodified-PII and glutamine--PII complexes, and the PII-bound enzyme fraction is $P_f/(K_{\rm eff}+P_f)$, where
\[
 K_{\rm eff}=K_P\frac{1+G/K_G}{1+G/(\alpha_1K_G)}.
\]
If $P_T=P_f+E_T P_f/(K_{\rm eff}+P_f)$, half occupancy implies $P_T^{1/2}=K_{\rm eff}+E_T/2$. Thus a measured enzyme total corrects the nominal halfpoint before normalization across glutamine conditions. This correction assumes one PII binding unit per enzyme in the stated topology and does not correct glutamine depletion; the latter must be prevented or modelled separately.

The normalized ratio $r(G)=K_{\rm eff}(G)/K_{\rm eff}(0)$ eliminates $K_P$. For $\alpha_1\ne1$,
\[
 \frac{G}{1-r(G)}=\frac{\alpha_1K_G}{1-\alpha_1}+\frac{G}{1-\alpha_1}.
\]
Two distinct nonzero glutamine levels plus the zero-glutamine baseline determine the line's slope $m$ and intercept $b$, hence $\alpha_1=1-1/m$ and $K_G=b(1-\alpha_1)/\alpha_1$. At $\alpha_1=1$, $r=1$ and $K_G$ is not identifiable from this observation. This is a noise-free identifiability statement; noisy titration data should be fitted with the original binding model because transforming ratios changes their errors and shares the baseline error. Numerical inversion at 1 and 10 mM checks the formula across the declared parameter range.

For prospective discrimination at 0 and 10 mM, use two disjoint effect-size classes $\alpha_1\in[0.136,0.2125]$ and $[0.544,0.850]$. In both, allow $K_G\in[7.8,31.2]$ mM independently. Monotonicity of the ratio gives log-ratio intervals $[-1.51931,-0.64158]$ and $[-0.38589,-0.04194]$. The minimum gap is 0.25570. These classes specify the alternatives being tested; neither is a confidence region from the published experiments. If each nominal $\alpha_1$ instead spans a factor of two, the classes touch and there is no strictly positive worst-case guarantee. The three-titration design then provides an alternative to an unjustified binary claim.

\textbf{Statistical assumptions and test.} Let each occupancy estimate have independent normal error with known SD $\sigma$, or let each independent halfpoint estimate have normal log error with known SD $\sigma$. Background-corrected fraction estimates are not clipped in this error model. Precision must be established in pilot measurements; the calculations do not include uncertain error variance, condition-dependent calibration bias, or fitted biological variability. A single GlnE halfpoint estimate is obtained from an entire titration, not one concentration point. For $n$ independent estimates per condition, the occupancy mean has variance $\sigma^2/n$, while the difference of two mean log halfpoints has variance $2\sigma^2/n$.

Use the least-favourable null boundary and the prespecified sign of the alternative. A one-sided size-0.05 normal test has worst-case power
\[
 \Phi\!\left(\frac{\sqrt n\,\Delta}{\sqrt v\,\sigma}-z_{0.95}\right),
 \qquad n\geq\frac{v\sigma^2(z_{0.95}+z_{0.80})^2}{\Delta^2},
\]
rounded upward. Applied to total-readout separation $\Delta\leq\epsilon$, it gives an optimistic necessary sample size $n\geq\sigma^2(z_{0.95}+z_{0.80})^2/\epsilon^2$ for one normally distributed fraction observable and two specified alternatives. With $\epsilon=0.01,\sigma=0.05$, this rounds up to 155; smaller actual separation or nuisance uncertainty can only worsen this bound. For the binding designs use $v=1$ for occupancy and $v=2$ for paired log halfpoints. Here $\Delta$ is the independently profiled interval gap, not a distance between two selected curves. At SD 0.20, five occupancy estimates yield worst-case power 0.862, and eight halfpoint estimates per condition yield 0.819. At occupancy SD 0.10, two estimates yield 0.965. The selected tests are independently verified by 100,000 synthetic null and alternative draws. These are prospective assay-error calculations, not experimentally established budgets. If independent full titrations are infeasible, a hierarchical joint fit needs a separate power analysis.

\textbf{Conditional kinetic comparison.} The capacity multiplier is $R_\lambda(z)=(1+2z+z^2)/(1+2z/\lambda+z^2)$. An unknown speed absorbs a single-condition amplitude. Paired capacities cancel a shared speed, but reciprocal doses are uninformative because $R_\lambda(z)=R_\lambda(1/z)$. Searching 101 doses over $0.02K_0$--$50K_0$ selects approximately $0.020K_0$ and $0.855K_0$, with minimum log-multiplier contrast 0.363 over factor-two $K$ variation. This contrast requires independently calibrated baseline capacities and $K$; independently refitting the baseline shape can restore overlap. The data label its sample sizes conditional, and Fig 10 excludes them from the profiled guarantees.

'''
t=(OLD/'S1_Appendix.tex').read_text();idx=t.index(r'\begin{figure}[p]');t=t[:idx]+SUPP+'\n\\clearpage\n'+t[idx:];(NEW/'S1_Appendix.tex').write_text(t)
cl=(OLD/'cover_letter.tex').read_text();idx=cl.index('Joint cascade comparisons')
cl=cl[:idx]+'''The study connects these identities to assay design. Explicit finite-rate enzyme--target models delineate when saturation preserves free-state equivalence, conservation bounds total-readout differences, and a substrate-bound regulatory extension shows when the identity can break. Complementary GlnD and GlnE binding designs remove nuisance scales or retain separated prediction ranges. Prospective error calculations specify precision requirements and identify cases that replication alone cannot resolve. The parameter ranges and assay-error assumptions are explicit scenarios, not physiological validation.\n\n'''+cl[idx:];(NEW/'cover_letter.tex').write_text(cl)
print(NEW)
