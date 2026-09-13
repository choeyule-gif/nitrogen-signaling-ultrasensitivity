"""Workspace authoring script; portable analyses do not depend on this file."""
from pathlib import Path
import shutil,json,csv,re

BASE=Path(__file__).resolve().parents[3]
ROOT=Path(__file__).resolve().parents[1]
OLD=BASE/'output/pdf/plos_model_comparison'
NEW=BASE/'output/pdf/plos_biochemical_revision';NEW.mkdir(exist_ok=True)
R=json.loads((ROOT/'results/summary.json').read_text())
for ext in ['*.sty','*.bst']:
    for p in OLD.glob(ext):shutil.copy2(p,NEW/p.name)
shutil.copytree(OLD/'figures',NEW/'figures',dirs_exist_ok=True)
for ext in ['pdf','png']:shutil.copy2(ROOT/f'figures/Fig7_finite_mechanisms.{ext}',NEW/f'figures/Fig7.{ext}')

CAPTION=r'''\textbf{Finite reaction models preserve ambiguity while kinetic constraints expose its scope.} (A) The independent-binding, sequential-activity model and its refractory counterpart have identical mean curves (solid and dashed). Points are the twelve digitized coordinates. The dotted curve imposes a UR activity half-range of 0.080 mM in a different, cooperative-binding realization. (B) Predicted state fractions at 10 mM for the independent-binding pair. (C) Monotone UT and UR capacities of that pair; UT curves are normalized to their zero-input limits and UR to its high-input limit. These are predictions, not additional kinetic observations. (D) Separate-site saturation breaks equivalence when effective affinities are held fixed and one gain is calibrated at 0.53 mM. Lines use $K_Q=0.28\,\mu$M; bands span 0.28--0.85 $\mu$M. The ordinate is the maximum absolute difference between unnormalized GS fractions over 0.02--10 mM. Rescaling the unknown unmodified-PII affinity restores equivalence (dashed baseline). Apparent-constant ranges motivate scenarios, not experimental uncertainty intervals. All calculations use the specified free-input and count-weighted reverse-recognition assumptions.'''

SECTION=r'''
\subsection*{A finite biochemical realization separates structural ambiguity from kinetic plausibility}

The population ambiguity does not require an arbitrary switch function. Consider three free enzyme states linked by $E_0+G\rightleftharpoons E_1$ and $E_1+G\rightleftharpoons E_2$. Their weights are $w=(1,G/K_1,G^2/(K_1K_2))$. State-specific UT and UR capacities, $\alpha_j$ and $\beta_j$, give $A(G)=\sum_j\alpha_jw_j$ and $B(G)=\sum_j\beta_jw_j$, with homogeneous mean $\theta=A/(A+B)$. These rates admit an explicit finite enzyme--target-complex realization for the free-target steady state (S18). They do not assign the two binding events to particular ACT domains.

For any $r\le\min_j\{\alpha_j/(\alpha_j+\beta_j)\}$, replacing the accessible population's capacities by
\[
\alpha_j^{\mathrm R}=(1-r)\alpha_j-r\beta_j,\qquad
\beta_j^{\mathrm R}=\beta_j
\]
preserves nonnegative rates and yields accessible mean $q=(\theta-r)/(1-r)$. Adding a fraction $r$ of fully modified, UR-inaccessible target recovers exactly the original mean. Thus the same finite ligand-binding network realizes both population models; only its catalytic coefficients and population composition differ. UR-inaccessibility is a model hypothesis, not a measured PII subpopulation.

An illustrative implementation uses two identical independent ligand-binding sites, $K_1=K/2$ and $K_2=2K$, with $\alpha_1=\alpha_2$ and $\beta_1=\beta_0$. UT decreases monotonically and UR increases monotonically. Fixing the UR specificity-fold scenario to $(6.5/0.82)/(2.7/2.3)=6.752$ from the rounded Mg-dependent kinetic entries in \cite{r17}, three fitted parameters give an RMSE of 0.0366 UMP groups per trimer. The inferred common-site $K$ is approximately 0.052 mM. The two population models predict approximately 0.004 and 0.146 fully modified trimers at 10 mM despite identical mean curves (Fig 7A--C). Independence of ligand binding in this construction does not imply independence of its state-dependent catalytic activities.

Finite realization is not kinetic validation. The predicted UT and UR half-range concentrations are approximately 0.022 and 0.126 mM, whereas the older study reports UT inhibition near 0.070--0.080 mM and UR activation near 0.080 mM. A second realization can impose the 0.080 mM UR half-range together with the 6.752-fold specificity change and retain an RMSE of 0.0378, but it requires strongly cooperative binding. Further imposing the older unliganded UT/UR specificity ratio increases that realization's RMSE to 0.540 (S19). Different assay conditions and the distinction between an apparent activity constant and a microscopic binding constant preclude treating these transfers as a joint kinetic fit. The successful partial transfers and failed stricter transfer are both reported.

The downstream ambiguity likewise need not rely on a very large gain. Set $P=T_{\rm PII}p_0$ and $Q=T_{\rm PII}\theta$. For two saturating regulatory inputs, a candidate ratio is $v=k[P/(K_P+P)]/[Q/(K_Q+Q)]$. Holding $K_P$ fixed generally breaks the gain-only equivalence. However, because $P_{\mathrm R}=P_{\mathrm H}/(1-r)^2$ and $Q_{\mathrm R}=Q_{\mathrm H}$, the single shared transformation
\[
K_{P,\mathrm R}=K_{P,\mathrm H}/(1-r)^2
\]
preserves both occupancies at every input and target concentration, without changing $k$. The fitted finite pair requires an approximately 37\% change in $K_P$. A finite receptor with separate binding sites realizes these occupancy classes when reverse recognition is proportional to modification count (S19). This also permits common glutamine regulation of receptor states. State-specific reverse recognition or appreciable depletion of free target can break the construction.

The measured apparent PII activation constants vary with glutamine and other regulators \cite{r71}; they do not independently fix the microscopic $K_P$ in the matched cascade assays. We therefore report fixed-affinity and rescaled-affinity scenarios separately (Fig 7D), without declaring either model physiologically established. The result replaces an unconstrained-gain argument with a specific test: measure unmodified-PII binding under the same regulatory conditions, together with modification-state fractions.

'''

METHOD=r'''
\subsection*{Finite reactions and biochemical constraint audit}

The independent-binding model fits two plateau fractions and one ligand scale, fixes $\alpha_1=\alpha_2$ and $\beta_1=\beta_0$, and uses the transported UR specificity-fold scenario described in S18. Three starts and leave-one-coordinate-out refitting assess numerical and predictive behavior. A separate positive rational family profiles effective coefficient ratios; these are not microscopic binding ratios when total catalytic capacities differ between enzyme states. The cooperative realization has two fitted parameters after fixing its effective coefficient ratio and transporting the UR fold and half-range. A stricter transfer also fixes the unliganded UT/UR specificity ratio. These scenarios are compared by coordinate RMSE, without assigning experimental likelihoods to cross-assay constraints. Primary-source tables, assay contexts, and transcription corrections are documented in S19 and S1 Data.

Finite reaction verification evaluates every enzyme, target, and complex balance at ten constructed free-input steady states. Saturation scenarios use 17 forward scales between 0.015 and 2.5 $\mu$M, 13 reverse scales between 0.28 and 0.85 $\mu$M, three target totals, and 401 glutamine inputs. Each fixed-affinity comparison matches one unnormalized GS response at 0.53 mM; it does not refit all three experimental cascade midpoints. The affinity-rescaling identity is verified for both the saturated input ratio and the receptor occupancy classes. Neither these scenarios nor the explicit complex verification constitutes a closed-total fit of the experimental cascade.

'''

for name in ['manuscript_review_with_figures.tex','manuscript_submission.tex']:
    t=(OLD/name).read_text()
    fig='\\begin{figure}[!ht]\n\\centering\n'
    if 'review' in name:fig+=r'\includegraphics[width=\textwidth]{figures/Fig7.pdf}'+'\n'
    fig+='\\begin{singlespace}\n\\caption{'+CAPTION+'}\n\\end{singlespace}\n\\label{fig7}\n\\end{figure}\n\n'
    marker=r'\subsection*{State measurements target the ambiguity left by mean titrations}'
    assert t.count(marker)==1;t=t.replace(marker,SECTION+fig+marker)
    old='We further show that a refractory PII fraction can preserve both the mean titration and the complete downstream response when an unknown selective-readout gain is calibrated.'
    new='A finite ligand-binding and catalytic-state construction realizes homogeneous and refractory populations with identical means. Their downstream equivalence can survive saturating recognition when an unknown PII affinity is rescaled; independent kinetic constraints can break it.'
    assert old in t;t=t.replace(old,new)
    t=t.replace('These equivalent models predict approximately 0.4\\% versus 15\\% fully uridylylated trimers at 10 mM glutamine.','Finite independent-binding alternatives predict approximately 0.3\\% versus 15\\% fully uridylylated trimers at 10 mM glutamine. Cross-assay kinetic transfers expose remaining limits to their biological interpretation.')
    t=t.replace('while its biological realizability and the fixed candidate shapes remain assumptions.','and a finite biochemical construction removes the need for an arbitrary switch in the population counterexample. Physiological parameter values, free-to-total observation mapping, and the fixed candidate shapes remain assumptions.')
    t=t.replace('The shared-parameter fits reported here remain reduced descriptions, not fully reconstituted mass-action models with independently measured microscopic rates.','The shared-parameter fits in Fig 6 remain reduced descriptions, not fully reconstituted models with independently measured microscopic rates. The separate finite construction in Fig 7 establishes realizability under a specified recognition law, but its failed strict kinetic transfer prevents interpreting it as a validated GlnD mechanism.')
    t=t.replace('Allowing the shared upstream ladder parameter to vary provides a stricter robustness check (Fig 6D).','Allowing the shared upstream ladder parameter to vary provides a mathematical sensitivity analysis (Fig 6D), rather than a biochemical plausibility test.')
    t=t.replace(r'\section*{Data and code availability}',METHOD+r'\section*{Data and code availability}')
    t=t.replace('A second realization can impose the 0.080 mM UR half-range', 'Its unliganded UT/UR specificity ratio is approximately 905, compared with 38.9 from the older kinetic entries. A second realization can impose the 0.080 mM UR half-range')
    t=t.replace('The successful partial transfers and failed stricter transfer are both reported.', 'The successful partial transfers and failed stricter transfer are both reported. Neither finite pair is a joint fit to all reported cascade coefficients.')
    t=t.replace('approximately 0.004 and 0.146 fully modified trimers', 'approximately 0.0034 and 0.146 fully modified trimers')
    abstract=r'''Sharp biochemical responses can arise from different mechanisms, limiting what aggregate titrations reveal about regulation. We analyze distributive modification cycles controlled by an effector-switched bifunctional enzyme. At steady state, the shared enzyme pool cancels from adjacent free-target ratios, and local sensitivity separates into modification-ladder and effector-switch contributions. Identical independent sites contribute a ladder factor of one regardless of site number. We construct competing explanations of published bacterial nitrogen-signaling measurements under common observation rules. A finite ligand-binding and catalytic-state model fits a twelve-point PII titration with an error of 0.0366 UMP groups per trimer and admits an exactly mean-matched refractory population. Independent ligand binding suffices for this construction. The two models predict approximately 0.3\% versus 15\% fully uridylylated trimers at 10 mM glutamine. Their downstream equivalence survives saturating recognition when an unknown PII affinity is rescaled, but independently fixed affinities or state-specific recognition can break it. Cross-assay kinetic transfers expose limits to physiological interpretation rather than validating either mechanism. Broader cascade comparisons show that amplification deficits and inhibition requirements depend on model assumptions. Modification-state and matched-condition binding measurements therefore target ambiguities that additional mean-response measurements can leave unresolved.'''
    t=re.sub(r'(\\section\*\{Abstract\}\s*).*?(\s*\\section\*\{Author summary\})',lambda m:m.group(1)+abstract+'\n\n'+m.group(2),t,flags=re.S)
    t=t.replace(r'\end{thebibliography}',r'\bibitem{r71} Jiang P, Mayo AE, Ninfa AJ. Escherichia coli glutamine synthetase adenylyltransferase (ATase, EC 2.7.7.49): kinetic characterization of regulation by PII, PII-UMP, glutamine, and alpha-ketoglutarate. Biochemistry. 2007;46:4133--4146. \href{https://doi.org/10.1021/bi0620510}{doi:10.1021/bi0620510}.'+'\n\n'+r'\end{thebibliography}')
    (NEW/name).write_text(t)

SUPP=r'''
\section*{S18 Finite ligand and catalytic states realizing the population ambiguity}

All quantities in this section refer to free target and buffered free glutamine unless otherwise stated. Three enzyme states have binding weights $w=(1,G/K_1,G^2/(K_1K_2))$ and nonnegative forward and reverse catalytic specificities $\alpha_j,\beta_j$. Their effective forward/reverse ratio is $A/B$, where $A=\sum_j\alpha_jw_j$ and $B=\sum_j\beta_jw_j$. Independent modification sites yield $\theta=A/(A+B)$ and binomial free-target state fractions.

One explicit finite realization uses $E_0+G\rightleftharpoons E_1$ and $E_1+G\rightleftharpoons E_2$, with ligand binding only to the free enzyme states. For each $j=0,1,2$ and $i=0,1,2$, include
\[
E_j+S_i\rightleftharpoons B_{ji}\longrightarrow E_j+S_{i+1},\qquad
E_j+S_{i+1}\rightleftharpoons C_{ji}\longrightarrow E_j+S_i.
\]
Choose forward specificities $(3-i)\alpha_j$ and reverse specificities $(i+1)\beta_j$. At steady state, eliminating each complex cancels its contribution to the free-enzyme balance. The ligand-binding ratios therefore retain weights $w$, and the target fluxes give the binomial ladder with input $A/B$. For example, in consistent nondimensional units, $k_{\rm off}=k_{\rm cat}=1$ and $k_{\rm on}=2$ times the desired specificity provide such a realization. Buffered nucleotide cofactors are absorbed into catalytic rates. This establishes the free-target result with explicit complexes; it does not identify closed-system total modification with free modification. The supplied verification constructs ten positive/nonnegative steady states and checks every species balance. Conserved totals are computed for each state and are not claimed to remain fixed across these constructed glutamine inputs.

Let $r\le\min_j \alpha_j/(\alpha_j+\beta_j)$ and set $\alpha_j^{\rm R}=(1-r)\alpha_j-r\beta_j$ and $\beta_j^{\rm R}=\beta_j$. Then $A_{\rm R}=(1-r)A-rB$ and $A_{\rm R}+B=(1-r)(A+B)$, so $q=A_{\rm R}/(A_{\rm R}+B)=(\theta-r)/(1-r)$. A separate fully modified species with free fraction $r$ and no GlnD-mediated turnover gives the same total free mean. The refractory species is assumed available to the stated downstream recognition law. It is neither an experimentally established species nor a dynamically exchanging state. Vanishing catalytic coefficients denote omitted reactions. For $r=\alpha_2/(\alpha_2+\beta_2)$, the accessible model has $\alpha_2^{\rm R}=0$.

For the independent-binding implementation, set $K_1=K/2$, $K_2=2K$, $\alpha_1=\alpha_2$ and $\beta_1=\beta_0$. Thus
\[
U(G)=\alpha_2+\frac{\alpha_0-\alpha_2}{(1+G/K)^2},\qquad
R(G)=\beta_0+(\beta_2-\beta_0)\left(\frac{G/K}{1+G/K}\right)^2.
\]
These capacities are monotone for $\alpha_0\ge\alpha_2$ and $\beta_2\ge\beta_0$. The refractory transformation preserves monotonicity because it subtracts an increasing reverse capacity from a decreasing forward capacity. The model has independent ligand binding with state-dependent catalytic gating; it does not claim that enzyme activities respond independently to each binding event.

The fit parameters are $q_0$, $q_2$, and $K$, where $q_j=\alpha_j/(\alpha_j+\beta_j)$. Fix the arbitrary common scale by $\alpha_0+\beta_0=1$. With transported specificity fold $F=6.752484$, set $d_2=F(1-q_0)/(1-q_2)$, $\alpha_2=d_2q_2$, and $\beta_2=d_2(1-q_2)$. Bounds are $0.9\le q_0\le0.999999$, $10^{-6}\le q_2\le0.499999$, and $0.001\le K\le10$ mM. Three initializations are supplied. Every leave-one-coordinate-out fit is independently optimized. Outputs report all fitted coefficients and errors. No microscopic binding data were included in the objective.

For comparison, positive rational means are parameterized as
\[
\theta(G)=\frac{q_0+q_1G/(s\sqrt C)+q_2(G/s)^2}{1+G/(s\sqrt C)+(G/s)^2}.
\]
Here $C$ is an effective coefficient ratio, not generally $K_1/K_2$. A profile restricts $C$ to upper bounds 0.25, 1, 10, 100, and 10000 and fits ordered $q_0\ge q_1\ge q_2$. The upper boundary is reached in every scenario; this is evidence of limited identification within that chosen family, not an estimate of binding cooperativity. For a second, fixed-$C=100$ construction, set $q_1=q_0$, fix $F$ as above, and impose activity half-range $h=0.080$ mM using
\[
d_2=\frac{h^2}{s^2+sh/\sqrt C},\qquad q_0=1-d_2(1-q_2)/F.
\]
Only $q_2$ and $s$ are fitted. Choose $d_1=1$; the actual binding constants are $K_1=s\sqrt C$ and $K_2=sd_2/\sqrt C$. Their ratio is approximately 4500 in this transferred-kinetics construction. This is a strong, unmeasured cooperativity requirement. The independent-binding fit demonstrates that such a requirement is not inherent to the population equivalence.

\section*{S19 What biochemical measurements constrain, and what remains unidentified}

Jiang, Peliska and Ninfa (1998), Table 2, reports Mg-dependent UT glutamine inhibition at 70--80 $\mu$M and UR activation at 80 $\mu$M. The Mg-dependent UR specificity ratio calculated from the rounded $k_{\rm cat}$ and $K_m$ entries is $(6.5/0.82)/(2.7/2.3)=6.752484$. The older manuscript's rounded specificity entries 7.93 and 1.17 give 6.778 instead; this small arithmetic difference is retained explicitly rather than treated as a new measurement. The Mn-dependent entries are UT inhibition 150 $\mu$M and UR activation 0.70 mM, with the footnoted cofactor restrictions. The repository had these two Mn values transposed and compensated by reversing the ratio in its original plotting script. Both labels and the ratio direction are corrected together; the numerical ceiling $2\sqrt{700/150}/(1+\sqrt{700/150})$ and the current manuscript's interpretation are unchanged.

These measurements are apparent activity constants and catalytic summaries under their stated substrate/cofactor conditions. They are not measurements of the microscopic sequential ligand-binding constants in S18. In the independent-binding fit, transporting only the UR specificity fold gives an RMSE of 0.0366, but UT and UR half-ranges of approximately 0.022 and 0.126 mM. In the fixed-$C$ alternative, transporting the UR fold and 0.080 mM half-range gives an RMSE of 0.0378. A stricter test also transports the unliganded specificity ratio $(137/3)/(2.7/2.3)=38.9012$, fixing $q_0=0.974938$. Refitting the remaining parameter in that fixed-$C$ topology gives an RMSE of 0.540. Thus the partial success must not be presented as agreement with the full enzyme-kinetic dataset. Nor does failure of this cross-assay transfer rule out every finite model.

Jiang, Mayo and Ninfa (2007), Tables 2--3, reports PII apparent activation constants from 0.015 to 2.5 $\mu$M in the selected rows without PII-UMP, as glutamine and other conditions change. Rows including PII-UMP reach approximately 6 and greater than 40 $\mu$M. Reported PII-UMP activation constants include 0.28, 0.35, and 0.85 $\mu$M under different conditions. These motivate a broad scenario envelope; pooling them does not produce an experimentally justified prior or confidence interval for $K_P$ or $K_Q$. The data file preserves selected row identifiers, units, and regulator conditions. The primary paper supports distinct PII and PII-UMP regulatory interactions and condition-dependent apparent constants. Our independent receptor sites are a simplifying candidate, not a reconstruction of all those interactions.

Write $P=T p_0$ and $Q=T\theta$. For a finite receptor with separate regulatory sites, the four aggregated occupancy classes have weights $(1,P/K_P,Q/K_Q,PQ/(K_PK_Q))$. The reverse input $Q$ is shorthand for three competing ligands $S_1,S_2,S_3$ with dissociation constants $3K_Q/i$ for state $i$ and equal effects after binding. The underlying receptor thus has eight microstates (forward site empty/bound, reverse site empty/bound to one of three states). This is a count-weighted recognition assumption. Adding identical glutamine regulation on both sides of the comparison preserves the following identity as long as the normalized input weights are preserved; it does not require PII-independent AT activity to be zero.

The population construction gives $P_{\rm R}=cP_{\rm H}$ with $c=(1-r)^{-2}$ and leaves $Q$ unchanged. Setting $K_{P,\rm R}=cK_{P,\rm H}$ preserves every occupancy class. More generally, multiplying every dissociation constant for unmodified PII by $c$ preserves the corresponding binding weights, including condition-specific ones. In the independent-binding fitted pair, $c$ is approximately 1.37. This is one shared transformation over all glutamine and target levels, not a separate calibration at each data point. It establishes an identifiable combination of population composition and affinity under the chosen observation model. It does not demonstrate that either parameter set is biologically correct.

If $K_P$ is known and held fixed, a separate-site saturated ratio $v=k[P/(K_P+P)]/[Q/(K_Q+Q)]$ is generally not preserved by a constant gain. The supplied 663 scenarios use three PII totals, 17 $K_P$ values, and 13 $K_Q$ values. At each setting, the homogeneous gain makes $Y=v/(1+v)=1/2$ at 0.53 mM and the refractory gain matches that response. The output is the maximum absolute GS-fraction difference across 401 inputs from 0.02 to 10 mM; it is not a response-range coefficient or a probability of experimental discrimination. Allowing the shared affinity transformation restores the full curve to numerical precision. State-specific reverse recognition, ligand/target sequestration, or independently fixed affinities can invalidate this result. No closed-total GlnE/PII conservation fit is claimed.

The measurements that would close the remaining gap are: (i) UT and UR initial-rate curves in the same glutamine, 2-OG, ATP/UTP, metal-ion, and target conditions as the mean titration; (ii) PII state fractions, including the high-input fully modified fraction; and (iii) state-resolved GlnE binding or activity measurements at those same conditions. The new construction makes their role explicit without assuming those missing measurements exist.

'''
SUPP=SUPP.replace('The repository had these two Mn values transposed and compensated by reversing the ratio in its original plotting script. Both labels and the ratio direction are corrected together; the numerical ceiling $2\\sqrt{700/150}/(1+\\sqrt{700/150})$ and the current manuscript\'s interpretation are unchanged.', 'Correcting transposed Mn metadata and the compensating plotting ratio preserves the numerical ceiling; S1 Data records both corrections.')
SUPP=SUPP.replace('The data file preserves selected row identifiers, units, and regulator conditions. The primary paper supports distinct PII and PII-UMP regulatory interactions and condition-dependent apparent constants. Our independent receptor sites are a simplifying candidate, not a reconstruction of all those interactions.', 'Source-row identifiers and conditions are retained in S1 Data. Independent receptor sites are a candidate simplification.')
SUPP=SUPP.replace('The measurements that would close the remaining gap are: (i) UT and UR initial-rate curves in the same glutamine, 2-OG, ATP/UTP, metal-ion, and target conditions as the mean titration; (ii) PII state fractions, including the high-input fully modified fraction; and (iii) state-resolved GlnE binding or activity measurements at those same conditions. The new construction makes their role explicit without assuming those missing measurements exist.', 'Discriminating measurements are matched-condition UT/UR initial rates, high-input PII state fractions, and state-resolved GlnE binding. These remain prospective tests.')
t=(OLD/'S1_Appendix.tex').read_text();marker=r'\begin{figure}[p]'
idx=t.index(marker);t=t[:idx]+SUPP+'\n\\clearpage\n'+t[idx:]
(NEW/'S1_Appendix.tex').write_text(t)

cl=(OLD/'cover_letter.tex').read_text()
cl=cl.replace('An exact population equivalence further demonstrates that a refractory PII fraction can preserve both the mean titration and the complete downstream response when an unknown selective-readout gain is calibrated.', 'A finite ligand-binding and catalytic-state construction realizes an exact population equivalence. It can persist downstream under saturating recognition when an unknown PII affinity is rescaled. Biochemical constraint tests distinguish this structural result from a validated kinetic mechanism.')
(NEW/'cover_letter.tex').write_text(cl)

rows=[
 ['10.1021/bi980667m','Table 2','Mg UT','Ki glutamine','70--80','uM','Mg-dependent; see source substrate conditions','apparent inhibition; not microscopic Kd'],
 ['10.1021/bi980667m','Table 2','Mg UR','Kact glutamine','80','uM','Mg-dependent; see source substrate conditions','apparent activation; cross-assay half-range scenario'],
 ['10.1021/bi980667m','Table 2','Mn UT','Ki glutamine','150','uM','No 2-oxoglutarate','corrected transposed repository entry'],
 ['10.1021/bi980667m','Table 2','Mn UR','Kact glutamine','700','uM','No 2-oxoglutarate or ATP','corrected transposed repository entry'],
 ['10.1021/bi980667m','Table 2','Mg UT','kcat; Km PII','137; 3.0','min^-1; uM','Different kinetic determinations in the source','used only in strict cross-assay scenario'],
 ['10.1021/bi980667m','Table 2','Mg UR no Gln','kcat; Km PII-UMP','2.7; 2.3','min^-1; uM','No glutamine','rounded inputs for specificity fold'],
 ['10.1021/bi980667m','Table 2','Mg UR with Gln','kcat; Km PII-UMP','6.5; 0.82','min^-1; uM','2.5 mM glutamine','rounded inputs for specificity fold'],
 ['10.1021/bi0620510','Table 2: 010605','PII activation','apparent Kact','1.00','uM','2-OG 0.05 mM; Gln absent; ATase 0.10 uM; PII-UMP absent','scenario input; not microscopic Kd'],
 ['10.1021/bi0620510','Table 2: 010705','PII activation','apparent Kact','2.50','uM','2-OG 1 mM; Gln absent; ATase 1.20 uM; PII-UMP absent','scenario upper input; not a universal bound'],
 ['10.1021/bi0620510','Table 2: 011005','PII activation','apparent Kact','0.015','uM','2-OG 0.05 mM; Gln 50 mM; ATase 0.01 uM; PII-UMP absent','scenario lower input; not a universal bound'],
 ['10.1021/bi0620510','Table 2: 011005','PII activation','apparent Kact','0.024','uM','2-OG 1 mM; Gln 50 mM; ATase 0.01 uM; PII-UMP absent','condition dependence'],
 ['10.1021/bi0620510','Table 2: 030705','PII activation','apparent Kact','>40','uM','2-OG 1 mM; Gln absent; ATase 2 uM; PII-UMP 10 uM','excluded from no-PII-UMP scenario envelope'],
 ['10.1021/bi0620510','Table 3: 022305','PII-UMP activation','apparent Kact','0.28','uM','2-OG 1 mM; PII absent; Gln absent; ATase 0.02 uM','scenario input; not microscopic Kd'],
 ['10.1021/bi0620510','Table 3: 022305','PII-UMP activation','apparent Kact','0.35','uM','2-OG 0.02 mM; PII absent; Gln absent; ATase 0.02 uM','condition dependence'],
 ['10.1021/bi0620510','Table 3: 042105','PII-UMP activation','apparent Kact','0.85','uM','2-OG 1 mM; PII absent; Gln 50 mM; ATase 0.03 uM','scenario input; not microscopic Kd'],
]
with (ROOT/'data/primary_source_constraints.csv').open('w') as f:
    w=csv.writer(f);w.writerow(['doi','source_locator','activity','quantity','value','units','conditions','use_and_limit']);w.writerows(rows)
print(NEW)
