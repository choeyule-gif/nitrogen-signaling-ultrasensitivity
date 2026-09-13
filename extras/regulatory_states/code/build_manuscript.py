"""Author the integrated regulatory-state revision from the preceding draft."""
from pathlib import Path
import re,shutil,json
BASE=Path(__file__).resolve().parents[3];ROOT=Path(__file__).resolve().parents[1]
OLD=BASE/'output/pdf/plos_biochemical_revision';NEW=BASE/'output/pdf/plos_regulatory_states'
NEW.mkdir(exist_ok=True)
for glob in ['*.sty','*.bst']:
    for p in OLD.glob(glob):shutil.copy2(p,NEW/p.name)
shutil.copytree(OLD/'figures',NEW/'figures',dirs_exist_ok=True)
shutil.copy2(OLD/'figures/Fig7.pdf',NEW/'figures/FigS10_finite_population.pdf')
for n,src in [(7,'Fig7_GlnD_states'),(8,'Fig8_GlnE_states')]:
    for ext in ['pdf','png']:shutil.copy2(ROOT/f'figures/{src}.{ext}',NEW/f'figures/Fig{n}.{ext}')

GLND=r'''
\subsection*{GlnD catalytic states distinguish ligand binding from modification output}

GlnD's separate NT and HD catalytic domains and its ACT-dependent glutamine response motivate a finite regulatory-state description \cite{r18}. The domain experiments support distinct catalytic and sensing functions; they do not identify two independent glutamine-binding sites or a particular sequential activation order. We therefore ask which features of a finite binding-and-catalysis model the titration can identify, without assigning the mathematical states to individual ACT domains.

Let three free enzyme states have weights $w=(1,G/K_1,G^2/(K_1K_2))$ and catalytic specificities $a_j,b_j$ for UT and UR. The independent-target steady state depends on $A/B$, where $A=\sum_j a_jw_j$ and $B=\sum_j b_jw_j$. These specificities can be realized by explicit enzyme--target complexes for the free-target steady state (S18). A reference model with two identical independent ligand sites, $K_1=K/2$ and $K_2=2K$, sets $a_1=a_2$ and $b_1=b_0$. The three-parameter fit to the twelve coordinates gives an RMSE of 0.0366 UMP groups per trimer and $K=0.0520$ mM. The UR specificity-fold scenario is retained from S18; the coefficients are not independently measured microscopic rates.

Now change the first and second binding constants in opposite directions and rescale the singly liganded state's two catalytic coefficients:
\[
K_1\mapsto\lambda K_1,\qquad K_2\mapsto K_2/\lambda,
\qquad (a_1,b_1)\mapsto\lambda(a_1,b_1).
\]
The singly liganded weight decreases by $\lambda$, so both $A$ and $B$ remain exactly unchanged. Consequently, every PII modification-state fraction, not only its mean, is preserved. For $\lambda=1,2,4$, the binding ratios $K_1/K_2$ are 0.25, 1, and 4. All three constructions retain monotonically decreasing UT and increasing UR capacities (Fig 7A; S20). They differ in enzyme occupancy and catalytic speed, rather than requiring a refractory target population.

At $G=K$, the singly liganded enzyme fractions are $1/2$, $1/3$, and $1/5$. In the rapid-binding, unsaturated conversion model, both opposing capacities are multiplied by the same factor, giving PII relaxation-rate ratios of 1, $4/3$, and $8/5$ (Fig 7B--D). The mean number of bound ligands is nevertheless one at this input for all three models. An assay of average ligand occupancy should therefore sample concentrations flanking $K$; resolving singly liganded enzyme or measuring relaxation can use the central region directly. This distinction prevents confusing an informative enzyme-state measurement with an uninformative average at a symmetry point.

'''

GLNE=r'''
\subsection*{A published GlnE regulatory model separates equilibrium output from catalytic speed}

For GlnE, we use the six-state regulatory topology proposed by Jiang, Mayo and Ninfa \cite{r71}, rather than introducing an unconstrained downstream switch. Its states are free enzyme and complexes with glutamine, PII, PII-UMP, glutamine plus PII, and PII plus PII-UMP. With $g=G/K_G$, $p=P/K_P$, and $u=U/K_U$, their common binding polynomial is
\[
Z_E=1+g+p+u+pg/\alpha_1+pu/\alpha_2.
\]
Here $\alpha_1$ and $\alpha_2$ describe the binding interactions of PII with glutamine and PII-UMP. The three AT-active states and the single AR-active state give
\[
v_{\rm AT}=k_P\frac{p+\beta g+\beta'pg/\alpha_1}{Z_E},\qquad
v_{\rm AR}=k_R\frac{u}{Z_E}.
\]
The term $\beta g$ explicitly retains AT activity without unmodified PII, as observed experimentally. The source topology and its regulatory interpretation are retained; the illustrative parameter set is not a global reconstitution of every historical assay.

The steady-state drive is therefore
\[
v=\frac{k_P}{k_R}\frac{p+\beta g+\beta'pg/\alpha_1}{u}.
\]
Two mechanistic ambiguities follow. First, $\alpha_2$ changes both rates but cancels from their ratio. Second, scaling $\alpha_1$ and $\beta'$ by the same positive factor preserves the ratio: weaker binding synergy can be compensated by greater catalytic activity of the doubly occupied state. These conclusions apply at every free input, without a new gain for each condition. Coupling this GlnE model to any member of the GlnD family therefore gives the same steady GS curve under the stated conversion law (Fig 8A).

The equivalent mechanisms predict distinct PII-binding curves and relaxation rates. At $G=10$ mM and $U=0.3\,\mu$M, multiplying $\alpha_1$ and $\beta'$ by four shifts the predicted PII-binding halfpoint from approximately 0.090 to 0.208 $\mu$M (Fig 8B). At $P=1\,\mu$M under those conditions, the reduced GS relaxation rate increases 2.08-fold (Fig 8C). These are prospective model predictions, not measured changes. They specify how joint binding and initial-rate measurements can distinguish affinity-mediated synergy from catalytic synergy.

The model also provides a concentration test of the PII-independent route. For the illustrative recognition law $P=T(1-\theta)^3$ and $U=T\theta$, at fixed $G$ and free-state fraction $\theta$, the drive takes the form $v(T)=a(G)+b(G)/T$. Its inverse-concentration term comes entirely from glutamine-activated AT without unmodified PII. For independent GS sites, $v=Y/(1-Y)$, so GS modification odds are affine in $1/T$ (Fig 8D). This prediction requires shared kinetic parameters and controlled free inputs; it is not inferred from the three historical assays that also varied GlnD concentration. S21 additionally shows that the direct glutamine contribution to the local log-drive sensitivity lies between zero and one. That is a bound on a partial elasticity, not on a fitted cascade Hill coefficient.

'''

CAP7=r'''\textbf{GlnD mechanisms with identical PII states predict different enzyme occupancies and kinetics.} (A) Three finite binding-and-catalysis constructions preserve the same fitted mean and every binomial PII state fraction. Points are the digitized coordinates. (B) Singly liganded enzyme fractions differ, reaching 0.5, 0.333, and 0.2 at the common binding scale. (C) PII relaxation-rate ratios in the rapid-binding, unsaturated conversion reduction; both capacities share the same multiplier. (D) An illustrative glutamine decrease from 0.5 to 0.05 mM gives different PII time courses with the same endpoint. Time is in reference-model units, not experimentally estimated minutes. The baseline is fitted once; the other models are obtained by the exact transformation in S20.'''
CAP8=r'''\textbf{Binding and catalytic regulation in the published six-state GlnE topology.} (A) Coupled GS steady-state curves are identical after a shared binding/catalysis transformation or a change in the PII--PII-UMP interaction parameter. One reference gain sets $Y=0.5$ at 0.53 mM glutamine and $T=5\,\mu$M; it is unchanged in every variant. (B) PII-binding isotherms at fixed $G=10$ mM and $U=0.3\,\mu$M distinguish binding from catalytic synergy. (C) Reduced GS relaxation rates at separately controlled $P=1\,\mu$M and $U=0.3\,\mu$M differ along a glutamine titration. These are not the coupled inputs in A. (D) The PII-independent route gives a drive affine in inverse free-target total at fixed glutamine and modification fraction. Its analytic large-pool limit is subtracted to display the concentration-dependent term. Lines are predictions, not additional observations. S21 separates source-reported fitted constants from chosen parameters and states the rapid-binding and independent-GS assumptions.'''

METHOD=r'''
\subsection*{Regulatory-state constructions and kinetic predictions}

The GlnD reference fits $q_0$, $q_2$, and $K$ using the same twelve unweighted coordinates and three initializations. Its coefficient transformation generates the $\lambda=2,4$ alternatives without additional fitting. A 701-point input grid evaluates occupancy and catalytic-rate predictions; a separate design file evaluates $K/3,K,3K$. The GlnE topology is transcribed from Fig 12 of \cite{r71}. We retain the source-reported fit example $K_G=15.6$ mM, $K_P=0.18\,\mu$M and $\alpha_1=0.17$, use $K_U=0.35\,\mu$M as a selected apparent-scale scenario, and explicitly choose $\beta=1$, $\beta'=10$, and $\alpha_2=4$. Variants 2.17 and 8.85 for $\alpha_2$ are motivated by source fits under different conditions, not by a confidence interval. A single reference gain sets $Y=0.5$ at the middle-condition midpoint; none of the transformed models is recalibrated. We do not treat these illustrations as joint fits to the three cascade coefficients.

Validation includes 100 independent stationary linear solves for PII ladders, 250 six-state receptor solves, finite-difference checks of the local sensitivity decomposition, and coupled-curve comparisons across three PII totals and twelve mechanism combinations. The step-response calculation uses buffered glutamine, a free target pool, fast regulatory binding, and unsaturated independent-site conversion. The reference PII and GS relaxation rates at the final input set shared time scales. Source constants, chosen values, transformations, initial conditions, and complete prediction grids accompany S1 Data.

'''

ABSTRACT=r'''Sharp biochemical responses can arise from different mechanisms, limiting what aggregate titrations reveal about regulation. We analyze distributive modification cycles controlled by an effector-switched bifunctional enzyme. At steady state, the shared enzyme pool cancels from adjacent free-target ratios, and local sensitivity separates into modification-ladder and effector-switch contributions. Identical independent sites contribute a ladder factor of one regardless of site number. We construct competing explanations of published bacterial nitrogen-signaling measurements under common observation rules. Finite GlnD binding-and-catalysis models preserve a fitted PII titration and every modification-state fraction while predicting different enzyme occupancies and relaxation rates. Applying the published six-state GlnE regulatory topology reveals a second ambiguity: binding synergy and catalytic synergy can compensate exactly, while another regulatory interaction changes both opposing rates without changing their ratio. The combined mechanisms therefore preserve steady cascade responses but yield distinct kinetic predictions. The same GlnE model predicts an inverse-concentration contribution from PII-independent glutamine regulation under controlled free inputs. Population alternatives and broader cascade comparisons further delimit what modification-state measurements can identify. These results distinguish structural identifiability from physiological parameter validation and show how enzyme occupancy, initial-rate, and state measurements complement mean titrations in resolving biochemical regulation.'''

AUTHOR=r'''Cells regulate proteins by adding and removing chemical tags. A sharp change in the average number of tags is often used to infer how the regulating enzymes work. We show how different enzyme mechanisms can produce the same average response in bacterial nitrogen signaling. Some alternatives differ in the distribution of tags across proteins; others preserve even that distribution. In the latter case, the enzymes themselves occupy different ligand-bound states and drive different response speeds. We construct finite GlnD models and analyze an experimentally motivated GlnE regulatory model to make these ambiguities explicit. The calculations suggest complementary measurements: counting protein modification states, measuring enzyme binding, and following initial rates or relaxation. The contribution is a way to determine which observations distinguish a proposed mechanism, while keeping experimentally established regulation separate from illustrative parameter choices.'''

for name in ['manuscript_review_with_figures.tex','manuscript_submission.tex']:
    t=(OLD/name).read_text()
    def fig(n,cap):
        s='\\begin{figure}[!ht]\n\\centering\n'
        if 'review' in name:s+=f'\\includegraphics[width=\\textwidth]{{figures/Fig{n}.pdf}}\n'
        return s+'\\begin{singlespace}\n\\caption{'+cap+'}\n\\end{singlespace}\n'+f'\\label{{fig{n}}}\n'+'\\end{figure}\n\n'
    start=t.index(r'\subsection*{A finite biochemical realization')
    end=t.index(r'\subsection*{State measurements target',start)
    t=t[:start]+GLND+fig(7,CAP7)+GLNE+fig(8,CAP8)+t[end:]
    t=re.sub(r'(\\section\*\{Abstract\}\s*).*?(\s*\\section\*\{Author summary\})',lambda m:m.group(1)+ABSTRACT+'\n'+m.group(2),t,flags=re.S)
    t=re.sub(r'(\\section\*\{Author summary\}\s*).*?(\s*\\clearpage)',lambda m:m.group(1)+AUTHOR+'\n'+m.group(2),t,flags=re.S)
    t=t.replace('The separate finite construction in Fig 7 establishes realizability under a specified recognition law, but its failed strict kinetic transfer prevents interpreting it as a validated GlnD mechanism.', 'The finite regulatory-state constructions in Figs 7 and 8 give explicit kinetic and occupancy predictions. Their coefficients remain conditional model parameters, with the earlier population construction and cross-assay constraints retained in S18--S19 and S10 Fig.')
    marker='The exact equivalence in Eq. (13) explains why adding mean-response measurements alone can leave a specific mechanistic ambiguity unresolved.'
    t=t.replace(marker,marker+' The GlnD family in Fig 7 preserves every PII state fraction, so those alternatives additionally require enzyme-occupancy or kinetic measurements. State resolution is informative for specific ambiguities, not a universal identification guarantee.')
    marker='This equivalence makes the choice of observation, not only the precision of a Hill fit, central to mechanistic inference.'
    t=t.replace(marker,marker+' Finite regulatory-state transformations extend the argument beyond population composition: binding and catalytic changes can preserve target states while altering enzyme occupancy and response speed.')
    t=t.replace('In a model of bacterial nitrogen signaling, a small population', 'In a model of bacterial nitrogen signaling, a small population')
    idx=t.index(r'\section*{Data and code availability}')
    t=t[:idx]+METHOD+t[idx:]
    t=t.replace('The new matched-model analysis, verification outputs, and figure-generation code', 'The new matched-model and regulatory-state analyses, verification outputs, and figure-generation code')
    # Retain the preceding methods/audit for reproducibility, but identify the SI figure.
    t=t.replace('The affinity-rescaling identity is verified', 'For S10 Fig, the affinity-rescaling identity is verified')
    idx=t.rfind(r'\end{document}')
    t=t[:idx]+r'''\paragraph*{S10 Fig. Finite population realizations and biochemical constraint scenarios.}
The preceding finite homogeneous/refractory comparison and affinity-rescaling calculation, retained with S18--S19. Main Figs 7 and 8 develop catalytic-state and receptor-state alternatives that do not require a refractory population.

'''+t[idx:]
    (NEW/name).write_text(t)

SUPP=r'''
\section*{S20 GlnD occupancy--catalysis transformations}

Use $E_0,E_1,E_2$ with the finite enzyme--target-complex realization in S18. Let $D_1=1+2z+z^2$, $z=G/K$, and $A=a_0+2a_1z+a_2z^2$, $B=b_0+2b_1z+b_2z^2$. The reference sets $a_1=a_2$ and $b_1=b_0$. For $\lambda>0$, change the middle binding weight to $2z/\lambda$ and middle catalytic coefficients to $\lambda a_1,\lambda b_1$. Equivalently, $K_1=\lambda K/2$ and $K_2=2K/\lambda$. Then $A$ and $B$ are unchanged, while $D_\lambda=1+2z/\lambda+z^2$. Both reduced capacities acquire the same positive multiplier
\[
\frac{U_\lambda}{U_1}=\frac{R_\lambda}{R_1}
=\frac{D_1}{D_\lambda}.
\]
The free-target ratio, mean, and all binomial modification-state fractions are unchanged. The finite-complex stationary identity concerns free targets. The displayed kinetic multiplier additionally uses a fast-binding, unsaturated-conversion reduction in which the regulated free-enzyme ensemble provides the common rate scale. It should not be applied unchanged to a strongly sequestered closed-total assay.

Ordered state activities ensure monotone capacities. For the fitted reference, the interval $1\le\lambda\le\min(b_2/b_0,a_0/a_1)$ preserves $a_0\ge\lambda a_1\ge a_2$ and $b_0\le\lambda b_1\le b_2$. Here the upper bound is the transported UR-fold scenario 6.752484, so $\lambda=1,2,4$ are all admissible. These are constructed alternatives, not three separately optimized explanations. The reference is refitted to the same twelve coordinates using three starts, with $q_0\in[0.9,0.999999]$, $q_2\in[10^{-6},0.499999]$, and $K\in[0.001,10]$ mM. The reference catalytic scale is arbitrary; only the ratios enter the mean fit. Complete coefficients are reported in S1 Data.

At $G=K$, the singly liganded enzyme fraction is $1/(1+\lambda)$ and the capacity multiplier is $2\lambda/(1+\lambda)$. Thus the $\lambda=1$ and 4 models have fractions 0.5 and 0.2 and a speed ratio of 1.6. Both have mean ligand occupancy exactly one at that input. At $G=K/3$, their mean ligand occupancies are 0.5 and approximately 0.304; at $G=3K$, they are 1.5 and approximately 1.696. These flanking inputs are useful for an average-binding assay, whereas a central single-state or rate measurement addresses the same ambiguity differently. The fitted $K$ is a conditional model scale, not a directly measured glutamine dissociation constant.

For the step-response illustration, $\dot\theta=s_D[U(G)(1-\theta)-R(G)\theta]$ and $\dot Y=s_E[v_{\rm AT}(1-Y)-v_{\rm AR}Y]$. The initial state is the steady state at $G=0.5$ mM, and the final free glutamine is 0.05 mM. $T=5\,\mu$M, and the downstream recognition law is $P=T(1-\theta)^3$, $U=T\theta$. Shared scales $s_D,s_E$ make the reference relaxation rates at the final steady state unity; transformed models retain those same scales. This is a prediction in nondimensional reference units, not a fit to experimental time courses. The source functions and generated file include both the PII and GS trajectories.

\section*{S21 Identifiable combinations in the six-state GlnE model}

The state activities and topology follow Fig 12 of Jiang, Mayo and Ninfa (2007). The ordered state weights are $(1,g,p,u,pg/\alpha_1,pu/\alpha_2)$ for $E$, $E$--Gln, $E$--PII, $E$--PII-UMP, $E$--Gln--PII, and $E$--PII--PII-UMP. AT activities are $k_P(0,\beta,1,0,\beta',0)$, and AR activities are $k_R(0,0,0,1,0,0)$. Glutamine and PII-UMP do not occupy the enzyme simultaneously in this topology. The recognition terms $P,U$ can be externally controlled free inputs. In the coupled illustrations, they are assigned the specified unmodified-trimer and count-weighted reverse inputs. This recognition assignment is distinct from the source topology itself.

The common polynomial cancels from $v=v_{\rm AT}/v_{\rm AR}$. Accordingly, $\alpha_2$ is absent from this steady drive and $\alpha_1,\beta'$ enter only through $\beta'/\alpha_1$. For $c>0$, the transformation $(\alpha_1,\beta')\mapsto(c\alpha_1,c\beta')$ multiplies both catalytic rates by $Z_E/Z'_E$ but leaves their ratio unchanged. Changing $\alpha_2$ has the same common-rate effect with a different occupancy dependence. These identities preserve any downstream steady-state construction depending only on that drive, not only the independent-site GS example. Kinetic predictions require the conversion assumptions specified in the main Methods. Summing the three states bound to unmodified PII gives a binding isotherm with halfpoint
\[
K_{P,\mathrm{eff}}=K_P\frac{1+g+u}{1+g/\alpha_1+u/\alpha_2}.
\]
This prediction depends on $\alpha_1$ separately from $\beta'$, providing a complementary measurement. It is a predicted binding halfpoint, not an assumption that every published apparent activation constant equals a microscopic dissociation constant.

A six-state continuous-time generator independently verifies the binding probabilities. Choose edges $E\leftrightarrow E$--Gln, $E\leftrightarrow E$--PII, $E\leftrightarrow E$--PII-UMP, $E$--Gln$\leftrightarrow E$--Gln--PII, $E$--PII$\leftrightarrow E$--Gln--PII, $E$--PII$\leftrightarrow E$--PII--PII-UMP, and $E$--PII-UMP$\leftrightarrow E$--PII--PII-UMP. Unit reverse rates and forward rates $(g,p,u,p/\alpha_1,g/\alpha_1,u/\alpha_2,p/\alpha_2)$ satisfy detailed balance with the stated weights. Stationary probabilities from linear solves agree with normalized weights in 250 random positive cases. These generator rates verify the equilibrium topology; they are not inferred binding-time constants.

Parameter provenance is kept separate from calculation. The source's A1/A2 fitting example reports $K_G=15.6$ mM, $K_P=0.18\,\mu$M and $\alpha_1=0.17$. The selected $K_U=0.35\,\mu$M is an apparent activation-scale scenario from Table 3. Values $\beta=1$, $\beta'=10$, and $\alpha_2=4$ are chosen illustrations. The source reports $\alpha_2$ estimates 2.17--8.85 from other global fits; these motivate variants but are not a shared-condition confidence interval. The reference $k_P/k_R$ is set once so that the independent-GS mean equals 0.5 at $G=0.53$ mM and $T=5\,\mu$M. No other midpoint or response coefficient is fitted. The transformations retain this gain. Thus the steady-curve identities are exact mathematical consequences, while the numerical occupancy and time-scale differences are conditional predictions.

For fixed $G,\theta$ and the coupled recognition law, write
\[
v(T)=\frac{k_PK_U}{k_RK_P\theta}(1-\theta)^3
\left(1+\frac{\beta'G}{\alpha_1K_G}\right)
+\frac{k_PK_U\beta G}{k_RK_G\theta}\frac1T.
\]
Only the PII-independent AT state supplies the inverse-$T$ term. The remaining pathway terms are concentration invariant at fixed free-state fractions. With independent GS sites, $v=Y/(1-Y)$; for a different downstream observation law the affine prediction applies to the effective drive rather than automatically to measured GS odds. In Fig 8D, the analytic $T\to\infty$ intercept is subtracted only to display the concentration-dependent contribution. The three historical cascade assays do not by themselves test this relation because free inputs were not held fixed and GlnD concentrations also changed.

Let the normalized positive contributions to the AT numerator be $w_P,w_G,w_{PG}$ from $p,\beta g,\beta'pg/\alpha_1$. At fixed $P,U$,
\[
\left.\frac{\partial\log v}{\partial\log G}\right|_{P,U}
=w_G+w_{PG}\in[0,1].
\]
Along the coupled input path, define $a=-d\log\theta/d\log G$ and $b=3\theta a/(1-\theta)$. Then
\[
\frac{d\log v}{d\log G}=a+(w_P+w_{PG})b+(w_G+w_{PG}).
\]
The first two terms are mediated by PII and the last is the direct glutamine contribution. A centered logarithmic finite difference checks this identity over 701 inputs. This is a local drive decomposition; endpoint-normalized response-range coefficients and nonindependent GS ladder factors require separate calculations.

The GlnD transformation and the two GlnE transformations combine without changing any steady target output under these assumptions. Thirty-six combinations of target total and mechanism are verified numerically. They define complementary experiments: resolve GlnD occupancy or PII relaxation to test the upstream family; measure GlnE binding and AT/AR initial rates across controlled ligand inputs to separate affinity and catalytic synergy; and vary free target total while monitoring its state fractions to test the PII-independent route. Parameter constraints and assay conservation must accompany a physiological interpretation.

'''

t=(OLD/'S1_Appendix.tex').read_text();idx=t.index(r'\begin{figure}[p]')
t=t[:idx]+SUPP+'\n\\clearpage\n'+t[idx:]
oldfig=r'''\begin{figure}[p]\centering\includegraphics[width=\textwidth,height=0.62\textheight,keepaspectratio]{figures/FigS10_finite_population.pdf}\begin{singlespace}\caption{\textbf{Finite population alternatives and affinity-rescaling scenarios.} The prior finite homogeneous/refractory comparison is retained with S18--S19. (A) Mean fits. (B) Distinct PII state predictions. (C) Reduced catalytic capacities. (D) Fixed-affinity saturation can distinguish the populations, whereas rescaling an unknown forward affinity restores the stipulated occupancy equivalence. The binding and kinetic assumptions differ from the six-state GlnE development in main Fig 8.}\end{singlespace}\end{figure}\clearpage
'''
t=t.replace(r'\end{document}',oldfig+'\n'+r'\end{document}');(NEW/'S1_Appendix.tex').write_text(t)

cl=(OLD/'cover_letter.tex').read_text()
start=cl.index('The main advance is');end=cl.index('Joint cascade comparisons',start)
cl=cl[:start]+'''The main advance is a constructive analysis of observational ambiguity with finite enzyme states. GlnD binding and catalytic changes can preserve the full PII modification-state distribution while changing enzyme occupancy and response speed. Analysis of the published six-state GlnE topology identifies a corresponding compensation between binding synergy and catalytic synergy. The two constructions preserve steady cascade curves and specify complementary binding, initial-rate, and concentration measurements. These are mechanistic identifiability results under explicit kinetic assumptions, not claims that the illustrative parameters establish a unique physiological mechanism.\n\n'''+cl[end:]
(NEW/'cover_letter.tex').write_text(cl)
print(NEW)
