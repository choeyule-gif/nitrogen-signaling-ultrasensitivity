"""Assumption boundaries and nuisance-resistant prospective designs.
All uncertainty/noise ranges here are declared scenarios, not experimental CIs.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, least_squares
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import regulatory_core as rc
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results'; FIG=ROOT/'figures'
for d in [OUT,FIG]:d.mkdir(exist_ok=True)
rng=np.random.default_rng(20260914)
x,y=np.loadtxt(ROOT/'data/jiangninfa2011_fig2B_digitised.csv',delimiter=',',skiprows=1).T
v,sse=rc.fit(x,y); a,b,K=rc.base_coefficients(v)
scale=1/sum(rc.glnd(K,v)[:2]);a*=scale;b*=scale

def coeff(lam):return a*np.array([1,lam,1]),b*np.array([1,lam,1])
def gen(G,lam,gamma):
 Q=np.zeros((3,3))
 for i,on in [(0,G/(lam*K/2)),(1,G/(2*K/lam))]:Q[i+1,i]=gamma*on;Q[i,i+1]=gamma
 Q[np.diag_indices(3)]=-Q.sum(axis=0)
 return Q

def explicit_ss(G,lam,T,E,hp,hm):
 aa,bb=coeff(lam);occ=rc.glnd(G,v,lam)[2];q=float(occ@aa/(occ@(aa+bb)))
 H0,H1=float(occ@hp),float(occ@hm);H=(1-q)*H0+q*H1
 # Stable positive root of H*x^2+(1+H*(E-T))*x-T=0.
 c=1+H*(E-T);free=2*T/(c+np.sqrt(c*c+4*H*T))
 ef=E/(1+H*free);ee=ef*occ;cp=hp*ee*free*(1-q);cm=hm*ee*free*q
 state=np.r_[ee,cp,cm,free*(1-q),free*q]
 return state,q,float((state[10]+cm.sum())/T)

def explicit_rhs(G,lam,T,E,hp,hm,gamma):
 aa,bb=coeff(lam);kp,km=aa/hp,bb/hm
 # koff=kcat; kon=2*h*kcat. All positive. a,b are effective specificities.
 Q=gen(G,lam,gamma)
 def f(t,z):
  ee,cp,cm=z[:3],z[3:6],z[6:9];x0,x1=z[9:]
  fp=2*aa*ee*x0-2*kp*cp;fm=2*bb*ee*x1-2*km*cm
  return np.r_[Q@ee-fp-fm,fp,fm,-(2*aa*ee*x0-kp*cp).sum()+(km*cm).sum(),-(2*bb*ee*x1-km*cm).sum()+(kp*cp).sum()]
 return f

def save(name,header,z):np.savetxt(OUT/name,z,delimiter=',',header=header,comments='')
# Full reaction balance, conservation, symmetric-saturation identity and bound.
res=conserve=bounderr=symerr=0.; rows=[]
for j in range(300):
 G=K*np.exp(rng.uniform(-4,4));T=10**rng.uniform(-1,1);eps=10**rng.uniform(-4,0);E=T*eps
 hp=10**rng.uniform(-2,2,3);hm=10**rng.uniform(-2,2,3);gamma=10**rng.uniform(-3,3)
 ys=[]
 for lam in [1,4]:
  z,q,yt=explicit_ss(G,lam,T,E,hp,hm);f=explicit_rhs(G,lam,T,E,hp,hm,gamma);dz=f(0,z)
  res=max(res,np.max(np.abs(dz))/(1+np.max(z)))
  conserve=max(conserve,abs(z[:9].sum()-E),abs(z[3:].sum()-T),abs(dz[:9].sum()),abs(dz[3:].sum()))
  _,qq,yy=explicit_ss(G,lam,T,E,hp,hp);symerr=max(symerr,abs(qq-yy));ys.append(yt)
 bounderr=max(bounderr,abs(ys[0]-ys[1])-eps)
assert res<1e-8 and conserve<1e-8 and symerr<1e-12 and bounderr<1e-12
# Independent convergence from an all-free enzyme/target initial condition.
conv=[]
for lam in [1,4]:
 hp=np.array([.2,2,10]);hm=np.array([5,.3,.1]);T=1.;E=.1;G=K;gamma=.3
 z,q,yt=explicit_ss(G,lam,T,E,hp,hm);z0=np.zeros(11);z0[0]=E;z0[9]=T
 sol=solve_ivp(explicit_rhs(G,lam,T,E,hp,hm,gamma),[0,20000],z0,method='BDF',rtol=1e-9,atol=1e-11)
 assert sol.success
 err=float(np.max(abs(sol.y[:,-1]-z)));assert err<1e-6;conv.append(err)
# Allow ligand exchange while substrate remains bound, with binding-cycle
# ratios consistent with hp/hm. Catalysis can now sustain regulatory currents.
exchange=[];exchange_res=0.;exchange_ode=[]
etas=np.r_[0.,np.geomspace(1e-4,100,25)]
for Gfac in [.3,1,3,10,30]:
 for eta in etas:
  ys=[]
  for lam in [1,4]:
   hp=np.array([.2,2,10]);hm=np.array([5,.3,.1]);T=1.;E=.1;G=K*Gfac
   z0,q0,yt=explicit_ss(G,lam,T,E,hp,hm)
   aa,bb=coeff(lam);kp,km=aa/hp,bb/hm
   rhop=np.array([.1,1,10]);rhom=np.array([10,1,.1])
   Qe=gen(G,lam,1.)
   def base(t,z):
    ee,cp,cm=z[:3],z[3:6],z[6:9];x0,x1=z[9:]
    fp=(1+rhop)*(aa*ee*x0-kp*cp);fm=(1+rhom)*(bb*ee*x1-km*cm)
    return np.r_[Qe@ee-fp-fm,fp,fm,-((1+rhop)*aa*ee*x0-rhop*kp*cp).sum()+(km*cm).sum(),-((1+rhom)*bb*ee*x1-rhom*km*cm).sum()+(kp*cp).sum()]
   def complexQ(h):
    Q=np.zeros((3,3))
    for j in [0,1]:Q[j+1,j]=eta*Qe[j+1,j]*h[j+1]/h[j];Q[j,j+1]=eta
    Q[np.diag_indices(3)]=-Q.sum(axis=0)
    return Q
   Qp,Qm=complexQ(hp*(1+rhop)/rhop),complexQ(hm*(1+rhom)/rhom)
   def fex(z):
    dz=base(0,z);dz[3:6]+=Qp@z[3:6];dz[6:9]+=Qm@z[6:9];return dz
   def residual(logz):
    z=np.exp(logz);dz=fex(z);dz[2]=z[:9].sum()-E;dz[10]=z[3:].sum()-T;return dz
   fit=least_squares(residual,np.log(z0),xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=1200)
   zz=np.exp(fit.x);err=float(np.max(abs(fex(zz))));exchange_res=max(exchange_res,err)
   assert fit.success and err<1e-7 and abs(zz[:9].sum()-E)<1e-7 and abs(zz[3:].sum()-T)<1e-7
   if Gfac==10 and eta==1:
    check=solve_ivp(lambda t,z:fex(z),[0,20000],z0,method='BDF',rtol=1e-9,atol=1e-11)
    assert check.success;ce=float(max(abs(check.y[:,-1]-zz)));assert ce<1e-6;exchange_ode.append(ce)
   qf=zz[10]/(zz[9]+zz[10]);qt=(zz[10]+zz[6:9].sum())/T
   ys.append((qf,qt));exchange.append([Gfac,eta,lam,qf,qt,q0,err])
save('substrate_bound_ligand_exchange.csv','G_over_K,bound_exchange_rate,lambda,free_modified_fraction,total_modified_fraction,no_exchange_free_fraction,residual',exchange)
# A disclosed asymmetric-affinity scenario, plus worst case envelope from random affinities.
epsgrid=np.geomspace(1e-4,1,41);ggrid=K*np.geomspace(.02,50,81);gap=np.zeros((41,81));h0=np.array([.2,2,10]);h1=np.array([5,.3,.1])
for i,eps in enumerate(epsgrid):
 for j,G in enumerate(ggrid):
  zs=[explicit_ss(G,l,1,eps,h0,h1) for l in [1,4]];gap[i,j]=abs(zs[0][2]-zs[1][2]);rows.append([eps,G/K,*[z[2] for z in zs],gap[i,j],eps])
save('total_observation_boundary.csv','enzyme_target_ratio,G_over_K,total_fraction_lambda1,total_fraction_lambda4,gap,rigorous_pairwise_bound',rows)
# Finite regulator binding: exact steady output, non-universal transient.
g0,g1=K*10,K;times=np.linspace(0,12,241);dyn=[];steadyerrs=[];lagerrs=[]
for gamma in [.01,.1,1,10,100]:
 sols=[]
 for lam in [1,4]:
  aa,bb=coeff(lam);e0=rc.glnd(g0,v,lam)[2];q0=float(rc.theta(g0,v,lam));Q=gen(g1,lam,gamma)
  def fun(t,z):return np.r_[Q@z[:3],(z[:3]@aa)*(1-z[3])-(z[:3]@bb)*z[3]]
  sol=solve_ivp(fun,[0,times[-1]],np.r_[e0,q0],t_eval=times,method='BDF',rtol=1e-9,atol=1e-11)
  assert sol.success;sols.append(sol.y[3]);dyn.extend(np.column_stack([np.full(len(times),gamma),np.full(len(times),lam),times,sol.y[3]]))
  ee=rc.glnd(g1,v,lam)[2];qq=float(rc.theta(g1,v,lam));steadyerrs.append(np.max(abs(fun(0,np.r_[ee,qq]))))
  rate=float(ee@(aa+bb));fast=qq+(q0-qq)*np.exp(-rate*times)
  lagerrs.append([gamma,lam,float(np.max(abs(sol.y[3]-fast)))])
save('finite_binding_dynamics.csv','binding_rate,lambda,time_reference_units,PII_fraction',dyn)
save('rapid_binding_error.csv','binding_rate,lambda,max_transient_fraction_error',lagerrs)
# Nuisance-resistant design. Rate scales are unknown but SHARED across doses.
def R(z,lam):return (1+2*z+z*z)/(1+2*z/lam+z*z)
# Search candidate pair with K uncertain by factor 2. log contrasts only.
doses=np.geomspace(.02,50,101);unc=np.geomspace(.5,2,31)
best=(-1,None,None,None)
for i,lo in enumerate(doses):
 for hi in doses[i+1:]:
  delta=np.log(R(lo/unc,4)/R(hi/unc,4));score=float(np.min(abs(delta)))
  if score>best[0]:best=(score,lo,hi,delta)
Dmin,lo,hi,_=best
# Glutamine-free vs glutamine-containing GlnE binding removes KP; U=0 removes alpha2.
# halfpoint Kapp=KP*(1+G/KG)/(1+G/(alpha1*KG)).
def kapp(G,KG,KP,alpha):return KP*(1+G/KG)/(1+G/(alpha*KG))
kggrid=15.6*np.geomspace(.5,2,31);a1grid=.17*np.geomspace(.5,2,31)
def econ(G,kg,al,c=4):return np.log((kapp(G,kg,1,c*al)/kapp(0,kg,1,c*al))/(kapp(G,kg,1,al)/kapp(0,kg,1,al)))
e_delta=np.array([econ(10,kg,al) for kg in kggrid for al in a1grid]);emin=float(e_delta.min())
# Composite nuisance envelopes must vary independently under each hypothesis.
# A known-scale template difference alone is NOT a minimax discrimination guarantee.
occ1=lambda z,l:(2*z/l)/(1+2*z/l+z*z)
omin=float(min(occ1(1/unc,1)));omax=float(max(occ1(1/unc,4)));ogap=omin-omax
assert ogap>0
# GlnE: explicit disjoint effect-size classes; these ranges are scenarios, not CIs.
eal=.17*np.geomspace(.8,1.25,21)
eenv=[]
for c in [1,4]:
 vals=np.array([np.log(kapp(10,kg,1,c*al)) for kg in kggrid for al in eal]);eenv.append([float(vals.min()),float(vals.max())])
egap=eenv[1][0]-eenv[0][1];assert egap>0
# For broader alpha1 factor-two classes, the classes touch and no guarantee exists.
# Three binding titrations at G=0,1,10 identify alpha1 and KG ideally, independent of KP.
inversion=0.
for kg in kggrid:
 for al in [.136,.17,.2125,.544,.68,.85]:
  Gs=np.array([1.,10.]);rr=kapp(Gs,kg,1,al);slope,intercept=np.polyfit(Gs,Gs/(1-rr),1)
  ahat=1-1/slope;ghat=intercept*(1-ahat)/ahat
  inversion=max(inversion,abs(ahat-al),abs(ghat-kg))
assert inversion<1e-9
zcrit=norm.ppf(.95);zpower=norm.ppf(.8);design=[]
# variance_multiplier=1 for one occupancy fraction; =2 for paired log estimates.
for assay,delta,mult in [('GlnD_occupancy',ogap,1),('GlnE_paired_halfpoint',egap,2),('GlnD_capacity_conditional',Dmin,2)]:
 for sigma in [.05,.1,.2,.3]:
  n=int(np.ceil(mult*sigma*sigma*(zcrit+zpower)**2/(delta*delta)))
  power=float(norm.cdf(np.sqrt(n)*delta/(np.sqrt(mult)*sigma)-zcrit))
  design.append([assay,sigma,n,delta,power,mult])
with (OUT/'noise_design.csv').open('w') as f:
 f.write('assay,measurement_SD,n_per_condition,min_contrast,power,variance_multiplier\n')
 for row in design:f.write(','.join(map(str,row))+'\n')
mc=[]
for assay,sigma,n,delta,powr,mult in design:
 if sigma!=.2:continue
 sd=np.sqrt(mult/n)*sigma;null=rng.normal(0,sd,100000);alt=rng.normal(delta,sd,100000)
 fpr=float(np.mean(null/sd>zcrit));pwr=float(np.mean(alt/sd>zcrit));assert abs(fpr-.05)<.005 and abs(pwr-powr)<.006
 mc.append({'assay':assay,'false_positive_rate':fpr,'power':pwr,'analytic_power':powr})
# Continuous alternatives: near-equivalent models require unbounded precision.
ld=np.linspace(1.02,4,151);detect=[]
for lam in ld:
 delta=max(0.,omin-float(max(occ1(1/unc,lam))))
 n=np.inf if delta==0 else int(np.ceil(.1**2*(zcrit+zpower)**2/delta**2));detect.append([lam,delta,n])
save('continuous_alternative_resolution.csv','lambda,worst_log_contrast,n_per_condition_fractionSD0p1',detect)
# Ligand depletion correction for isolated GlnE binding (U=0, G buffered).
# Bound B solves Ptotal=Pfree+Etotal*Pfree/(Keff+Pfree).
def bound(Ptot,Et,keff):
 c=Et+keff-Ptot;pf=2*Ptot*keff/(c+np.sqrt(c*c+4*Ptot*keff))
 return pf/(keff+pf)
dep=[]
for E in [.001,.01,.1,1]:
 for G in [0,10]:
  for c in [1,4]:
   ke=kapp(G,15.6,.18,.17*c);nom=ke+E/2
   assert abs(bound(nom,E,ke)-.5)<1e-12
   dep.append([E,G,c,ke,nom,nom-E/2])
save('GlnE_depletion_correction.csv','GlnE_total_uM,G_mM,synergy_multiplier,free_halfpoint_uM,nominal_halfpoint_uM,corrected_halfpoint_uM',dep)
# Single-input nuisance re-fit identities, including symmetric-dose blindness.
single=float(abs(R(1,4)*(1/R(1,4))-1));sym=float(abs(np.log(R(.1,4)/R(10,4))))
# Capacity rate contrasts cease to be guaranteed with unknown dose-specific scales.
report={'scope':'Mechanistic boundary and prospective design, no new experimental observations or physiological confidence intervals.',
 'GlnD':{'K_mM':K,'full_mass_action_random_cases':300,'stationary_residual':float(res),'conservation_error':float(conserve),'symmetric_affinity_total_identity_error':float(symerr),'pairwise_total_gap_bound':'min(1,Etotal/Ttotal), for one target unit per complex, same totals and free fraction','bound_excess':float(bounderr),'independent_ODE_convergence_errors':conv,'bound_ligand_exchange_steady_residual':exchange_res,'bound_exchange_independent_ODE_errors':exchange_ode,'bound_ligand_exchange_scope':'Ligand exchange in substrate-bound complexes; thermodynamic binding ratios h*(1+rho)/rho; koff/kcat rho+=(.1,1,10), rho-=(10,1,.1); buffered G; G/K=.3,1,3,10,30','finite_binding_stationary_residual':float(max(steadyerrs)),'single_condition_unknown_scale_identity_error':single,'symmetric_dose_contrast':sym},
 'design':{'criterion':'maximin absolute log contrast on a declared finite grid, shared unknown scale canceled','GlnD_doses_over_nominal_K':[float(lo),float(hi)],'GlnD_doses_mM':[float(lo*K),float(hi*K)],'GlnD_K_uncertainty_factors':[.5,2],'GlnD_min_log_contrast':Dmin,'GlnE_doses_mM':[0,10],'GlnE_U_uM':0,'GlnE_KG_uncertainty_factors':[.5,2],'GlnE_pointwise_log_contrast_not_profiled':emin,'GlnE_profiled_envelopes':eenv,'GlnE_independent_alpha1_factors':[.8,1.25],'GlnE_min_profiled_log_contrast':egap,'GlnD_occupancy_null_min':omin,'GlnD_occupancy_alternative_max':omax,'GlnD_occupancy_min_gap':ogap,'GlnE_three_titration_inversion_error':float(inversion),'noise_model':'normal error in absolute occupancy fraction, or independent normal log-capacity/log-halfpoint estimates; known SD, shared scale cancels; no dose-specific bias','test':'one-sided alpha .05 at worst composite-null boundary for occupancy and GlnE; target power .80; GlnD capacity numbers CONDITIONAL on calibrated baseline shape, not profiled K','noise_design':design,'Monte_Carlo_checks':mc},
 'limits':['Finite binding rates leave the free-enzyme-exchange buffered steady states unchanged but affect transients; ligand exchange in substrate-bound complexes can break that identity.','Monovalent target-unit saturation model does not assert complete trimer sequestration stoichiometry.','State-dependent substrate affinity is held fixed across symmetry alternatives; symmetric affinities preserve total output even at saturation.','Paired capacity protocol requires pre-equilibrated enzyme and controlled linear-response/free-target conditions.','GlnE binding halfpoint uses U=0 and buffered glutamine, with total PII depletion corrected using measured GlnE total.','Unknown dose-specific biases or unconstrained functional regulation can destroy guaranteed discrimination.','Overlapping nuisance envelopes cannot be resolved by replication alone; occupancy guarantees require a separated lambda alternative.','GlnD capacity contrast alone is conditional on independently calibrated baseline K and catalytic shape; its K envelope is not a refitted-template guarantee.','GlnE alpha1 class ranges are explicitly specified effect-size hypotheses, not physiological confidence intervals.']}
(OUT/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
colors=['#477C9B','#B77A75','#7B9778','#A59062','#887EAA']
def finish(fig,axs,name):
 for ax,L in zip(axs.flat,'ABCD'):ax.text(-.17,1.04,L,transform=ax.transAxes,fontweight='bold',fontsize=12)
 fig.savefig(FIG/(name+'.pdf'),dpi=600);fig.savefig(FIG/(name+'.png'),dpi=600);plt.close(fig)
fig,axs=plt.subplots(2,2,figsize=(7.3,5.8),layout='constrained')
ax=axs[0,0]
for lam,c in zip([1,4],colors):
 q=[];yt=[]
 for G in ggrid:
  _,qq,yy=explicit_ss(G,lam,1,.3,h0,h1);q.append(qq);yt.append(yy)
 ax.plot(ggrid/K,yt,color=c,label=f'Total, λ={lam}')
ax.plot(ggrid/K,q,'k--',lw=1,label='Free, both models');ax.set(xscale='log',xlabel='Glutamine / K',ylabel='Modified target fraction');ax.legend(fontsize=7)
ax=axs[0,1];m=ax.pcolormesh(ggrid/K,epsgrid,gap,shading='auto',cmap='viridis',rasterized=True);ax.set(xscale='log',yscale='log',xlabel='Glutamine / K',ylabel='Enzyme / target total');fig.colorbar(m,ax=ax,label='Total-fraction difference')
ax=axs[1,0]
for gamma,c in zip([.01,.1,1,10,100],colors):
 z=np.array(dyn);z=z[(z[:,0]==gamma)&(z[:,1]==4)];ax.plot(z[:,2],z[:,3],color=c,label=f'Binding rate {gamma:g}')
ax.set(xlabel='Time (reference units)',ylabel='PII modification fraction');ax.legend(fontsize=7)
ax=axs[1,1];ex=np.array(exchange)
for Gfac,col in zip([1,10,30],colors):
 x1=ex[(ex[:,0]==Gfac)&(ex[:,2]==1)];x4=ex[(ex[:,0]==Gfac)&(ex[:,2]==4)]
 ax.semilogx(x1[1:,1],abs(x1[1:,3]-x4[1:,3]),color=col,label=f'G/K={Gfac:g}')
ax.set(xlabel='Ligand-exchange rate in target complexes',ylabel='Free-fraction difference');ax.legend(fontsize=7)
finish(fig,axs,'Fig9_assumption_boundary')
fig,axs=plt.subplots(2,2,figsize=(7.3,5.8),layout='constrained')
ax=axs[0,0];gg=np.geomspace(.05,20,250)
for lam,c in zip([1,4],colors):
 curves=np.array([occ1(gg/u,lam) for u in unc]);ax.fill_between(gg,curves.min(axis=0),curves.max(axis=0),color=c,alpha=.22)
 ax.plot(gg,occ1(gg,lam),color=c,label=f'λ={lam}')
ax.axvline(1,color='0.6',ls=':',lw=.8)
ax.set(xscale='log',xlabel='Glutamine / nominal K',ylabel='Singly liganded GlnD fraction');ax.legend(fontsize=7)
ax=axs[0,1];gg=np.linspace(0,10,200)
for c,col in zip([1,4],colors):
 cv=np.array([kapp(gg,kg,1,c*al) for kg in kggrid for al in eal]);ax.fill_between(gg,cv.min(axis=0),cv.max(axis=0),color=col,alpha=.2)
 ax.plot(gg,kapp(gg,15.6,1,.17*c),color=col,label=f'α₁ class × {c}')
ax.set(xlabel='Glutamine (mM), PII-UMP = 0',ylabel='PII halfpoint / zero-Gln halfpoint');ax.legend(fontsize=7)
ax=axs[1,0]
for j,(label,col) in enumerate(zip(['GlnD occupancy (fraction SD)','GlnE halfpoints (log SD)'],colors)):
 ds=design[j*4:(j+1)*4];ax.plot([z[1] for z in ds],[z[2] for z in ds],'o-',c=col,label=label)
ax.set(xlabel='SD of assay estimate (units in legend)',ylabel='Replicates for 80% power',ylim=(0,None));ax.legend(fontsize=7)
ax=axs[1,1];dd=np.array(detect);good=np.isfinite(dd[:,2]);ax.semilogy(dd[good,0],dd[good,2],color=colors[0]);ax.axvspan(1,1.25,color='0.9');ax.text(1.06,3,'Overlapping\nenvelopes',rotation=90,fontsize=7)
ax.set(xlabel='GlnD alternative λ',ylabel='Replicates at fraction SD = 0.10',title='K uncertain by a factor of two')
finish(fig,axs,'Fig10_discriminating_design')
print(json.dumps(report,indent=2))
