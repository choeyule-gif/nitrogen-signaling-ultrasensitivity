"""Rebuild the revision figures and numerical summaries.
Input points are the author-supplied digitization of published measurements,
not an experimental raw-data file. Run from any directory.
"""
from pathlib import Path
import json, math
import numpy as np
from scipy.optimize import least_squares, brentq
from scipy.special import gammaln, logsumexp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'data'; F=ROOT/'diagnostic_figures'; F.mkdir(exist_ok=True)
xy=np.loadtxt(D/'jiangninfa2011_fig2B_digitised.csv',delimiter=',',skiprows=1)
x,y=xy.T; N=len(x)
BOUNDS=([0,1.5,-8,.1],[1.5,3,5,8])
def hill(xx,p):
 lo,hi,ls,h=p
 return lo+(hi-lo)/(1+(np.asarray(xx)/np.exp(ls))**h)
def fit(xx,yy,h=None,start=None):
 p0=np.array([.45,2.97,np.log(.56),2.]) if start is None else np.array(start)
 if h is None:
  r=least_squares(lambda p:hill(xx,p)-yy,p0,bounds=BOUNDS,xtol=1e-10,ftol=1e-10,gtol=1e-10)
  return r.x,float(r.fun@r.fun)
 r=least_squares(lambda q:hill(xx,[*q,h])-yy,p0[:3],bounds=(np.array(BOUNDS[0])[:3],np.array(BOUNDS[1])[:3]),xtol=1e-10,ftol=1e-10,gtol=1e-10)
 return np.r_[r.x,h],float(r.fun@r.fun)
p,sse=fit(x,y); fits={};rows=[]
for h in [1,2,3,None]:
 pp,ss=fit(x,y,h); k=5 if h is None else 4
 aic=N*np.log(ss/N)+2*k;aicc=aic+2*k*(k+1)/(N-k-1)
 fits[str(h)]=(pp,ss);rows.append([0 if h is None else h,ss,k,aic,aicc])
rows=np.array(rows);rows[:,3]-=rows[:,3].min();rows[:,4]-=rows[:,4].min()
np.savetxt(D/'model_comparison.csv',rows,delimiter=',',header='exponent_0_means_free,SSE,k_including_variance,delta_AIC,delta_AICc',comments='')
def hprof(h):return N*np.log(fit(x,y,h)[1]/sse)
hci=[brentq(lambda h:hprof(h)-3.841458820694,1,p[3]),brentq(lambda h:hprof(h)-3.841458820694,p[3],3)]
def swingfit(w):
 def fun(q):
  low,ls,h=q;t=low/3;hi=3*w*t/(1-t+w*t)
  return hill(x,[low,hi,ls,h])-y
 r=least_squares(fun,[p[0],p[2],p[3]],bounds=([1e-8,-8,.1],[1.5,5,8]),xtol=1e-10,ftol=1e-10,gtol=1e-10)
 return float(r.fun@r.fun)
tlo,thi=p[:2]/3
swing=(thi/(1-thi))/(tlo/(1-tlo))
wlo=np.exp(brentq(lambda z:N*np.log(swingfit(np.exp(z))/sse)-3.841458820694,np.log(10),np.log(swing)))
# An alternative observation model: a refractory fraction with upper plateau fixed at 3.
r=least_squares(lambda q:hill(x,[q[0],3,q[1],q[2]])-y,[p[0],p[2],p[3]],bounds=([0,-8,.1],[1.5,5,8]))
p_ref=np.array([r.x[0],3,r.x[1],r.x[2]])

def response(l,pp,k,reverse='linear'):
 t=hill(l,pp)/3; fw=(1-t)**3;rv=t if reverse=='linear' else 1-fw
 return k*fw/(k*fw+rv)
def limits(pp,k,reverse):return response(0,pp,k,reverse),response(1e15,pp,k,reverse)
def calibration(pp,smid,reverse='linear'):
 def fun(logk):
  k=np.exp(logk);a,b=limits(pp,k,reverse)
  return response(smid,pp,k,reverse)-(a+b)/2
 return np.exp(brentq(fun,-28,28))
def score(pp,k,reverse='linear'):
 # Analytical inversion: solve for upstream theta at each downstream range threshold.
 tl,th=pp[:2]/3;a,b=limits(pp,k,reverse);ts=[]
 for q in [.1,.9]:
  target=a+q*(b-a)
  def fun(t):
   fw=(1-t)**3;rv=t if reverse=='linear' else 1-fw
   return k*fw/(k*fw+rv)-target
  t=brentq(fun,tl+1e-13,th-1e-13)
  ls=pp[2]+np.log((th-t)/(t-tl))/pp[3]
  ts.append(ls)
 return np.log(81)/(ts[1]-ts[0])
pii=np.array([.5,5.,36.]);mid=np.array([.600,.530,.225]);observed=np.array([5.23,6.46,4.48])
ks=np.array([calibration(p,m) for m in mid]);pred=np.array([score(p,k) for k in ks]);res=observed/pred
# Coordinate perturbation sensitivity analysis, not experimental confidence intervals.
rng=np.random.default_rng(20260909); B=2000; ps=[];preds=[];sw=[];failed=0
for i in range(B):
 xx=x*np.exp(rng.normal(0,.06,N));yy=y+rng.normal(0,.06,N)
 try:
  pp,ss=fit(xx,yy,start=p)
  vals=np.array([score(pp,calibration(pp,m)) for m in mid])
  ps.append(pp);preds.append(vals)
  l,h=pp[:2]/3;sw.append((h/(1-h))/(l/(1-l)) if h<1 else np.inf)
 except (ValueError,RuntimeError,FloatingPointError):failed+=1
ps=np.array(ps);preds=np.array(preds);sw=np.array(sw)
np.savez_compressed(D/'coordinate_perturbations.npz',parameters=ps,composite_coefficients=preds,regulatory_swing=sw)
q_h=np.quantile(ps[:,3],[.025,.5,.975]);q_res=np.quantile(observed/preds,[.025,.5,.975],axis=0)
np.savetxt(D/'cascade_comparison.csv',np.c_[pii,mid,observed,ks,pred,res,q_res[0],q_res[2]],delimiter=',',header='PII_uM,midpoint_mM,reported_coefficient,kappa,calculated_range_coefficient,reported_over_calculated,perturbation_p025,perturbation_p975',comments='')
# State distributions at the same mean modification (illustrative, not fitted).
n=12;ii=np.arange(n+1);logbin=gammaln(n+1)-gammaln(ii+1)-gammaln(n-ii+1)
def dist(j):
 z=logbin+j*(ii-n/2)**2;return np.exp(z-logsumexp(z))
j=brentq(lambda j:4*np.sum(dist(j)*(ii-n/2)**2)/n-2.2,0,.3)
probs=np.c_[ii,dist(0),dist(j)];np.savetxt(D/'illustrative_state_distributions.csv',probs,delimiter=',',header='modified_sites,independent,nonbinomial_nu_2p2',comments='')
# Formal/specific checks.
for nn in [1,3,12]:
 z=np.linspace(-8,8,51);i=np.arange(nn+1)
 for zz in z:
  ww=np.exp(gammaln(nn+1)-gammaln(i+1)-gammaln(nn-i+1)+i*zz);ww/=ww.sum();mu=ww@i;var=ww@(i-mu)**2
  assert abs(nn*var/(mu*(nn-mu))-1)<1e-8
for m,k in zip(mid,ks):
 a,b=limits(p,k,'linear');assert abs(response(m,p,k)-(a+b)/2)<1e-10
# Local logit slope for the physical mean fraction.
def nlocal(l,pp):
 t=hill(l,pp)/3;q=(l/np.exp(pp[2]))**pp[3]
 deriv=(pp[1]-pp[0])/3*pp[3]*q/(1+q)**2
 return deriv/(t*(1-t))
S=np.exp(p[2]);nloc=float(nlocal(S,p))
summary=dict(parameters=dict(Umin=p[0],Umax=p[1],S_mM=S,h=p[3]),SSE=sse,profile_h_interval=hci,coordinate_h_quantiles=q_h.tolist(),coordinate_assumptions=dict(B=B,seed=20260909,log_x_sd=.06,y_sd=.06,usable=len(ps),failed=failed),plateau_swing=swing,profile_swing_lower=wlo,UT_inhibition_lower=wlo/(7.93/1.17),refractory_h=p_ref[3],refractory_LR=N*np.log(float(r.fun@r.fun)/sse),local_slope_at_S=nloc,nu2_example_j=j,cascade_predicted=pred.tolist(),cascade_ratio=res.tolist(),cascade_ratio_perturbation_q=q_res.tolist(),source='Author-supplied by-eye digitization of Jiang and Ninfa 2011 Fig. 2B; no original experimental replicates')
(D/'analysis_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
# Shared scientific figure style.
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.bbox':'tight'})
blue='#2166ac';orange='#d95f02';green='#1b7837';purple='#762a83';gray='#666666'
def save(fig,name):
 fig.savefig(F/(name+'.pdf'));fig.savefig(F/(name+'.png'),dpi=220);plt.close(fig)
def letter(ax,t):ax.text(-.19,1.20,t,transform=ax.transAxes,weight='bold',fontsize=12,va='top')
def box(ax,x,y,w,h,text,fc='white'):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.012',fc=fc,ec='#777777',lw=.8));ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=9)
def arr(ax,a,b,color=gray):ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=12,lw=1.3,color=color))
# Figure 1.
fig=plt.figure(figsize=(7.1,6.0));gs=fig.add_gridspec(2,2,height_ratios=[1,1.08],hspace=.52,wspace=.42)
a=fig.add_subplot(gs[0,:]);a.set(xlim=(0,1),ylim=(0,1));a.axis('off');letter(a,'A')
box(a,.02,.40,.14,.25,'Glutamine');box(a,.25,.40,.14,.25,'GlnD\nUT / UR','#e7f0f8');box(a,.48,.40,.14,.25,'PII\n0–3 UMP');box(a,.71,.40,.14,.25,'ATase\nAT / AR','#e7f0f8')
arr(a,(.16,.525),(.245,.525));arr(a,(.39,.525),(.475,.525));arr(a,(.62,.525),(.705,.525));arr(a,(.85,.525),(.90,.525));a.text(.96,.525,'GS\n0–12 AMP',ha='center',va='center',fontsize=9)
a.text(.32,.23,'Effector switch',color=blue,ha='center');a.text(.55,.23,'Modification ladder',color=orange,ha='center');a.text(.78,.23,'Collective readout',color=green,ha='center');a.text(.50,.87,'Two cycles connect glutamine to glutamine-synthetase modification',ha='center',fontsize=11)
b=fig.add_subplot(gs[1,0]);b.axis('off');letter(b,'B');b.text(.02,.87,r'$J_i^+=e\,a_i s_i$',fontsize=13);b.text(.02,.66,r'$J_i^-=e\,b_i f(l)s_{i+1}$',fontsize=13);b.text(.02,.38,r'$\frac{s_{i+1}}{s_i}=\frac{a_i/b_i}{f(l)}$',fontsize=16);b.text(.02,.09,r'$n_H^{\rm loc}=\nu\,\eta$',fontsize=18,color=blue)
c=fig.add_subplot(gs[1,1]);letter(c,'C');xx=np.logspace(-2,2,300)
for nn,ls,col in [(1,'-',blue),(3,'--',orange),(12,':',green)]:c.plot(xx,1/(1+xx),ls,lw=2.5,label=f'n = {nn}',color=col)
c.set(xscale='log',xlabel=r'Free input $l/\hat K$',ylabel=r'Mean modification $\theta$',ylim=(-.03,1.03));c.legend(frameon=False);c.set_title('Identical mean response')
save(fig,'Fig1_architecture')
# Figure 2.
fig,axs=plt.subplots(1,3,figsize=(7.2,2.95));fig.subplots_adjust(wspace=.56)
a,b,c=axs
for ax,lab in zip(axs,'ABC'):letter(ax,lab)
v=np.linspace(.55,3,300)
for val in [1,2,3,4,6]:a.plot(v,val/v,color='#bdbdbd',lw=.9);a.text(min(2.72,val/.65),val/min(2.72,val/.65),f'{val}',fontsize=7,color=gray)
a.axvline(1,color=blue,ls='--',lw=1);a.axvline(3,color=purple,lw=1.5);a.set(xlim=(.5,3.1),ylim=(.5,3.5),xlabel=r'Ladder factor $\nu$',ylabel=r'Switch factor $\eta$',title='Local coordinates only')
rat=np.logspace(0,4,300);mx=2*np.sqrt(rat)/(1+np.sqrt(rat));b.plot(rat,mx,color=orange);b.scatter([80/70,700/150],2*np.sqrt([80/70,700/150])/(1+np.sqrt([80/70,700/150])),color='black',s=15);b.annotate('Mg',xy=(80/70,1.035),xytext=(1.8,1.15),fontsize=8);b.annotate('Mn',xy=(700/150,1.367),xytext=(9,1.3),fontsize=8);b.axhline(2,ls=':',color=gray);b.set(xscale='log',ylim=(.98,2.04),xlabel=r'Affinity ratio $A_{UT}/A_{UR}$',ylabel=r'Maximum $\eta$',title='Reciprocal switch')
xx=np.logspace(-2,1.3,400);c.plot(xx,nlocal(xx,p),color=blue,label=r'$n_H^{\rm loc}(l)$');c.axhline(p[3],color=orange,ls='--',label=r'$n_H^{\rm fit}$');c.scatter([S],[nloc],color=blue,s=16);c.set(xscale='log',xlabel='Glutamine (mM)',ylabel='Response coefficient',title='Local and fitted slopes',ylim=(0,2.3));c.legend(frameon=False,loc='lower center');save(fig,'Fig2_sensitivity')
# Figure 3.
fig,axs=plt.subplots(2,2,figsize=(7.1,5.6));fig.subplots_adjust(hspace=.58,wspace=.40)
for ax,lab in zip(axs.flat,'ABCD'):letter(ax,lab)
a,b,c,d=axs.flat;xx=np.logspace(-2,1.25,350);curves=np.array([hill(xx,pp) for pp in ps]);band=np.quantile(curves,[.025,.975],axis=0)
a.fill_between(xx,*band,color='#d8e7f4',label='Coordinate perturbations');a.plot(xx,hill(xx,p),color=blue,label='Free-exponent fit');a.scatter(x,y,s=24,facecolors='white',edgecolors='black',zorder=5,label='Reconstructed points');a.set(xscale='log',xlabel='Glutamine (mM)',ylabel='UMP per trimer',ylim=(0,3.15));a.legend(frameon=False,fontsize=7,loc='upper right')
sigma=np.sqrt(fits['2'][1]/(N-3))
for h,col in [(1,gray),(2,blue),(3,orange)]:b.plot(x,(y-hill(x,fits[str(h)][0]))/sigma,'o-',ms=3,lw=1,color=col,label=f'h = {h}')
b.axhline(0,color='black',lw=.6);b.set(xscale='log',xlabel='Glutamine (mM)',ylabel='Residual / common scale');b.legend(frameon=False,ncol=3,fontsize=7)
c.hist(ps[:,3],bins=35,density=True,color='#d8e7f4',edgecolor='white');c.axvline(p[3],color=blue);c.axvline(q_h[0],color=gray,ls=':');c.axvline(q_h[2],color=gray,ls=':');c.set(xlabel='Fitted exponent',ylabel='Perturbation density',title='Assumed coordinate errors');c.text(.03,.93,f'95% range: {q_h[0]:.2f}–{q_h[2]:.2f}',transform=c.transAxes,va='top',fontsize=8)
ww=np.logspace(1.5,5,110);loss=np.array([N*np.log(swingfit(w)/sse) for w in ww]);d.plot(ww,loss,color=purple);d.axhline(3.841458820694,color=gray,ls='--');d.axvline(wlo,color=purple,ls=':');d.set(xscale='log',xlabel=r'Regulatory swing $R_{fwd}R_{rev}$',ylabel='Profile likelihood loss',ylim=(0,12));d.text(.05,.9,f'Nominal lower bound: {wlo:.0f}',transform=d.transAxes,fontsize=8);save(fig,'Fig3_titration')
# Figure 4.
fig,axs=plt.subplots(1,3,figsize=(7.2,3.0));fig.subplots_adjust(wspace=.57)
for ax,lab in zip(axs,'ABC'):letter(ax,lab)
a,b,c=axs;xx=np.logspace(-2,1,400);colors=[blue,orange,green]
for t,m,k,col in zip(pii,mid,ks,colors):
 aa,bb=limits(p,k,'linear');a.plot(xx,(response(xx,p,k)-aa)/(bb-aa),color=col,label=f'{t:g} μM');a.scatter([m],[.5],s=13,color=col)
a.set(xscale='log',xlabel='Glutamine (mM)',ylabel='Normalized GS response',title='Calibrated responses');a.legend(frameon=False,fontsize=7)
ind=np.arange(3);b.bar(ind-.19,pred,.35,color=blue,label='Calculated');b.bar(ind+.19,observed,.35,color='white',edgecolor='black',hatch='///',label='Reported');b.set(xticks=ind,xticklabels=['0.5','5','36'],xlabel='PII (μM)',ylabel='Response coefficient',ylim=(0,8.5));b.legend(frameon=False,fontsize=7)
c.errorbar(pii,res,yerr=np.maximum(0,np.array([res-q_res[0],q_res[2]-res])),fmt='o-',color=orange,capsize=3);c.axhline(1,color=gray,ls='--');c.set(xscale='log',xlabel='PII (μM)',ylabel='Reported / calculated',ylim=(.65,2));c.set_title('Input-curve uncertainty');save(fig,'Fig4_cascade')
# Figure 5.
fig,axs=plt.subplots(1,3,figsize=(7.2,3.25));fig.subplots_adjust(wspace=.60)
for ax,lab in zip(axs,'ABC'):letter(ax,lab)
a,b,c=axs;xx=np.logspace(-2,2,300)
a.plot(xx,xx/(1+xx),color=blue,label='One ligand');a.plot(xx,2*xx**2/(1+xx**2),color=orange,label='Two ligands');a.set(xscale='log',xlabel='Scaled glutamine',ylabel='Bound ligands / enzyme',title='Measure occupancy');a.legend(frameon=False,fontsize=7,loc='lower right');a.text(.5,-.40,'Illustrative binding laws',transform=a.transAxes,ha='center',fontsize=7,color=gray)
b.axis('off');b.set_title('Vary enzyme and target');b.set(xlim=(0,1),ylim=(0,1));b.text(.5,.91,'Target level',ha='center');b.text(.34,.74,'1×',ha='center');b.text(.74,.74,'4×',ha='center');b.text(-.08,.52,'E 1×',fontsize=8);b.text(-.08,.22,'E 2×',fontsize=8)
for bx in [.21,.61]:
 for by in [.08,.39]:box(b,bx,by,.27,.25,'Gln',fc='#eef3f8')
b.text(.50,-.15,'Titrate glutamine at each pair\nRead midpoint, range and slope',ha='center',va='top',fontsize=8)
c.bar(ii-.2,probs[:,1],.4,color=blue,label=r'Independent: $\nu=1$');c.bar(ii+.2,probs[:,2],.4,color=orange,label=r'Example: $\nu=2.2$');c.set(xlabel='Modified sites on GS',ylabel='State probability',title=r'Same mean: $\theta=0.5$',xticks=[0,3,6,9,12]);c.legend(frameon=False,fontsize=7);c.text(.5,-.40,'Illustrative distributions, not a fit',transform=c.transAxes,ha='center',fontsize=7,color=gray);save(fig,'Fig5_experiments')
print('Five figures written.',flush=True)
