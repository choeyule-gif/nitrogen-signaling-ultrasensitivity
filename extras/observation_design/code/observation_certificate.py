"""A11 certificate anchored by an observable PII count threshold.
Uniform in downstream catalytic gain, and no independent PII-site assumption.
The reported titration motivates a gate; the gate must be checked in the assay.
"""
from pathlib import Path
import json
import numpy as np
from scipy.special import logsumexp
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results';OUT.mkdir(exist_ok=True)

def certificate(Pt,St,Dt,Et,m_upper,G,KG,KP,alpha,beta,beta_prime,nu=12):
 eps=(Dt+2*Et)/Pt;es=Et/St
 if eps>=1:return dict(bound=1.,p0_lower=0.,epsilon_P=eps,epsilon_S=es)
 p0=max(0.,1-m_upper/(1-eps))
 # Effective concentration ell = B/A; infinity recovers the generic bound.
 ell=np.inf if p0==0 else beta*(G/KG)*KP/(p0*(1+beta_prime*(G/KG)/alpha))
 logratio=-np.log1p(-eps) if not np.isfinite(ell) else np.log1p(ell/(Pt*(1-eps)))-np.log1p(ell/Pt)
 return dict(bound=float(min(1.,np.tanh(nu*logratio/4)+es)),generic_bound=float(min(1.,np.tanh(nu*(-np.log1p(-eps))/4)+es)),p0_lower=p0,ell_upper_uM=float(ell) if np.isfinite(ell) else None,epsilon_P=eps,epsilon_S=es,log_drive_ratio=float(logratio))

def mean(n,logr,coef):
 j=np.arange(n+1);w=np.log(coef)+j*logr;pi=np.exp(w-logsumexp(w));return float(pi@j/n)

def main():
 args=dict(Pt=5.,St=.5,Dt=.5,Et=.02,m_upper=.67,G=2.,KG=15.6,KP=.18,alpha=.17,beta=175/7.5,beta_prime=1315/7.5)
 results={}
 for label,changes in [('source_fit_units',{}),('literal_SI_KP_units',{'KP':180.}),('apparent_scale_sensitivity',{'KP':1.}),('original_E_concentration',{'Et':.1}),('no_PII_observation',{'m_upper':3.})]:results[label]=certificate(**(args|changes))
 r=results['source_fit_units'];tolerance=.05;L=args['Pt']*(1-r['epsilon_P']);U=args['Pt'];R=np.exp(4*np.arctanh(tolerance-r['epsilon_S'])/12)
 ell_limit=(R-1)/(1/L-R/U)
 # Independent positive ladder evaluation over unknown catalytic gains.
 rng=np.random.default_rng(9914);worst=0.
 for _ in range(5000):
  n=int(rng.integers(1,13));coef=np.exp(rng.uniform(-8,8,n+1));t1,t2=rng.uniform(L,U,2);ell=rng.uniform(0,r['ell_upper_uM']);k=np.exp(rng.uniform(-10,10))
  gap=abs(mean(n,np.log(k*(1+ell/t1)),coef)-mean(n,np.log(k*(1+ell/t2)),coef))
  bound=np.tanh(n*r['log_drive_ratio']/4)
  worst=max(worst,gap-bound);assert gap<=bound+1e-11
 # Independent mixtures verify the free zero-count lower bound from measured totals.
 mixture_error=0.
 for _ in range(5000):
  pf=rng.dirichlet(np.ones(4));pb=rng.dirichlet(np.ones(4));s=rng.uniform(0,r['epsilon_P']);mt=(1-s)*(pf@np.arange(4))+s*(pb@np.arange(4))
  lower=max(0,1-mt/(1-r['epsilon_P']));mixture_error=max(mixture_error,lower-pf[0]);assert lower<=pf[0]+1e-12
 report=dict(parameters=args,scenarios=results,ell_limit_for_five_percent=ell_limit,nominal_coefficient_headroom=ell_limit/r['ell_upper_uM'],random_ladder_checks=5000,random_mixture_checks=5000,max_ladder_violation=worst,max_mixture_violation=mixture_error,scope='A prospective observation gate m_total<=.67, not a confidence interval inferred from rounded digitization. GlnE source-fit coefficients follow main-text KP units; literal SI units shown separately. Source assay rates are not transported as a validated full-cascade physiological fit.')
 (OUT/'observation_certificate.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
