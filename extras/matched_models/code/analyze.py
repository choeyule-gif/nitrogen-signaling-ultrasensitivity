"""Matched-data mechanistic alternatives. No new experimental observations.

Run: python code/analyze.py. Inputs and generated results are separate.
Exact mean-matched ladders are structural counterexamples, not finite-dimensional
biochemical fits. Finite-switch fits are separately reported. Cascade summaries
are compared without invented likelihoods or replicate confidence intervals.
"""
from pathlib import Path
import json, math, argparse
import numpy as np
from scipy.optimize import least_squares, brentq, differential_evolution
from scipy.special import expit, logsumexp, gammaln
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results'
OUT.mkdir(exist_ok=True)
SEED=20260913

def savecsv(name, header, a):
    np.savetxt(OUT/name, np.asarray(a), delimiter=',', header=header, comments='')

def hill(x, p):
    top, bottom, logmid, h=p
    x=np.asarray(x)
    return bottom+(top-bottom)*expit(-h*(np.log(np.maximum(x,1e-300))-logmid))

class Ladder:
    def __init__(self,n=3,J=0):
        self.n=n; self.i=np.arange(n+1)
        self.logc=gammaln(n+1)-gammaln(self.i+1)-gammaln(n-self.i+1)+J*((self.i-n/2)**2-(n/2)**2)
        z=np.linspace(-32,32,16001)
        q=self.mean(z)
        uq,ix=np.unique(q,return_index=True)
        self.inv=PchipInterpolator(uq, z[ix], extrapolate=True)
    def probs(self,z):
        w=np.asarray(z)[...,None]*self.i+self.logc
        return np.exp(w-logsumexp(w,axis=-1,keepdims=True))
    def mean(self,z):
        return self.probs(z)@self.i/self.n
    def inverse(self,q):
        return self.inv(np.clip(q,1e-12,1-1e-12))
    def exact_inverse(self,q):
        return brentq(lambda z:float(self.mean(z))-q,-40,40,xtol=1e-13)

LO=np.array([1.50001,0.00001,np.log(.005),.1])
HI=np.array([2.999999,1.49999,np.log(30),8.])

def finite_switch(x,p,lad):
    """Saturating effective u with exponent h, midpoint of mean response fixed at S.
    Equivalent to a basal saturating f=1/u, but not a claim of ligand stoichiometry.
    """
    q0,qf=p[0]/3,p[1]/3
    z0,zf,zm=lad.inverse(np.array([q0,qf,(q0+qf)/2]))
    u0,uf,um=np.exp([z0,zf,zm])
    logt=np.log((u0-um)/(um-uf))+p[3]*(np.log(np.maximum(x,1e-300))-p[2])
    u=uf+(u0-uf)*expit(-logt)
    return 3*lad.mean(np.log(u))

def fit(x,y,lad=None,start=None,multistart=5):
    fun=hill if lad is None else lambda x,p:finite_switch(x,p,lad)
    rng=np.random.default_rng(SEED)
    starts=[np.array([2.968,.454,np.log(.56),2.018]) if start is None else np.array(start)]
    for _ in range(multistart-1):
        starts.append(np.array([rng.uniform(2.7,2.999),rng.uniform(.1,.7),np.log(rng.uniform(.3,.8)),rng.uniform(.8,3)]))
    fits=[least_squares(lambda p:fun(x,p)-y,np.clip(s,LO+1e-8,HI-1e-8),bounds=(LO,HI),xtol=1e-11,ftol=1e-11,gtol=1e-11,max_nfev=1200) for s in starts]
    best=min(fits,key=lambda z:np.sum(z.fun*z.fun))
    return best.x,float(best.fun@best.fun),bool(best.success)

def distribution(x,p,J=0,refractory=False):
    q=hill(x,p)/3
    if refractory:
        r=p[1]/3
        qa=np.clip((q-r)/(1-r),0,1)
        probs=np.stack([(1-qa)**3,3*qa*(1-qa)**2,3*qa**2*(1-qa),qa**3],axis=-1)*(1-r)
        probs[...,3]+=r
        return probs
    if abs(J)<1e-12:
        return np.stack([(1-q)**3,3*q*(1-q)**2,3*q*q*(1-q),q**3],axis=-1)
    lad=Ladder(3,J)
    return lad.probs(lad.inverse(q))

def prepare_input(p,J=0,refractory=False,grid=3001):
    lx=np.linspace(np.log(1e-8),np.log(1e6),grid)
    x=np.exp(lx)
    q=hill(x,p)/3
    prob=distribution(x,p,J,refractory)
    base=np.log(np.maximum(prob[:,0],1e-300)/q)
    b0=np.log(float(distribution(0,p,J,refractory)[0])/(p[0]/3))
    bf=np.log(float(distribution(1e100,p,J,refractory)[0])/(p[1]/3))
    return lx,base,b0,bf

def cascade_score(inp,mids,kind='baseline',pars=(),return_curves=False):
    lx,base,b0,bf=inp
    if kind=='direct':
        R,KG=np.exp(pars)
        logphi=np.logaddexp(0,np.log(R)+lx-np.log(KG))-np.logaddexp(0,lx-np.log(KG))
        base=base+logphi; bf=bf+np.log(R)
    if kind=='power':
        exponent=float(np.exp(pars[0]))
    else: exponent=1.
    if kind in ('gs_ladder','gs_asym'):
        n=12; i=np.arange(n+1)
        lc=gammaln(n+1)-gammaln(i+1)-gammaln(n-i+1)+pars[0]*((i-n/2)**2-(n/2)**2)
        if kind=='gs_asym':lc=lc+pars[1]*(i-n/2)**3/(n/2)
        def response(z):
            w=np.asarray(z)[...,None]*i+lc
            return np.exp(w-logsumexp(w,axis=-1,keepdims=True))@i/n
    else:
        response=lambda z:expit(exponent*np.asarray(z))
    ans=[]; kappas=[]; curves=[]
    for m in mids:
        bm=np.interp(np.log(m),lx,base)
        def cal(k):
            a,b,c=response(np.array([k+b0,k+bm,k+bf]))
            # scale to avoid spurious roots at numerical saturation
            return (b-a)/(c-a)-.5 if c-a>1e-15 else np.nan
        ks=np.linspace(-35,35,141)
        rr=response(ks[:,None]+np.array([b0,bm,bf]))
        span=rr[:,2]-rr[:,0]
        vals=np.full(len(ks),np.nan)
        good=span>1e-15
        vals[good]=(rr[good,1]-rr[good,0])/span[good]-.5
        cross=np.flatnonzero(np.isfinite(vals[:-1]) & np.isfinite(vals[1:]) & (vals[:-1]*vals[1:]<0))
        if len(cross)!=1: raise ValueError(f'calibration roots {len(cross)}')
        ix=cross[0]; k=brentq(cal,ks[ix],ks[ix+1],xtol=1e-11)
        end=response(np.array([k+b0,k+bf])); y=(response(k+base)-end[0])/(end[1]-end[0])
        if np.min(np.diff(y))<-1e-10: raise ValueError('nonmonotone response')
        l10=np.interp(.1,y,lx); l90=np.interp(.9,y,lx)
        ans.append(np.log(81)/(l90-l10)); kappas.append(np.exp(k));curves.append(y)
    if return_curves: return np.array(ans),kappas,np.array(curves)
    return np.array(ans)

def main(bootstrap=300):
    data=np.loadtxt(ROOT/'data/jiangninfa2011_fig2B_digitised.csv',delimiter=',',skiprows=1)
    x,y=data.T; p,sse,success=fit(x,y,multistart=12)
    report={'seed':SEED,'n_coordinates':len(x),'reference':{'parameters':p.tolist(),'SSE':sse,'success':success},'finite_switch':[], 'matched':{},'cascade':[]}
    xplot=np.geomspace(.005,30,600)
    rows=[]
    for J in [-.7,0,.7]:
        lad=Ladder(3,J)
        fp,ss,ok=fit(x,y,lad,multistart=8)
        k=5;n=len(x)
        aicc=n*np.log(ss/n)+2*k+2*k*(k+1)/(n-k-1)
        report['finite_switch'].append({'J':J,'a':float(np.exp(-2*J)),'parameters':fp.tolist(),'SSE':ss,'RMSE':float(np.sqrt(ss/n)),'AICc':aicc,'success':ok})
        rows.append(finite_switch(xplot,fp,lad))
    savecsv('finite_switch_curves.csv','G_mM,J_minus_0p7,J_0,J_plus_0p7',np.column_stack([xplot]+rows))
    # Profile fixed ladder coupling, counting all fits as a profile not selected models.
    profile=[]
    for J in np.linspace(-1.5,1.5,61):
        fp,ss,ok=fit(x,y,Ladder(3,J),multistart=2)
        profile.append([J,ss,len(x)*np.log(ss/sse),*fp,ok])
    savecsv('finite_switch_profile.csv','J,SSE,loss_vs_reference,top,bottom,logmid,h,success',profile)
    for name,J,refr in [('independent',0,False),('extreme_enriched',.7,False),('central_enriched',-.7,False),('refractory',0,True)]:
        probs=distribution(xplot,p,J,refr)
        mean=probs@np.arange(4)
        variance=probs@(np.arange(4)**2)-mean**2
        xc=float(np.exp(p[2])*((p[0]-1.5)/(1.5-p[1]))**(1/p[3]))
        pc=distribution(xc,p,J,refr)
        report['matched'][name]={'J':J,'refractory':refr,'max_mean_error':float(np.max(abs(mean-hill(xplot,p)))), 'G_at_mean_half_mM':xc,'p_at_mean_half':pc.tolist(),'variance_at_mean_half':float(pc@(np.arange(4)-1.5)**2)}
        savecsv('matched_'+name+'.csv','G_mM,U,p0,p1,p2,p3,variance',np.column_stack([xplot,mean,probs,variance]))
    # Conditional curves all use the same mean reference. Midpoints are nuisance calibrations.
    tab=np.loadtxt(ROOT/'data/cascade_comparison.csv',delimiter=',',skiprows=1)
    pii,mids,obs=tab[:,:3].T
    for name,m in report['matched'].items():
        inp=prepare_input(p,m['J'],m['refractory'])
        for kind in ['baseline','direct','power','gs_ladder']:
            bounds={'direct':[(0,np.log(100)),(np.log(.001),np.log(100))],'power':[(np.log(.3),np.log(5))],'gs_ladder':[(-.15,.8)]}
            if kind=='baseline': pars=np.array([]); objective=float(np.sum(np.log(cascade_score(inp,mids)/obs)**2))
            else:
                def obj(v):
                    try:return float(np.sum(np.log(cascade_score(inp,mids,kind,v)/obs)**2))
                    except ValueError:return 1e6
                result=differential_evolution(obj,bounds[kind],seed=SEED,popsize=7,maxiter=45,tol=1e-7,polish=True)
                pars=result.x;objective=result.fun
            pred,kappa,curves=cascade_score(inp,mids,kind,pars,True)
            entry={'upstream':name,'kind':kind,'shared_parameters':pars.tolist(),'objective_log_summary_SSE':objective,'predicted':pred.tolist(),'reported':obs.tolist(),'midpoints_mM':mids.tolist(),'PII_uM':pii.tolist(),'kappa':kappa,'max_abs_relative_error':float(np.max(abs(pred/obs-1)))}
            report['cascade'].append(entry)
            savecsv(f'cascade_{name}_{kind}.csv','G_mM,Y_0p5,Y_5,Y_36',np.column_stack([np.exp(inp[0]),curves.T]))
            print(name,kind,np.round(pred,4),'max relative error',round(entry['max_abs_relative_error'],4),flush=True)
    # Joint uncertainty scenarios: same coordinate draws for all structural alternatives.
    rng=np.random.default_rng(SEED)
    draws=[];distdraw=[];casdraw=[]
    for b in range(bootstrap):
        xb=x*np.exp(rng.normal(0,.06,len(x)));yb=y+rng.normal(0,.06,len(y))
        pb,ss,ok=fit(xb,yb,start=p,multistart=1)
        if not ok: raise RuntimeError('perturbed fit failed')
        draws.append(pb)
        for j,(name,m) in enumerate(report['matched'].items()):
            q=distribution(.7,pb,m['J'],m['refractory'])
            distdraw.append([b,j,*q])
            inp=prepare_input(pb,m['J'],m['refractory'],grid=1201)
            pred=cascade_score(inp,mids)
            casdraw.append([b,j,*pred])
    savecsv('coordinate_draws.csv','top,bottom,logmid,h',draws)
    savecsv('distribution_draws_at_0p7mM.csv','draw,model,p0,p1,p2,p3',distdraw)
    savecsv('cascade_baseline_draws.csv','draw,model,n_0p5,n_5,n_36',casdraw)
    report['coordinate_sensitivity']={'draws':bootstrap,'sigma_logx':.06,'sigma_U':.06,'h_percentiles':np.percentile(np.array(draws)[:,3],[2.5,50,97.5]).tolist()}
    (OUT/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print('DONE',OUT/'summary.json',flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--draws',type=int,default=300)
    main(ap.parse_args().draws)
