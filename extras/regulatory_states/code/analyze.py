"""Develop GlnD catalytic-state and published six-state GlnE models.

Buffered free inputs and rapid regulatory binding. No new experiments.
The GlnE topology is Jiang, Mayo & Ninfa (2007), Figure 12. Its catalytic
illustration is not a joint fit to all historical kinetic measurements.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results';FIG=ROOT/'figures'
OUT.mkdir(exist_ok=True);FIG.mkdir(exist_ok=True)
SEED=20260913
F=(6.5/.82)/(2.7/2.3)

def base_coefficients(v):
    q0,q2,K=v[0],v[1],np.exp(v[2])
    d2=F*(1-q0)/(1-q2)
    return np.array([q0,d2*q2,d2*q2]),np.array([1-q0,1-q0,d2*(1-q2)]),K

def glnd(G,v,lam=1.):
    a,b,K=base_coefficients(v); a=a.copy();b=b.copy();a[1]*=lam;b[1]*=lam
    g=np.asarray(G)/K
    w=np.stack([np.ones_like(g),2*g/lam,g*g],axis=-1)
    occ=w/w.sum(axis=-1,keepdims=True)
    return occ@a,occ@b,occ

def theta(G,v,lam=1.):
    a,b,_=glnd(G,v,lam);return a/(a+b)

def fit(x,y):
    starts=[[.995,.15,np.log(.08)],[.98,.1,np.log(.2)],[.999,.2,np.log(.03)]]
    fits=[least_squares(lambda v:3*theta(x,v)-y,s,bounds=([.9,.000001,np.log(.001)],[.999999,.499999,np.log(10)]),xtol=1e-12,gtol=1e-12,ftol=1e-12,max_nfev=1800) for s in starts]
    z=min(fits,key=lambda z:z.fun@z.fun);assert z.success
    return z.x,float(z.fun@z.fun)

def glne(G,P,U,par):
    """E, E-Gln, E-PII, E-PIIUMP, E-Gln-PII, E-PII-PIIUMP."""
    g,p,u=np.broadcast_arrays(np.asarray(G)/par['KG'],np.asarray(P)/par['KP'],np.asarray(U)/par['KU'])
    w=np.stack([np.ones_like(g),g,p,u,p*g/par['a1'],p*u/par['a2']],axis=-1)
    occ=w/w.sum(axis=-1,keepdims=True)
    at=par['k']*(occ[...,2]+par['beta']*occ[...,1]+par['betap']*occ[...,4])
    ar=occ[...,3]
    return at,ar,occ

def receptor_stationary(G,P,U,par):
    g,p,u=G/par['KG'],P/par['KP'],U/par['KU']
    Q=np.zeros((6,6))
    for i,j,on in [(0,1,g),(0,2,p),(0,3,u),(1,4,p/par['a1']),(2,4,g/par['a1']),(2,5,u/par['a2']),(3,5,p/par['a2'])]:
        Q[j,i]+=on;Q[i,j]+=1
    Q[np.diag_indices(6)]=-Q.sum(axis=0)
    M=Q.copy();M[-1]=1;b=np.zeros(6);b[-1]=1
    q=np.linalg.solve(M,b)
    return q,float(np.max(np.abs(Q@q)))

def save(name,head,data):
    np.savetxt(OUT/name,np.asarray(data),delimiter=',',header=head,comments='')

def main():
    x,y=np.loadtxt(ROOT/'data/jiangninfa2011_fig2B_digitised.csv',delimiter=',',skiprows=1).T
    v,sse=fit(x,y);G=np.geomspace(.002,30,701);a0,b0,K=base_coefficients(v)
    variants=[1.,2.,4.];drows=[];models=[]
    for lam in variants:
        a,b,occ=glnd(G,v,lam)
        models.append({'lambda':lam,'K1_mM':lam*K/2,'K2_mM':2*K/lam,'K1_over_K2':lam*lam/4,'alpha':(a0*np.array([1,lam,1])).tolist(),'beta':(b0*np.array([1,lam,1])).tolist(),'E1_fraction_at_G_equal_K':float(glnd(K,v,lam)[2][1]),'catalytic_speed_ratio_at_G_equal_K':2*lam/(lam+1)})
        drows.extend(np.column_stack([np.full(len(G),lam),G,3*a/(a+b),a,b,occ]))
        assert np.all(np.diff(a)<=1e-12) and np.all(np.diff(b)>=-1e-12)
    save('glnd_state_families.csv','lambda,G_mM,mean_UMP,UT_capacity,UR_capacity,E0,E1,E2',drows)
    design=[]
    for lam in variants:
        for gg in [K/3,K,3*K]:
            aa,bb,oc=glnd(gg,v,lam)
            design.append([lam,gg,float(oc@np.arange(3)),float(oc[1])])
    save('glnd_occupancy_design.csv','lambda,G_mM,mean_bound_glutamine,E1_fraction',design)
    mean_error=max(float(np.max(np.abs(theta(G,v,l)-theta(G,v)))) for l in variants)
    # Values KG,KP,a1 are published fitted values from the A1/A2 analysis,
    # not a globally validated parameter set. KU is a selected table value.
    # beta,betap and a2 are explicitly chosen illustrations.
    par={'KG':15.6,'KP':.18,'KU':.35,'a1':.17,'a2':4.,'beta':1.,'betap':10.,'k':1.}
    T=5.;gm=.53;thm=theta(gm,v);Pm=T*(1-thm)**3;Um=T*thm
    am,bm,_=glne(gm,Pm,Um,par);par['k']=float(bm/am)
    # Two symmetry directions: binding/catalytic synergy and silent double-bound state.
    epars=[dict(par),dict(par,a1=4*par['a1'],betap=4*par['betap']),dict(par,a2=2.17),dict(par,a2=8.85)]
    names=['reference','binding_catalysis_tradeoff','PII_UMP_interaction_2p17','PII_UMP_interaction_8p85']
    th=theta(G,v);P=T*(1-th)**3;U=T*th
    erows=[];curves=[];rate_ratios=[]
    abase,bbase,_=glne(G,P,U,par)
    for name,pp in zip(names,epars):
        aa,bb,occ=glne(G,P,U,pp);curves.append(aa/(aa+bb));rate_ratios.append((aa+bb)/(abase+bbase))
        erows.extend([[name,*row] for row in np.column_stack([G,aa/(aa+bb),aa,bb,occ])])
    with (OUT/'glne_cascade.csv').open('w') as f:
        f.write('model,G_mM,GS_fraction,AT_capacity,AR_capacity,E,EG,EP,EU,EGP,EPU\n')
        for row in erows:f.write(','.join(map(str,row))+'\n')
    eerror=max(float(np.max(np.abs(c-curves[0]))) for c in curves)
    # Purified-protein condition, distinct from the coupled cascade scenario.
    Pfix,Ufix=1.,.3
    gp=np.geomspace(.02,10,301);a,b,oc=glne(gp,Pfix,Ufix,par)
    q1=[];q2=[];normrate=[]
    for pp in epars:
        aa,bb,oo=glne(gp,Pfix,Ufix,pp)
        q1.append(oo[:,4]);q2.append(oo[:,5]);normrate.append((aa+bb)/(a+b))
    save('glne_purified_predictions.csv','G_mM,EGP_reference,EGP_tradeoff,EGP_a2low,EGP_a2high,rate_ratio_reference,rate_ratio_tradeoff,rate_ratio_a2low,rate_ratio_a2high',np.column_stack([gp,*q1,*normrate]))
    pbind=np.geomspace(.001,10,301);bindcurves=[];binding_half=[]
    for pp in epars[:2]:
        aa,bb,oo=glne(10.,pbind,Ufix,pp);bindcurves.append(oo[:,2]+oo[:,4]+oo[:,5])
        gg,uu=10/pp['KG'],Ufix/pp['KU']
        keff=pp['KP']*(1+gg+uu)/(1+gg/pp['a1']+uu/pp['a2']);binding_half.append(float(keff))
        assert np.max(np.abs(bindcurves[-1]-pbind/(keff+pbind)))<1e-12
    save('glne_PII_binding_isotherms.csv','PII_uM,bound_fraction_reference,bound_fraction_tradeoff',np.column_stack([pbind,*bindcurves]))
    # Explicit concentration signature: all microscopic parameters shared.
    tt=np.geomspace(.5,36,120);ct=[];conc_errors=[]
    for g in [.3,.53,.8]:
        t=float(theta(g,v));pp=tt*(1-t)**3;uu=tt*t
        aa,bb,_=glne(g,pp,uu,par);odds=aa/bb
        interp=np.polyfit(1/tt,odds,1);conc_errors.append(float(np.max(np.abs(np.polyval(interp,1/tt)-odds))))
        lim=par['k']*par['KU']/t*((1-t)**3/par['KP']*(1+par['betap']*g/(par['a1']*par['KG'])))
        ct.extend(np.column_stack([np.full(len(tt),g),tt,odds,np.full(len(tt),lim),odds-lim]))
    save('concentration_signature.csv','G_mM,PII_uM,GS_odds,large_pool_limit,odds_above_large_pool_limit',ct)
    # Exact local sensitivity decomposition; finite difference is independent check.
    h=1e-5
    eta_theta=-(np.log(theta(G*np.exp(h),v))-np.log(theta(G*np.exp(-h),v)))/(2*h)
    eta_p=3*th/(1-th)*eta_theta
    g=G/par['KG'];p=P/par['KP'];u=U/par['KU'];term=np.array([p,par['beta']*g,par['betap']*p*g/par['a1']]).T
    fraction=term/term.sum(axis=1,keepdims=True)
    direct=fraction[:,1]+fraction[:,2]
    via_pii=eta_theta+(fraction[:,0]+fraction[:,2])*eta_p
    def ratio(gx):
        t=theta(gx,v);aa,bb,_=glne(gx,T*(1-t)**3,T*t,par);return aa/bb
    num=(np.log(ratio(G*np.exp(h)))-np.log(ratio(G*np.exp(-h))))/(2*h)
    deriv_error=float(np.max(np.abs(num-direct-via_pii)))
    assert direct.min()>=0 and direct.max()<=1 and deriv_error<1e-6
    save('local_sensitivity.csv','G_mM,total_log_ratio_sensitivity,via_PII,direct_glutamine',np.column_stack([G,num,via_pii,direct]))
    # Cascaded dynamical illustration: common rate units, fixed free-target pool.
    gstart,gstop=.5,.05
    ts=np.linspace(0,8,501);init_theta=float(theta(gstart,v));ai,bi,_=glne(gstart,T*(1-init_theta)**3,T*init_theta,par);init_y=float(ai/(ai+bi))
    ad,bd,_=glnd(gstop,v);sd=1/float(ad+bd)
    tf=float(theta(gstop,v));ae,be,_=glne(gstop,T*(1-tf)**3,T*tf,par);se=1/float(ae+be)
    dynamics=[];solutions=[]
    for label,lam,pp in [('reference',1,par),('GlnD_occupancy',4,par),('GlnE_occupancy',1,epars[1]),('both',4,epars[1])]:
        aa,bb,_=glnd(gstop,v,lam)
        def rhs(t,z):
            zz=float(np.clip(z[0],0,1));at,ar,_=glne(gstop,T*(1-zz)**3,T*zz,pp)
            return [sd*(aa*(1-z[0])-bb*z[0]),se*(at*(1-z[1])-ar*z[1])]
        sol=solve_ivp(rhs,[0,ts[-1]],[init_theta,init_y],t_eval=ts,rtol=1e-9,atol=1e-11)
        assert sol.success and sol.y.min()>-1e-10 and sol.y.max()<1+1e-10
        solutions.append(sol.y)
        dynamics.extend([[label,*row] for row in np.column_stack([ts,sol.y.T])])
    with (OUT/'step_responses.csv').open('w') as f:
        f.write('model,time_reference_units,PII_fraction,GS_fraction\n')
        for row in dynamics:f.write(','.join(map(str,row))+'\n')
    # Random independent equilibrium generator solves and gauge checks.
    rng=np.random.default_rng(SEED);stationary_error=0.;balance_error=0.;identity_error=0.
    for j in range(250):
        pp=dict(par,a1=np.exp(rng.uniform(-4,1)),a2=np.exp(rng.uniform(-1,3)),betap=np.exp(rng.uniform(-2,4)))
        gg,ppii,uu=np.exp(rng.uniform(-4,2,3));mult=np.exp(rng.uniform(-1,1))
        at,ar,oc=glne(gg,ppii,uu,pp);q,res=receptor_stationary(gg,ppii,uu,pp)
        a2,b2,_=glne(gg,ppii,uu,dict(pp,a1=pp['a1']*mult,betap=pp['betap']*mult))
        stationary_error=max(stationary_error,float(np.max(np.abs(q-oc))));balance_error=max(balance_error,res)
        identity_error=max(identity_error,float(abs(at/(at+ar)-a2/(a2+b2))))
    coupled_error=0.;glnd_stationary_error=0.
    for total in [.5,5,36]:
        thref=theta(G,v);a,b,_=glne(G,total*(1-thref)**3,total*thref,par);ref=a/(a+b)
        for lam in variants:
            thl=theta(G,v,lam)
            for pp in epars:
                a,b,_=glne(G,total*(1-thl)**3,total*thl,pp)
                coupled_error=max(coupled_error,float(np.max(np.abs(a/(a+b)-ref))))
    for j in range(100):
        gg=np.exp(rng.uniform(-5,2));lam=rng.uniform(1,F);a,b,_=glnd(gg,v,lam)
        Q=np.zeros((4,4))
        for i in range(3):Q[i+1,i]=(3-i)*a;Q[i,i+1]=(i+1)*b
        Q[np.diag_indices(4)]=-Q.sum(axis=0);M=Q.copy();M[-1]=1;target=np.array([0.,0.,0.,1.]);sol=np.linalg.solve(M,target)
        t=float(theta(gg,v));exact=np.array([(1-t)**3,3*t*(1-t)**2,3*t*t*(1-t),t**3])
        glnd_stationary_error=max(glnd_stationary_error,float(np.max(np.abs(sol-exact))))
    assert max(mean_error,eerror,stationary_error,balance_error,identity_error,coupled_error,glnd_stationary_error)<1e-10
    report={'seed':SEED,'scope':'Free-input reaction models; state binding at rapid equilibrium; no closed-total fit or absolute biological time prediction.',
      'glnd':{'fitted_parameters_q0_q2_logK':v.tolist(),'RMSE':float(np.sqrt(sse/len(x))),'fitted_parameters_count':3,'transported_UR_fold':F,'family':models,'mean_invariance_error':mean_error,'lambda_monotonicity_interval':[1,F]},
      'glne':{'topology_source':'Jiang, Mayo & Ninfa 2007, Fig 12, DOI 10.1021/bi0620510','parameter_scenarios':dict(zip(names,epars)),'parameter_provenance':'KG=15.6 mM, KP=0.18 uM, alpha1=0.17: source reported joint-fit example; KU=0.35 uM: source Table 3 selected apparent scale; beta=1 and beta_prime=10: chosen illustration; alpha2 variants 2.17,4,8.85: illustrative range motivated by source fits, not CI. k calibrated once to GS=0.5 at G=0.53mM,T=5uM.', 'steady_response_invariance_error':eerror,'purified_condition':{'PII_uM':Pfix,'PII_UMP_uM':Ufix,'G10_rate_ratio_tradeoff':float(normrate[1][-1]),'G10_EGP_reference':float(q1[0][-1]),'G10_EGP_tradeoff':float(q1[1][-1]),'G10_binding_halfpoints_uM':binding_half,'binding_halfpoint_ratio':binding_half[1]/binding_half[0]},'concentration_line_error':max(conc_errors),'local_decomposition_error':deriv_error,'direct_local_range':[float(direct.min()),float(direct.max())]},
      'verification':{'random_six_state_networks':250,'stationary_linear_solve_error':stationary_error,'species_balance_error':balance_error,'random_synergy_identity_error':identity_error,'coupled_36_scenarios_error':coupled_error,'random_GlnD_stationary_solves':100,'GlnD_stationary_solve_error':glnd_stationary_error},
      'dynamics':{'G_start_mM':gstart,'G_end_mM':gstop,'PII_uM':T,'time_units':'reference relaxation scales set at final input; illustrative, not minutes','max_GS_separation':float(np.max(np.abs(solutions[-1][1]-solutions[0][1])))}}
    (OUT/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    colors=['#477C9B','#B77A75','#7B9778','#A59062']
    def finish(fig,axs,name):
        for ax,l in zip(axs.flat,'ABCD'):ax.text(-.16,1.04,l,transform=ax.transAxes,fontweight='bold',fontsize=12)
        fig.savefig(FIG/(name+'.pdf'));fig.savefig(FIG/(name+'.png'),dpi=600);plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(7.3,5.7),layout='constrained')
    ax=axs[0,0];ax.scatter(x,y,c='black',s=15,label='Digitized mean',zorder=5)
    for lam,c,ls in zip(variants,colors,['-','--',':']):ax.plot(G,3*theta(G,v,lam),c=c,ls=ls,label=f'λ = {lam:g}')
    ax.set(xscale='log',xlim=(.015,15),xlabel='Free glutamine (mM)',ylabel='UMP groups per trimer');ax.legend(fontsize=7)
    ax=axs[0,1]
    for lam,c in zip(variants,colors):ax.plot(G,glnd(G,v,lam)[2][:,1],c=c,label=f'λ = {lam:g}')
    ax.set(xscale='log',xlim=(.002,10),xlabel='Free glutamine (mM)',ylabel='Singly liganded GlnD fraction');ax.legend(fontsize=7)
    ax=axs[1,0];ar,br,_=glnd(G,v)
    for lam,c in zip(variants,colors):
        aa,bb,_=glnd(G,v,lam);ax.plot(G,(aa+bb)/(ar+br),c=c,label=f'λ = {lam:g}')
    ax.set(xscale='log',xlim=(.002,10),xlabel='Free glutamine (mM)',ylabel='PII relaxation-rate ratio');ax.legend(fontsize=7)
    ax=axs[1,1]
    for i,label in [(0,'Reference'),(1,'GlnD occupancy changed')]:ax.plot(ts,solutions[i][0],c=colors[i],label=label)
    ax.set(xlabel='Time (reference units)',ylabel='Mean PII modification fraction');ax.legend(fontsize=7)
    finish(fig,axs,'Fig7_GlnD_states')
    fig,axs=plt.subplots(2,2,figsize=(7.3,5.7),layout='constrained')
    ax=axs[0,0]
    for i,(cv,c,ls) in enumerate(zip(curves,colors,['-','--',':','-.'])):ax.plot(G,cv,c=c,ls=ls,label=['Reference','α₁ and β′ × 4','α₂ = 2.17','α₂ = 8.85'][i])
    ax.set(xscale='log',xlim=(.02,10),xlabel='Free glutamine (mM)',ylabel='GS adenylylated fraction');ax.legend(fontsize=7)
    ax=axs[0,1]
    for i in [0,1]:ax.plot(pbind,bindcurves[i],c=colors[i],label=['Reference','α₁ and β′ × 4'][i])
    ax.axhline(.5,color='0.7',ls=':',lw=.7)
    ax.set(xscale='log',xlabel='Free unmodified PII (μM)',ylabel='PII-bound GlnE fraction');ax.legend(fontsize=7)
    ax=axs[1,0]
    for i,c in enumerate(colors):ax.plot(gp,normrate[i],c=c,label=['Reference','α₁ and β′ × 4','α₂ = 2.17','α₂ = 8.85'][i])
    ax.set(xscale='log',xlabel='Free glutamine (mM)',ylabel='GS relaxation-rate ratio');ax.legend(fontsize=7)
    ax=axs[1,1];ct=np.array(ct)
    for gg,c in zip([.3,.53,.8],colors):
        sub=ct[ct[:,0]==gg];ax.plot(1/sub[:,1],sub[:,4],c=c,label=f'{gg:g} mM Gln')
    ax.set(xlabel='1 / free PII total (μM⁻¹)',ylabel='GS odds above large-pool limit');ax.legend(fontsize=7)
    finish(fig,axs,'Fig8_GlnE_states')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
