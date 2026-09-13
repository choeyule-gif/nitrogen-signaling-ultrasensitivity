"""Finite biochemical counterexample and explicit transportability stress tests.

No new observations, no microscopic constants inferred from apparent constants.
Run with Python, NumPy, SciPy, Matplotlib. Outputs are deterministic.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import least_squares, brentq
from scipy.special import comb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
FIG=ROOT/'figures'; FIG.mkdir(exist_ok=True)
SEED=20260913

def unpack(p):
    q0,q2,t,logs,logc=p
    return np.array([q0,q2+t*(q0-q2),q2]),np.exp(logs),np.exp(logc)

def mean(x,p):
    q,s,c=unpack(p)
    z=np.asarray(x)/s
    w=np.stack([np.ones_like(z),z/np.sqrt(c),z*z],axis=-1)
    return (w@q)/w.sum(axis=-1)

def fit(x,y,cap=100.,fixed_c=None,starts=12):
    lo=np.array([.500001,.000001,.000001,np.log(.005),np.log(.0025)])
    hi=np.array([.999999,.499999,.999999,np.log(30),np.log(cap)])
    if fixed_c is not None:
        lo[4]=np.log(fixed_c)-1e-10; hi[4]=np.log(fixed_c)+1e-10
    rng=np.random.default_rng(SEED)
    pp=[np.array([.989,.151,.5,np.log(.56),np.log(min(cap,10.))])]
    for j in range(starts-1):
        pp.append(np.array([rng.uniform(.92,.999),rng.uniform(.10,.20),rng.uniform(.05,.95),np.log(rng.uniform(.35,.75)),rng.uniform(lo[4],hi[4])]))
    fits=[least_squares(lambda p:3*mean(x,p)-y,np.clip(p,lo+1e-12,hi-1e-12),bounds=(lo,hi),xtol=1e-11,gtol=1e-11,ftol=1e-11,max_nfev=1800) for p in pp]
    best=min(fits,key=lambda z:z.fun@z.fun)
    return best.x,float(best.fun@best.fun),bool(best.success)

def fit_transport(x,y,C=100.,half=.080,fold=(6.5/.82)/(2.7/2.3)):
    """Stress test: transport UR half-range and specificity fold across assays.
    Two fitted parameters (q2,s); q1=q0 is a fixed double-occupancy topology.
    Values are not claimed as shared microscopic binding constants.
    """
    def transform(v):
        q2,s=v[0],np.exp(v[1])
        d2=half**2/(s*s+s*half/np.sqrt(C))
        q0=1-d2*(1-q2)/fold
        return np.array([q0,q2,1.,np.log(s),np.log(C)])
    starts=[[.13,np.log(.5)],[.05,np.log(.25)],[.2,np.log(1.)]]
    fits=[least_squares(lambda v:3*mean(x,transform(v))-y,v,bounds=([.000001,np.log(.1)],[.499999,np.log(3)]),xtol=1e-12,gtol=1e-12,ftol=1e-12,max_nfev=1500) for v in starts]
    z=min(fits,key=lambda z:z.fun@z.fun)
    return transform(z.x),float(z.fun@z.fun),bool(z.success)

def fit_independent(x,y,fold=(6.5/.82)/(2.7/2.3)):
    """Independent ligand binding, sequential activity gating, monotone rates.
    E1 has the low UT capacity of E2 and the basal UR capacity of E0.
    Three fitted parameters; UR specificity fold is a transported scenario.
    """
    def transform(v):
        q0,q2,K=v[0],v[1],np.exp(v[2])
        d2=fold*(1-q0)/(1-q2);d1=d2*q2+(1-q0)
        q1=d2*q2/d1;s=K/np.sqrt(d2);C=d2/(4*d1*d1)
        return np.array([q0,q2,(q1-q2)/(q0-q2),np.log(s),np.log(C)])
    starts=[[.995,.15,np.log(.08)],[.98,.1,np.log(.2)],[.999,.2,np.log(.03)]]
    fits=[least_squares(lambda v:3*mean(x,transform(v))-y,v,bounds=([.9,.000001,np.log(.001)],[.999999,.499999,np.log(10)]),xtol=1e-12,gtol=1e-12,ftol=1e-12,max_nfev=1800) for v in starts]
    z=min(fits,key=lambda z:z.fun@z.fun)
    return transform(z.x),float(z.fun@z.fun),bool(z.success)

def dist(theta,r=0):
    q=(np.asarray(theta)-r)/(1-r)
    assert np.min(q)>-1e-10 and np.max(q)<1+1e-10
    q=np.clip(q,0,1)
    d=np.stack([comb(3,i)*q**i*(1-q)**(3-i) for i in range(4)],axis=-1)*(1-r)
    d[...,3]+=r
    return d

def enzyme_rates(x,p,fold=6.5/.82/(2.7/2.3)):
    """A choice of catalytic scale makes UR specificity fold match Mg Table 2.
    Binding polynomial changes accordingly; mean curve remains unchanged.
    This is a cross-assay scenario, not a joint kinetic fit.
    """
    q,s,c=unpack(p)
    d2=fold*(1-q[0])/(1-q[2])
    # Equal UR capacities in E0 and E1 guarantee monotone UR activation;
    # monotone q_j then guarantees decreasing UT capacity as G increases.
    d=np.array([1,(1-q[0])/(1-q[1]),d2])
    alpha=d*q; beta=d*(1-q)
    r=q[2]; ar=(1-r)*alpha-r*beta
    assert np.min(ar)>-1e-14
    z=np.asarray(x)/s
    w=np.stack([np.ones_like(z),z/np.sqrt(c),z*z],axis=-1)/d
    den=w.sum(axis=-1)
    return w@alpha/den,w@beta/den,w@ar/den,alpha,beta,ar,d

def stationary(alpha,beta):
    Q=np.zeros((4,4))
    for i in range(3):
        Q[i+1,i]=(3-i)*alpha; Q[i,i+1]=(i+1)*beta
    for i in range(4):Q[i,i]=-Q[:,i].sum()
    M=Q.copy(); M[-1]=1.; b=np.array([0.,0.,0.,1.])
    return np.linalg.solve(M,b),Q

def save(name,header,data):
    np.savetxt(OUT/name,np.asarray(data),delimiter=',',header=header,comments='')

def main():
    x,y=np.loadtxt(ROOT/'data/jiangninfa2011_fig2B_digitised.csv',delimiter=',',skiprows=1).T
    profile=[]; cases=[]
    for cap in [.25,1,10,100,10000]:
        p,sse,ok=fit(x,y,cap)
        q,s,c=unpack(p)
        row={'cap_C':cap,'parameters':p.tolist(),'q_states':q.tolist(),'s_mM':s,'C':c,'SSE':sse,'RMSE':np.sqrt(sse/len(x)),'success':ok,'C_upper_boundary':bool(c>=cap*.999)}
        cases.append(row); profile.append([cap,c,sse,row['RMSE'],*q,s])
    # Fixed, finite C=100: no continuum inverse and no optimized cooperativity.
    freep,freesse,freeok=fit(x,y,cap=100,fixed_c=100)
    transportp,transportsse,transportok=fit_transport(x,y)
    p,sse,ok=fit_independent(x,y)
    q,s,c=unpack(p); r=q[2]
    loo=[]
    for i in range(len(x)):
        mask=np.arange(len(x))!=i
        pp,ss,suc=fit_independent(x[mask],y[mask])
        loo.append([x[i],y[i],3*mean(x[i],pp),int(suc)])
    xx=np.geomspace(.002,100,1001); theta=mean(xx,p)
    ph=dist(theta); pr=dist(theta,r)
    a,b,ar,alpha,beta,alphaR,d=enzyme_rates(xx,p)
    eqerr=float(np.max(np.abs(ph@np.arange(4)-pr@np.arange(4))))
    gain=(1-r)**2
    verr=float(np.max(np.abs(ph[:,0]/theta-gain*pr[:,0]/theta)))
    kineticerr=float(np.max(np.abs(a/(a+b)-theta)))
    revkineticerr=float(np.max(np.abs(r+(1-r)*ar/(ar+b)-theta)))
    staterr=0.; residual=0.
    for j in np.linspace(0,len(xx)-1,41,dtype=int):
        sh,Q=stationary(a[j],b[j]); sr,QR=stationary(ar[j],b[j])
        mr=sr*(1-r);mr[3]+=r
        staterr=max(staterr,float(np.max(np.abs(sh-ph[j]))),float(np.max(np.abs(mr-pr[j]))))
        residual=max(residual,float(np.max(np.abs(Q@sh))),float(np.max(np.abs(QR@sr))))
    assert max(eqerr,verr,kineticerr,revkineticerr,staterr,residual)<1e-10
    save('cooperativity_profile.csv','C_cap,C_fit,SSE,RMSE,q0,q1,q2,s_mM',profile)
    transferprofile=[]
    for cap in [.25,1,10,100,10000]:
        tp,ts,tok=fit_transport(x,y,cap)
        tq,tt,tc=unpack(tp)
        transferprofile.append([cap,ts,np.sqrt(ts/len(x)),*tq,tt,int(tok)])
    save('transported_kinetics_profile.csv','effective_C,SSE,RMSE,q0,q1,q2,s_mM,success',transferprofile)
    # Do not cherry-pick successful transfers: transporting the unliganded
    # UT/UR specificity ratio as well strongly worsens this minimal model.
    ratio0=(137./3.)/(2.7/2.3)
    fixed_q0=ratio0/(1+ratio0); fold0=(6.5/.82)/(2.7/2.3)
    def strict_parameters(v):
        q2=v[0]; d2=fold0*(1-fixed_q0)/(1-q2)
        s0=(-.08/10+np.sqrt((.08/10)**2+4*.08**2/d2))/2
        return np.array([fixed_q0,q2,1.,np.log(s0),np.log(100.)])
    strictfit=least_squares(lambda v:3*mean(x,strict_parameters(v))-y,[.15],bounds=([.000001],[.499999]),xtol=1e-12,gtol=1e-12,ftol=1e-12)
    strictp=strict_parameters(strictfit.x)
    save('strict_cross_assay_transfer.csv','G_mM,U_observed,U_strict_transfer',np.column_stack([x,y,3*mean(x,strictp)]))
    save('leave_one_out.csv','G_mM,U_observed,U_predicted,success',loo)
    save('finite_pair_curves.csv','G_mM,theta,p0_H,p1_H,p2_H,p3_H,p0_R,p1_R,p2_R,p3_R,UT_H,UR,UT_R',np.column_stack([xx,theta,ph,pr,a,b,ar]))
    # The finite model does not include a universal, unmeasured gain range.
    # Deliberately test a topology that breaks the common-denominator identity.
    # Kp/Kq here are effective half-saturation scenario inputs, not measured Kd.
    kpgrid=np.geomspace(.015,2.5,17); kqgrid=np.geomspace(.28,.85,13)
    gx=np.geomspace(.02,10,401); th=mean(gx,p); hp=dist(th)[:,0]; rp=dist(th,r)[:,0]
    gcal=.53; tcal=float(mean(gcal,p)); hcal=float(dist(tcal)[0]); rcal=float(dist(tcal,r)[0])
    saturation=[]; compensation_error=0.; receptor_error=0.
    for T in [.5,5,36]:
        for kp in kpgrid:
            for kq in kqgrid:
                vh=hp/th*(kq+T*th)/(kp+T*hp)
                vr=rp/th*(kq+T*th)/(kp+T*rp)
                vhcal=hcal/tcal*(kq+T*tcal)/(kp+T*hcal)
                vrcal=rcal/tcal*(kq+T*tcal)/(kp+T*rcal)
                kr=vhcal/vrcal
                # Shared monotone GS=expit(log v), H normalized to 0.5 at cal.
                yh=vh/(vhcal+vh); yr=kr*vr/(vhcal+kr*vr)
                saturation.append([T,kp,kq,kr,np.max(np.abs(yh-yr))])
                # An unknown forward affinity restores exact equivalence with
                # ONE shared transformation at all G and T; no gain changes.
                kpR=kp/gain
                vr_aff=rp/th*(kq+T*th)/(kpR+T*rp)
                compensation_error=max(compensation_error,float(np.max(np.abs(vh-vr_aff)/(1+np.abs(vh)))))
                # Explicit four-state finite receptor, two independent binding
                # sites: (empty, P only, Q only, P+Q). Both occupancies match.
                wh=np.array([np.ones_like(th),T*hp/kp,T*th/kq,T*T*hp*th/(kp*kq)]).T
                wr=np.array([np.ones_like(th),T*rp/kpR,T*th/kq,T*T*rp*th/(kpR*kq)]).T
                wh/=wh.sum(axis=1)[:,None];wr/=wr.sum(axis=1)[:,None]
                receptor_error=max(receptor_error,float(np.max(np.abs(wh-wr))))
    save('separate_site_saturation_scenarios.csv','PII_uM,Kp_uM,Kq_uM,gain_R_over_H,max_absolute_GS_fraction_difference',saturation)
    saturation=np.array(saturation)
    # Compare enzyme activity half-range locations, without claiming assay transfer.
    def halfpoints(which):
        a0,b0,ar0,*_=enzyme_rates(0.,p)
        aF,bF,arF,*_=enzyme_rates(1e12,p)
        f0=[a0,b0,ar0][which]; fF=[aF,bF,arF][which]
        return np.exp(brentq(lambda z:(enzyme_rates(np.exp(z),p)[which]-(f0+fF)/2),-20,20))
    # These are macroscopic sequential binding constants, NOT effective C.
    k1=s*np.sqrt(c)*d[1]; k2=s/np.sqrt(c)*d[2]/d[1]
    report={'seed':SEED,'n_digitized_coordinates':len(x),'scope':'Finite buffered-ligand unsaturated catalytic reaction model; not a full GlnD/GlnE reconstitution.',
      'profile':cases,'fixed_C100_unconstrained_RMSE':float(np.sqrt(freesse/len(x))),'transported_C100':{'parameters':transportp.tolist(),'fitted_parameter_count':2,'transported_UR_half_range_mM':.08,'RMSE':float(np.sqrt(transportsse/len(x))),'success':transportok},'independent_binding':{'parameters':p.tolist(),'fitted_parameter_count':3,'effective_C_not_binding_ratio':c,'SSE':sse,'RMSE':float(np.sqrt(sse/len(x))),'LOO_RMSE':float(np.sqrt(np.mean((np.array(loo)[:,2]-y)**2))),'success':ok,'q_states':q.tolist(),'r':r,'gain_R_over_H':gain,'state_predictions_G10':{'H':dist(mean(10.,p)).tolist(),'R':dist(mean(10.,p),r).tolist()}},
      'finite_reaction_realization':{'binding_K1_mM':k1,'binding_K2_mM':k2,'binding_K1_over_K2':k1/k2,'alpha_H':alpha.tolist(),'beta_shared':beta.tolist(),'alpha_R':alphaR.tolist(),'catalytic_rate_scale':'arbitrary common positive rate scale','UR_specificity_fold_scenario':float((6.5/.82)/(2.7/2.3)),'UT_half_range_mM':halfpoints(0),'UR_half_range_mM':halfpoints(1),'UT_R_half_range_mM':halfpoints(2),'interpretation':'UR fold transported as a scenario; activity half-range points are not assumed identical to reported Ki/Kact.'},
      'checks':{'mean_equivalence_error':eqerr,'cascade_ratio_equivalence_error':verr,'finite_rate_mean_error_H':kineticerr,'finite_rate_mean_error_R':revkineticerr,'stationary_solve_error':staterr,'stationary_residual':residual,'positive_rate_check':bool(min(alpha.min(),beta.min(),alphaR.min())>=-1e-14)},
      'separate_site_saturation':[{'PII_uM':T,'scenario_min_max_difference':float(saturation[saturation[:,0]==T,4].min()),'scenario_max_max_difference':float(saturation[saturation[:,0]==T,4].max())} for T in [.5,5,36]],
      'unknown_affinity_equivalence':{'Kp_R_over_Kp_H':1/gain,'example_Kp_H_uM':.1,'example_Kp_R_uM':.1/gain,'saturated_ratio_error':compensation_error,'four_occupancy_class_error':receptor_error,'scope':'Buffered free PII input; rescale every unmodified-PII dissociation constant by the same factor. Reverse recognition weighted by modification count. Total-bound conservation not included.'},
      'strict_cross_assay_transfer':{'unliganded_specificity_ratio':ratio0,'q0':fixed_q0,'parameters':strictp.tolist(),'RMSE':float(np.sqrt(np.mean(strictfit.fun**2))),'success':bool(strictfit.success),'interpretation':'Transporting all three kinetic summaries across assays fails in this fixed-C minimal topology; this is not validation of a GlnD mechanism.'},
      'uncertainty':'No replicate likelihood or experimental confidence intervals; C bounds and cross-assay apparent-constant envelopes are scenarios.'}
    assert compensation_error<1e-12 and receptor_error<1e-12
    (OUT/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    # Four filled panels, vector text and curves, raster companion at 600 dpi.
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    fig,axs=plt.subplots(2,2,figsize=(7.3,5.7),layout='constrained')
    colors=['#3F7193','#B87973','#75947A','#AA9768']
    ax=axs[0,0]; ax.scatter(x,y,s=17,c='black',label='Digitized mean',zorder=4)
    xp=np.geomspace(.015,15,400)
    ax.plot(xp,3*mean(xp,p),c=colors[0],label='Finite homogeneous')
    ax.plot(xp,3*(r+(1-r)*(mean(xp,p)-r)/(1-r)),c=colors[1],ls='--',label='Finite refractory')
    ax.plot(xp,3*mean(xp,transportp),c='0.5',ls=':',lw=1,label='Two-summary transfer')
    ax.set(xscale='log',xlabel='Free glutamine (mM)',ylabel='UMP groups per trimer',ylim=(0,3.15));ax.legend(fontsize=7,loc='lower left')
    ax=axs[0,1]; ix=np.arange(4); h=dist(mean(10,p)); rr=dist(mean(10,p),r)
    ax.plot(ix,h,'-o',ms=5,color=colors[0],label='Homogeneous');ax.plot(ix,rr,'-s',ms=5,color=colors[1],label='Refractory')
    ax.set(xticks=ix,xlabel='UMP groups per trimer',ylabel='Fraction at 10 mM');ax.legend(fontsize=7)
    a0,b0,ar0,*_=enzyme_rates(0.,p);aF,bF,arF,*_=enzyme_rates(1e12,p)
    ax=axs[1,0]; ax.plot(xx,a/a0,c=colors[0],label='UT, homogeneous');ax.plot(xx,ar/ar0,c=colors[1],ls='--',label='UT, refractory');ax.plot(xx,b/bF,c=colors[2],label='UR, both / maximum')
    ax.set(xscale='log',xlim=(.002,100),xlabel='Free glutamine (mM)',ylabel='Relative catalytic capacity');ax.legend(fontsize=7)
    ax=axs[1,1]
    for i,T in enumerate([.5,5,36]):
        sub=saturation[saturation[:,0]==T]; z=sub[:,4].reshape(len(kpgrid),len(kqgrid))
        ax.plot(kpgrid,z[:,0],color=colors[i],label=f'{T:g} μM PII')
        ax.fill_between(kpgrid,z.min(axis=1),z.max(axis=1),color=colors[i],alpha=.18)
    ax.axhline(0,color='black',ls='--',lw=1,label='Affinity rescaled (all T)')
    ax.set(xscale='log',xlabel='Assumed forward half-saturation (μM)',ylabel='Maximum GS-fraction difference'); ax.legend(fontsize=7)
    for ax,l in zip(axs.flat,'ABCD'):ax.text(-.15,1.04,l,transform=ax.transAxes,fontweight='bold',fontsize=12)
    fig.savefig(FIG/'Fig7_finite_mechanisms.pdf');fig.savefig(FIG/'Fig7_finite_mechanisms.png',dpi=600);plt.close(fig)
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
