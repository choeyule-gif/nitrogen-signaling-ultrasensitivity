"""Protein-conserving, buffered-glutamine 120-species cascade.
Concentrations in micromolar, arbitrary common time units. This is a specified
reaction topology, not a fitted physiological model. Regulatory exchange occurs
only on catalyst-free enzyme. Nucleotide donors/products are buffered.
"""
import numpy as np
from scipy.special import gammaln, logsumexp
from scipy.optimize import least_squares

K=52.044461275635986
A=np.array([.9988956385628978,.0012739724723029192,.0012739724723029192])
B=np.array([.0011043614371022414,.0011043614371022414,.007457183145744586])

def ladder(n,logr,J=0):
    i=np.arange(n+1)
    v=gammaln(n+1)-gammaln(i+1)-gammaln(n-i+1)+i*logr+J*(i*i-n*i)
    return np.exp(v-logsumexp(v))

class Cascade:
    def __init__(self,G=530.,lam=1.,sigma=1.,direct=1.,J=0.,Pt=5.,St=10.,Dt=.05,Et=.025,recognition=1.,alpha2=4.,productive=0.):
        self.G,self.Pt,self.St,self.Dt,self.Et=G,Pt,St,Dt,Et
        self.J,self.direct=J,direct
        self.K1,self.K2=K*.5*lam,2*K/lam
        self.a=A.copy();self.b=B.copy();self.a[1]*=lam;self.b[1]*=lam
        self.KG,self.KP,self.KU=15600.,.18,.35
        self.alpha1,self.alpha2=.17*sigma,alpha2
        self.recognition=recognition
        self.at=np.array([0,1.4877648555011858*direct,1.4877648555011858,0,0,0,1.4877648555011858*10*sigma,0,0,0])
        self.at[7:]=productive*1.4877648555011858
        self.ar=np.array([0,0,0,1,1,1,0,0,0,0.])
        # Regulatory states: E, EG, EP0, EU1, EU2, EU3, EGP0, EP0U1..3.
        self.pc=[[],[],[0],[1],[2],[3],[0],[0,1],[0,2],[0,3]]
        self.gc=np.array([0,1,0,0,0,0,1,0,0,0])
        self.names=[];self.cons=[];self.pm=[];self.sm=[];self.gl=[]
        def sp(name,c,p=0,s=0,g=0):
            z=len(self.names);self.names.append(name);self.cons.append(c);self.pm.append(p);self.sm.append(s);self.gl.append(g);return z
        self.P=np.array([sp('P'+str(i),[1,0,0,0],i) for i in range(4)])
        self.S=np.array([sp('S'+str(i),[0,1,0,0],s=i) for i in range(13)])
        self.D=np.array([sp('D'+str(j),[0,0,1,0],g=j) for j in range(3)])
        self.Dc=[]
        for j in range(3):
            for sign,inds in [(1,range(3)),(-1,range(1,4))]:
                for i in inds:
                    c=sp(f'DC_{j}_{i}_{sign}',[1,0,1,0],i,g=j)
                    mult=3-i if sign==1 else i
                    spec=(self.a[j] if sign==1 else self.b[j])*mult
                    h=(np.array([.2,2,10]) if sign==1 else np.array([5,.3,.1]))[j]*mult
                    self.Dc.append((c,j,i,sign,spec,h))
        self.E=np.array([sp('E'+str(j),[len(self.pc[j]),0,0,1],sum(self.pc[j]),g=self.gc[j]) for j in range(10)])
        self.Ec=[]
        for j in [1,2,6,3,4,5]+([7,8,9] if productive>0 else []):
            sign=1 if j in [1,2,6,7,8,9] else -1
            for i in (range(12) if sign==1 else range(1,13)):
                c=sp(f'EC_{j}_{i}',[len(self.pc[j]),1,0,1],sum(self.pc[j]),i,self.gc[j])
                mult=(12-i)*np.exp(J*(2*i+1-12)) if sign==1 else i
                spec=(self.at[j] if sign==1 else self.ar[j])*mult
                # direct=0: EG catalytic edges absent; species slots retained with zero abundance.
                h=.1*mult if spec>0 else 0.
                self.Ec.append((c,j,i,sign,spec,h))
        self.cons=np.array(self.cons).T;self.pm=np.array(self.pm);self.sm=np.array(self.sm);self.gl=np.array(self.gl)
        self.cols=[];self.rs=[];self.ks=[]
        def reaction(react,prod,k):
            col=np.zeros(len(self.names))
            for z in react:col[z]-=1
            for z in prod:col[z]+=1
            self.cols.append(col);self.rs.append(react);self.ks.append(k)
        def bind(x,y,lig,ratio):
            reaction([x]+([] if lig is None else [lig]),[y],ratio)
            reaction([y],[x]+([] if lig is None else [lig]),1.)
        bind(self.D[0],self.D[1],None,G/self.K1);bind(self.D[1],self.D[2],None,G/self.K2)
        e=self.E
        bind(e[0],e[1],None,G/self.KG);bind(e[0],e[2],self.P[0],1/self.KP)
        bind(e[1],e[6],self.P[0],1/(self.alpha1*self.KP));bind(e[2],e[6],None,G/(self.alpha1*self.KG))
        for i in range(1,4):
            bind(e[0],e[2+i],self.P[i],(i/3)**recognition/self.KU)
            bind(e[2],e[6+i],self.P[i],(i/3)**recognition/(self.alpha2*self.KU))
            bind(e[2+i],e[6+i],self.P[0],1/(self.alpha2*self.KP))
        for cc,enz,tar in [(self.Dc,self.D,self.P),(self.Ec,self.E,self.S)]:
            for c,j,i,sign,spec,h in cc:
                if spec==0:continue
                cat=spec/h
                reaction([enz[j],tar[i]],[c],2*spec)
                reaction([c],[enz[j],tar[i]],cat)
                reaction([c],[enz[j],tar[i+sign]],cat)
        self.N=np.array(self.cols).T;self.ks=np.array(self.ks)
        self.r1=np.array([r[0] for r in self.rs]);self.r2=np.array([r[1] if len(r)==2 else len(self.names) for r in self.rs])
        assert len(self.names)==(156 if productive>0 else 120)
        assert np.max(np.abs(self.cons@self.N))==0
    def rhs(self,t,y):
        yy=np.append(y,1.)
        return self.N@(self.ks*yy[self.r1]*yy[self.r2])
    def construct(self,T,V):
        G=self.G
        wd=np.array([1,G/self.K1,G*G/(self.K1*self.K2)])
        pi=ladder(3,np.log(wd@self.a/(wd@self.b)))
        p=T*pi;p0=p[0]/self.KP;u=p[1:]*(np.arange(1,4)/3)**self.recognition/self.KU;g=G/self.KG
        we=np.array([1,g,p0,*u,p0*g/self.alpha1,*(p0*u/self.alpha2)])
        r=(we@self.at)/(we@self.ar)
        si=ladder(12,np.log(r),self.J);s=V*si
        denD=wd.sum()+sum(wd[j]*h*p[i] for c,j,i,sg,a,h in self.Dc)
        denE=we.sum()+sum(we[j]*h*s[i] for c,j,i,sg,a,h in self.Ec)
        d=self.Dt*wd/denD;e=self.Et*we/denE
        y=np.zeros(len(self.names));y[self.P]=p;y[self.S]=s;y[self.D]=d;y[self.E]=e
        for c,j,i,sg,a,h in self.Dc:y[c]=d[j]*h*p[i]
        for c,j,i,sg,a,h in self.Ec:y[c]=e[j]*h*s[i]
        return y
    def steady(self,start=(.8,.8)):
        totals=np.array([self.Pt,self.St,self.Dt,self.Et])
        def fun(z):return (self.cons@self.construct(*(np.exp(z)*totals[:2])))[:2]/totals[:2]-1
        fit=least_squares(fun,np.log(start),xtol=1e-14,gtol=1e-14,ftol=1e-14)
        y=self.construct(*(np.exp(fit.x)*totals[:2]))
        if np.max(np.abs(fun(fit.x)))>1e-8:raise RuntimeError(str((self.G,self.Dt,self.Et,fit.x,fun(fit.x),fit.message)))
        return y
    def obs(self,y):
        pf=y[self.P]/y[self.P].sum();sf=y[self.S]/y[self.S].sum()
        pt=np.zeros(4)
        for i in range(4):pt[i]+=y[self.P[i]]
        for c,j,i,sg,a,h in self.Dc:pt[i]+=y[c]
        for j,idx in enumerate(self.E):
            for i in self.pc[j]:pt[i]+=y[idx]
        for c,j,k,sg,a,h in self.Ec:
            for i in self.pc[j]:pt[i]+=y[c]
        return dict(pfree=pf,ptotal=pt/self.Pt,GSfree=sf@np.arange(13)/12, GStotal=y@self.sm/(12*self.St),Pfree=y[self.P].sum(),boundG=y@self.gl)
