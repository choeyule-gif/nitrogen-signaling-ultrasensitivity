from pathlib import Path
import sys,json,ast,io
import numpy as np
from scipy.optimize import brentq,minimize_scalar,least_squares
from scipy.special import logsumexp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pypdf import PdfReader,PdfWriter,Transformation
from reportlab.pdfgen.canvas import Canvas
ROOT=Path(__file__).resolve().parents[1];F=ROOT/'figures';D=ROOT/'data/audit';FD=ROOT/'figdata';sys.path.insert(0,str(ROOT/'code'))
from esbm.model import nH_local_tot
plt.rcParams.update({'font.size':10,'axes.labelsize':11,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False})
blue='#4d88cf';orange='#d87978';teal='#33a99e'
def save(fig,name):
 fig.savefig(F/(name+'.pdf'),bbox_inches='tight');fig.savefig(F/(name+'.png'),dpi=160,bbox_inches='tight');plt.close(fig)
def letter(ax,ch):ax.text(-.15,1.08,ch,transform=ax.transAxes,fontweight='bold',fontsize=13)
def csvread(n):return np.genfromtxt(FD/n,delimiter=',',names=True,comments='#',skip_header=1)
def assemble(names,out):
 W,H=640,290*((len(names)+1)//2);wr=PdfWriter();pg=wr.add_blank_page(W,H)
 for i,n in enumerate(names):
  p=PdfReader(F/'matlab'/n).pages[0];b=p.cropbox;scale=min(308/float(b.width),266/float(b.height));p.add_transformation(Transformation().translate(-float(b.left),-float(b.bottom)));pg.merge_transformed_page(p,Transformation().scale(scale).translate((i%2)*320+6,H-(i//2+1)*290+6))
 buf=io.BytesIO();cv=Canvas(buf,pagesize=(W,H));cv.setFont('Helvetica-Bold',15)
 for i in range(len(names)):cv.drawString((i%2)*320+7,H-(i//2)*290-14,chr(65+i))
 cv.save();buf.seek(0);pg.merge_page(PdfReader(buf).pages[0]);wr.write(F/(out+'.pdf'))
assemble(['FigS1A.pdf','FigS1B.pdf'],'FigS1_invariants')
# S2: stable weighted moments; refine every interior maximum over log input.
la=np.linspace(-3,3,61);lb=la.copy();zz=np.linspace(-16,16,401);ii=np.arange(4.)
def nuval(a,b,z):
 logs=np.log([1,3*a,3*b,1]);w=logs[:,None]+ii[:,None]*np.atleast_1d(z);p=np.exp(w-logsumexp(w,axis=0));mu=(ii[:,None]*p).sum(0);un=((3-ii)[:,None]*p).sum(0);var=(((ii[:,None]-mu)**2)*p).sum(0);return 3*var/(mu*un)
def maximum(a,b):
 vv=nuval(a,b,zz);best=max(1.,vv.max())
 for i in np.where((vv[1:-1]>vv[:-2])&(vv[1:-1]>vv[2:]))[0]+1:
  best=max(best,-minimize_scalar(lambda z:-nuval(a,b,z)[0],bounds=(zz[i-1],zz[i+1]),method='bounded').fun)
 return best
Z=np.array([[maximum(10**a,10**b) for a in la] for b in lb]);LA,LB=np.meshgrid(la,lb)
fig=plt.figure(figsize=(6.2,4.5));ax=fig.add_subplot(projection='3d',computed_zorder=False);ax.plot_surface(LA,LB,Z,cmap='viridis',edgecolor='.45',linewidth=.15,alpha=.75,zorder=1)
for k,col in [('fit',blue),('range',orange)]:
 tab=csvread(f'FigS2_locus_{k}.csv');zs=[maximum(10**a,10**b) for a,b in zip(tab['log10_a'],tab['log10_b'])];ax.plot(tab['log10_a'],tab['log10_b'],np.array(zs)+.025,color=col,lw=2,label=k+' constraint',zorder=5,path_effects=[__import__('matplotlib.patheffects',fromlist=['Stroke']).Stroke(linewidth=3.8,foreground='white'),__import__('matplotlib.patheffects',fromlist=['Normal']).Normal()])
ax.scatter([0],[0],[1],edgecolor='k',facecolor='white',s=35);ax.set(xlabel=r'$\log_{10}a$',ylabel=r'$\log_{10}b$',zlabel=r'$\sup_u\nu$',zlim=(1,3.05));ax.view_init(25,-42);ax.legend(loc='upper left',fontsize=8);save(fig,'FigS2_ladder')
np.savez(D/'ladder_surface_checked.npz',log10a=la,log10b=lb,nu_sup=Z)
# S3: retain author MATLAB surface after recalculating its full grid.
tab=csvread('FigS3_plateau_surface.csv');print('S3 columns',tab.dtype.names,flush=True)
X=tab[tab.dtype.names[0]];Y=tab[tab.dtype.names[1]];old=tab[tab.dtype.names[2]]
new=np.array([nH_local_tot(x*y/1.2,y/1.2) for x,y in zip(X,Y)]);err=float(np.max(abs(new-old)));print('S3 max error',err,flush=True)
import shutil
shutil.copy2(F/'matlab/FigS3.pdf',F/'FigS3_concentrations.pdf')
# S4: interaction parameter w, with condition-specific local factors.
fit=json.loads((ROOT/'data/analysis_summary.json').read_text())['parameters'];lo,hi,S,h=fit['Umin']/3,fit['Umax']/3,fit['S_mM'],fit['h']
def theta(g):return lo+(hi-lo)/(1+(g/S)**h)
def pmean(w,k,t):
 i=np.arange(13.);logs=np.array([np.log(float(__import__('math').comb(12,int(j)))) for j in i])+i*np.log(k*(1-t)**3/t)+.5*i*(i-1)*np.log(w);p=np.exp(logs-logsumexp(logs));mu=p@i;return mu/12,12*(p@(i-mu)**2)/(mu*(p@(12-i)))
ws=csvread('FigS4_ratio_vs_nu2.csv')['w'];conditions=[(.5,.600,5.23),(5,.530,6.46),(36,.225,4.48)];mat=[]
for w in ws:
 row=[w]
 for pii,mid,obs in conditions:
  def midfun(lk):
   k=np.exp(lk);a=pmean(w,k,hi)[0];b=pmean(w,k,lo)[0];return (pmean(w,k,theta(mid))[0]-(a+b)/2)
  grid=np.linspace(-40,10,301);vals=np.array([midfun(z) for z in grid]);ix=np.where(vals[:-1]*vals[1:]<0)[0];ix=[i for i in ix if abs(vals[i])+abs(vals[i+1])>1e-10][0];lk=brentq(midfun,grid[ix],grid[ix+1]);k=np.exp(lk);a=pmean(w,k,hi)[0];b=pmean(w,k,lo)[0]
  lg=[]
  for q in [.1,.9]:
   target=a+q*(b-a);t=brentq(lambda t:pmean(w,k,t)[0]-target,lo+1e-14,hi-1e-14);lg.append(np.log(S)+np.log((hi-t)/(t-lo))/h)
  coeff=np.log(81)/(lg[1]-lg[0]);row.extend([k,pmean(w,k,theta(mid))[1],coeff,obs/coeff])
 mat.append(row)
mat=np.array(mat);np.savetxt(D/'downstream_ladder_checked.csv',mat,delimiter=',',header='w,'+','.join(f'{v}_PII{pii}' for pii,_,_ in conditions for v in ['Khat','nu_at_midpoint','range_coefficient','ratio']),comments='')
fig,axs=plt.subplots(1,2,figsize=(8,3.1),layout='constrained')
for j,((pii,_,_),col) in enumerate(zip(conditions,[blue,orange,teal])):
 axs[0].plot(mat[:,0],mat[:,4+4*j],color=col,label=f'{pii:g} μM');axs[1].plot(mat[:,0],mat[:,2+4*j],color=col,label=f'{pii:g} μM')
axs[0].axhline(1,color='.5',ls='--');axs[0].set(xlabel='Interaction parameter w',ylabel='Reported / calculated');axs[1].set(xlabel='Interaction parameter w',ylabel=r'Local $\nu_2$ at calibrated midpoint');axs[0].legend(fontsize=8)
for a,ch in zip(axs,'AB'):letter(a,ch)
save(fig,'FigS4_downstream')
# S5: corrected original A, retained original B, and isolated failures C and D.
source=json.loads((D/'counterexamples_parameters.json').read_text());replacement=json.loads((D/'replacement_condition_i.json').read_text())
from math import comb
cases=[dict(n=3,A=1.,c=[1,3,3,1],Bhat=[4.5,3,1.5],Ghat=[.8,1.6,2.4],TE=1.,TS=4.,TL=2.),source['condition_i'],source['condition_ii'][-1],replacement]
def params(d):
 n=d['n'];A=d['A'];bh=np.array(d['Bhat'][:n]);
 if 'c' in d:c=np.array(d['c']);gh=np.array(d['Ghat'])
 else:c=np.r_[1,np.cumprod(d['K'])];gh=d['rho']*bh/(A*np.array(d['K']))
 return n,A,bh,gh,c

def state(d,l):
 n,A,bh,gh,c=params(d);rr=c*l**(-np.arange(n+1));M=rr.sum();Q=bh@rr[:-1]+A*l*(gh@rr[1:]);aa=(1+A*l)*Q;bb=(1+A*l)*M+(d['TS']-d['TE'])*Q;e=2*d['TE']*M/(bb+np.sqrt(bb*bb+4*aa*d['TE']*M));s0=d['TS']/(M+e*Q);return l+A*e*l+A*e*s0*l*(gh@rr[1:])-d['TL']
fig,axs=plt.subplots(2,2,figsize=(8,6),layout='constrained');rootsall=[]
for j,(d,ax) in enumerate(zip(cases,axs.flat)):
 limits=[(.01,100),(3e-4,.5),(4e-4,6),(1e-4,600)][j];ll=np.logspace(*np.log10(limits),1800);rv=np.array([state(d,l) for l in ll]);idx=np.where(rv[:-1]*rv[1:]<0)[0];roots=[brentq(lambda l:state(d,l),ll[i],ll[i+1],xtol=1e-15) for i in idx];rootsall.append(roots)
 norm=1 if j<2 else d['TL'];ax.semilogx(ll,rv/norm,color=blue);ax.axhline(0,color='.5',ls='--');ax.plot(roots,np.zeros(len(roots)),'o',color=orange);ax.set(xlabel='Free effector l',ylabel='R(l)' if j<2 else r'$R(l)/T_L$');letter(ax,'ABCD'[j])
 if j==0:ax.set_ylim(-1.5,6)
 if j==1:ax.set_ylim(-.3,.35)
 if j==2:ax.set_yscale('symlog',linthresh=.05);ax.set_yticks([-1,-.1,0,.1,1])
 if j==3:ax.set_yscale('symlog',linthresh=.1);ax.set_yticks([-1,-.1,0,.1,1])
save(fig,'FigS5_uniqueness');(D/'figure_S5_roots.json').write_text(json.dumps(rootsall,indent=2))
# S6: recompute Re F from the direct resolvent; use the same threshold partition for ALL n.
p=ROOT/'code/ref_model.py';tree=ast.parse(p.read_text());ns={'np':np};exec(compile(ast.Module(body=[x for x in tree.body if isinstance(x,ast.FunctionDef)],type_ignores=[]),str(p),'exec'),ns)
fig=plt.figure(figsize=(6.5,4.5));ax=fig.add_subplot(projection='3d');check=[]
for n in [6,4,3,2,1]:
 L,q,pi,phi,m=ns['chain'](n);om=np.logspace(-3,3,300);thr=ns['circ_threshold'](L,q,pi,m);Fvals=np.array([(phi@(pi*(-L@np.linalg.solve(1j*w*np.eye(m)-L,phi)))).real for w in om]);cut=np.searchsorted(om,thr)
 ax.plot(np.log10(om[:cut+1]),np.full(cut+1,n),Fvals[:cut+1],color=orange,lw=2);ax.plot(np.log10(om[cut:]),np.full(len(om)-cut,n),Fvals[cut:],color=blue,lw=2);ax.scatter([-3],[n],[1-(pi@phi)**2],facecolor='white',edgecolor=blue,s=25);ax.scatter([np.log10(thr)],[n],[0],color='.35',marker='|',s=80)
 check.append(dict(n=n,min_ReF=float(Fvals.min()),F0=float(1-(pi@phi)**2)))
ax.set(xlabel=r'$\log_{10}\omega$',ylabel='n',zlabel=r'$\mathrm{Re}\,F$',zlim=(0,1.05),yticks=[1,2,3,4,6]);ax.view_init(27,-51);save(fig,'FigS6_stability')
# S7: restore direct-route and concentration-discrimination calculations.
def baseline(g,k):
 t=theta(g);return k*(1-t)**3/(t+k*(1-t)**3)
def limits(k):return baseline(0,k),baseline(1e15,k)
def calib(mid):return np.exp(brentq(lambda z:baseline(mid,np.exp(z))-sum(limits(np.exp(z)))/2,-30,30))
def score(fn):
 a=fn(0);b=fn(1e18);lg=[brentq(lambda z:fn(np.exp(z))-(a+q*(b-a)),-35,40) for q in [.1,.9]];return np.log(81)/(lg[1]-lg[0])
k=calib(.53);base=score(lambda g:baseline(g,k));kgs=np.logspace(-2,2,81)*S;gain=[]
for kg in kgs:
 def fn(g):
  phi=(1+3.7*g/kg)/(1+g/kg);return baseline(g,k*phi)
 gain.append(score(fn)/base)
fig,axs=plt.subplots(1,3,figsize=(10,3),layout='constrained');ks=np.logspace(-8,8,140);axs[0].semilogx(ks,[score(lambda g:baseline(g,k))/h for k in ks],color=blue);axs[0].set(xlabel=r'$\kappa$',ylabel='Composite range coefficient / h')
axs[1].semilogx(kgs/S,gain,color=blue);axs[1].axhline(1,ls='--',color='.5');axs[1].set(xlabel=r'$K_G/S$',ylabel='Direct-route gain (R = 3.7)')
tab=csvread('Fig5BC_titrations.csv')
for lab,col,sty in [('1x1x',blue,'-'),('2x1x',orange,'--'),('1x4x',teal,':'),('2x4x','#d69b41','-.')]:
 axs[2].semilogx(tab['d_TL_'+lab],tab['d_theta_'+lab],color=col,ls=sty,label=lab)
axs[2].set(xlabel=r'Total effector $T_L$',ylabel=r'Total modified fraction $\theta$');axs[2].legend(fontsize=7)
for a,ch in zip(axs,'ABC'):letter(a,ch)
save(fig,'FigS7_alternatives');np.savetxt(D/'direct_route_checked.csv',np.c_[kgs,gain],delimiter=',',header='KG,gain',comments='')
(D/'figure_checks.json').write_text(json.dumps(dict(S3_max_error=err,S3_min=float(new.min()),S3_max=float(new.max()),S3_interior_min=float(new[X<=.01].min()),S2_min=float(Z.min()),S2_max=float(Z.max()),S7_direct_gain=[min(gain),max(gain)],S6=check),indent=2))
print('Seven supplementary figures prepared.',flush=True)
