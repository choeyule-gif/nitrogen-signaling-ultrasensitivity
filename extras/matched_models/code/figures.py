"""Vector scientific figures; all numbers loaded from reproducible analyses."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze import ROOT,OUT,distribution,hill
F=ROOT/'figures';F.mkdir(exist_ok=True)
r=json.loads((OUT/'summary.json').read_text());p=np.array(r['reference']['parameters'])
names=list(r['matched']);colors=['#427BA8','#8764A1','#42978B','#D28B50']
labels=['Independent','End-state enriched','Middle-state enriched','Refractory mixture']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'legend.fontsize':7,'xtick.labelsize':7,'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,'lines.linewidth':1.7})
def panel(ax,letter,title):
    ax.set_title(title,loc='left',pad=10)
    ax.text(-.14,1.08,letter,transform=ax.transAxes,weight='bold',fontsize=12,va='bottom')
def save(fig,name):
    fig.savefig(F/(name+'.pdf'),bbox_inches='tight')
    fig.savefig(F/(name+'.png'),dpi=220,bbox_inches='tight')
    fig.savefig(F/(name+'.tif'),dpi=600,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'})
    plt.close(fig)

fig,axs=plt.subplots(2,2,figsize=(7.35,5.65),layout='constrained')
a=axs[0,0];panel(a,'A','Finite-switch fits to the same data')
d=np.loadtxt(ROOT/'data/jiangninfa2011_fig2B_digitised.csv',delimiter=',',skiprows=1)
c=np.loadtxt(OUT/'finite_switch_curves.csv',delimiter=',',skiprows=1)
for j,ix in enumerate([2,0,1]):a.semilogx(c[:,0],c[:,j+1],color=colors[ix],label=labels[ix])
a.scatter(*d.T,s=15,color='#252525',zorder=5,label='Digitized coordinates')
a.set(xlim=(.015,15),ylim=(0,3.1),xlabel='Glutamine (mM)',ylabel='Mean UMP groups per trimer')
a.legend(frameon=False,loc='lower left',fontsize=6.7)
b=axs[0,1];g=r['design']['best_G_mM'];panel(b,'B',f'Mean-matched states at {g:.2f} mM')
for j,(name,m) in enumerate(r['matched'].items()):
    prob=distribution(g,p,m['J'],m['refractory'])
    b.plot(np.arange(4),prob,marker=['o','s','^','D'][j],ms=5,color=colors[j],label=labels[j])
b.set(xticks=np.arange(4),xticklabels=['0','1','2','3'],xlabel='UMP groups per trimer',ylabel='Trimer fraction',ylim=(0,.65))
b.legend(frameon=False,fontsize=6.5)
cax=axs[1,0];panel(cax,'C','High-input state predictions')
for j,name in enumerate(names):
    q=np.array(r['state_prediction_percentiles'][name+'_10']);v=q[1,3]
    cax.errorbar(j,v,yerr=[[v-q[0,3]],[q[2,3]-v]],fmt=['o','s','^','D'][j],ms=7,color=colors[j],ecolor=colors[j],elinewidth=2,capsize=4)
cax.set(xticks=range(4),xticklabels=['Independent','End-state\nenriched','Middle-state\nenriched','Refractory'],ylabel=r'PII-UMP$_3$ fraction at 10 mM',ylim=(-.006,.23))
cax.text(.03,.97,'Points: medians; intervals: coordinate perturbations',transform=cax.transAxes,va='top',fontsize=6.5)
dax=axs[1,1];panel(dax,'D','Input selection by state separation')
dat=np.loadtxt(OUT/'state_design.csv',delimiter=',',skiprows=1)
for k in range(2,dat.shape[1]):dax.semilogx(dat[:,0],dat[:,k],color='#CDCFD2',linewidth=.7,zorder=1)
dax.semilogx(dat[:,0],dat[:,1],color='#2F5268',linewidth=2,label='Minimum across six pairs')
dax.axvline(g,color='#555555',ls=':',lw=1)
dax.scatter([g],[r['design']['minimum_pair_TV']],color='#2F5268',s=22,zorder=4)
dax.set(xlabel='Glutamine (mM)',ylabel='Total-variation distance',xlim=(.02,10),ylim=(0,float(dat[:,2:].max())*1.12))
dax.legend(frameon=False,loc='upper left',fontsize=6.5)
save(fig,'Fig5_matched_models')

fig,axs=plt.subplots(2,2,figsize=(7.35,6.05),layout='constrained')
a=axs[0,0];panel(a,'A','Equal curves after gain calibration')
ca=np.loadtxt(OUT/'cascade_independent_baseline.csv',delimiter=',',skiprows=1)
cr=np.loadtxt(OUT/'cascade_refractory_baseline.csv',delimiter=',',skiprows=1)
a.semilogx(ca[:,0],ca[:,2],color=colors[0],label='Independent')
ii=np.flatnonzero((cr[:,0]>.035)&(cr[:,0]<8))[::12]
a.semilogx(cr[ii,0],cr[ii,2],ls='none',marker='o',ms=3,mfc='none',mec=colors[3],label='Refractory; rescaled '+r'$\kappa$')
a.scatter([.53],[.5],marker='x',s=35,color='black',label='Calibrated midpoint')
a.set(xlim=(.03,10),ylim=(-.03,1.03),xlabel='Glutamine (mM)',ylabel='Normalized GS response')
a.legend(frameon=False,fontsize=6.8,loc='upper left')
b=axs[0,1];panel(b,'B','Baseline depends on the PII ladder')
obs=np.array(r['cascade'][0]['reported'])
for j,name in enumerate(names[:3]):
    item=next(t for t in r['cascade'] if t['upstream']==name and t['kind']=='baseline')
    b.plot(np.arange(3),item['predicted'],marker='o',ms=4,color=colors[j],label=['Independent','Endpoint-rich','Intermediate-rich'][j])
b.plot(np.arange(3),obs,'ks',ms=5,label='Reported',ls='none')
b.set(xticks=range(3),xticklabels=['0.5','5','36'],xlabel=r'PII ($\mu$M)',ylabel='Response-range coefficient',ylim=(2,8))
b.legend(frameon=False,fontsize=6.4,loc='upper left',ncol=2,columnspacing=.7,handlelength=1.5)
c=axs[1,0];panel(c,'C','Joint downstream refits')
arr=[];yl=[]
kn={'baseline':'Baseline (0)','direct':'Direct route (2)','power':'Transfer power (1)','gs_ladder':'GS ladder (1)'}
for nm,tag in [('independent','Ind.'),('central_enriched','Middle')]:
    for kd in kn:
        t=next(t for t in r['cascade'] if t['upstream']==nm and t['kind']==kd)
        arr.append(np.array(t['predicted'])/obs);yl.append(tag+': '+kn[kd])
im=c.imshow(arr,vmin=.5,vmax=1.5,cmap='BrBG',aspect='auto')
for i in range(8):
    for j in range(3):c.text(j,i,f'{arr[i][j]:.2f}',ha='center',va='center',fontsize=7,color='white' if arr[i][j]<.66 else '#222222')
c.set(xticks=range(3),xticklabels=['0.5','5','36'],yticks=range(8),yticklabels=yl,xlabel=r'PII ($\mu$M)')
c.tick_params(axis='y',labelsize=6.5)
cb=fig.colorbar(im,ax=c,fraction=.046,pad=.025);cb.set_label('Predicted / reported',fontsize=7)
d=axs[1,1];panel(d,'D','Joint upstream and downstream fits')
for j,t in enumerate(r['joint_shared_ladder']):
    lab={'baseline':'Ladder only','direct':'+ direct route','power':'+ transfer power','gs_ladder':'+ GS ladder'}[t['kind']]
    d.plot(range(3),t['predicted'],marker=['o','s','^','D'][j],ms=4,color=colors[j],label=lab)
d.plot(range(3),obs,'kx',ms=6,mew=1.5,label='Reported')
d.set(xticks=range(3),xticklabels=['0.5','5','36'],xlabel=r'PII ($\mu$M)',ylabel='Response-range coefficient',ylim=(3.7,7.1))
d.legend(frameon=False,fontsize=6.6,loc='lower center')
save(fig,'Fig6_cascade_comparison')
print(F)
