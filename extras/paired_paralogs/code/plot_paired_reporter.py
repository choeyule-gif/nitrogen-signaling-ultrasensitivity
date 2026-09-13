from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from native_trimers import data
from paired_reporter import free_k
R=Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':'sans-serif','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.8,'legend.frameon':False})
f,axs=plt.subplots(2,2,figsize=(10.4,7.6),layout='constrained')
d=data('WT');cols=['#0690F6','#FF7B02']
for j,name in enumerate(['GlnB','GlnK']):
 z=d[d[:,1]==j];axs[0,0].errorbar(z[:,0],z[:,2],yerr=z[:,3],fmt='o-',color=cols[j],ms=4,lw=1,label=name)
axs[0,0].set(xlabel='Time from nitrogen upshift (min)',ylabel='Paired-peptide modified fraction',ylim=(0,1.08));axs[0,0].legend()
x=np.linspace(.02,.98,301)
for rho,c in zip([1,3,9],['#2BD302','#0690F6','#FF7B02']):axs[0,1].plot(x,free_k(x,rho),color=c,label=fr'$\rho={rho}$')
axs[0,1].set(xlabel='Free GlnB modified fraction',ylabel='Predicted free GlnK fraction',xlim=(0,1),ylim=(0,1));axs[0,1].legend()
from reporter_intervals import calibration,sequestration
rho_interval=calibration((.48,.52),(.73,.77))
error=np.linspace(0,.15,151)
bounds=np.array([sequestration((.48,.52),(.35,.39),rho_interval,e) for e in error])
axs[1,0].fill_between(error,bounds[:,0],bounds[:,1],color='#0690F6',alpha=.18)
axs[1,0].plot(error,bounds,color='#0690F6',lw=1.5)
axs[1,0].set(xlabel='Allowed error in free-pool relation',ylabel='Synthetic bound-fraction interval',ylim=(0,1))
# Measured fraction pairs must be distinguished from estimated bound fractions.
for j in range(0,len(d),2):
 z=d[j:j+2]
 if len(z)==2 and z[0,0]==z[1,0]:axs[1,1].plot(z[0,2],z[1,2],'o',color='#0690F6',ms=5)
for s,c in zip([0,.3,.6],['#2BD302','#FF7B02','#0E56FF']):axs[1,1].plot(x,(1-s)*free_k(x,3),color=c,label=fr'$s={s}$')
axs[1,1].set(xlabel='Measured GlnB fraction',ylabel='Measured GlnK fraction',xlim=(0,1),ylim=(0,1));axs[1,1].legend(title=r'Scenario: $\rho=3$')
for a,lab,title in zip(axs.flat,'ABCD',['Independent in vivo measurements','Calibration of the free-pool relation','Synthetic assay uncertainty','Conditional sequestration curves']):a.set_title(title,loc='left',fontsize=11);a.text(-.13,1.06,lab,transform=a.transAxes,fontweight='bold',fontsize=14)
f.savefig(R/'results/paired_reporter.png',dpi=600)
f.savefig(R/'results/paired_reporter.pdf')
