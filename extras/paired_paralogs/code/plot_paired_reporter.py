from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from native_trimers import data
from paired_reporter import free_k
R=Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.8,'legend.frameon':False})
f,axs=plt.subplots(2,2,figsize=(7.5,5.8),layout='constrained')
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
# State-count bounds need no stationarity or independent modification sites.
with (R/'results/state_bounds.csv').open() as fh: states=list(csv.DictReader(fh))
for estimator,col,label in [('paired_peptide_fraction','#0690F6','Paired peptides'),('modified_over_total','#FF7B02','Separate total')]:
 z=sorted([v for v in states if v['strain']=='WT' and v['protein']=='GlnK' and v['estimator']==estimator],key=lambda v:float(v['time_min']))
 t=[float(v['time_min']) for v in z];upper=[float(v['p0_hi']) for v in z]
 axs[1,1].plot(t,upper,'o-',color=col,lw=1.3,ms=3,label=label+' + uncertainty')
 if estimator=='paired_peptide_fraction':
  axs[1,1].plot(t,[float(v['point_p0_max']) for v in z],'--',color=col,lw=1,label='Paired-peptide point ceiling')
axs[1,1].set(xlabel='Time from nitrogen upshift (min)',ylabel='Upper limit on bound GlnK fraction',ylim=(0,1.06));axs[1,1].legend(fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.27),borderaxespad=0)
for a,lab,title in zip(axs.flat,'ABCD',['Independent in vivo measurements','Calibration of the free-pool relation','Synthetic assay uncertainty','Observed modification limits binding']):a.set_title(title,loc='left',fontsize=9);a.text(-.13,1.06,lab,transform=a.transAxes,fontweight='bold',fontsize=11)
f.savefig(R/'results/paired_reporter.png',dpi=600)
f.savefig(R/'results/paired_reporter.pdf')
