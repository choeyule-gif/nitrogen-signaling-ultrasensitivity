"""Review figure for the new independently calibrated assay benchmark."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
BLUE='#0690F6';ORANGE='#FF7B02';GREEN='#2BD302';GRAY='#797979'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,
 'axes.spines.right':False,'axes.linewidth':.65,'pdf.fonttype':42,'lines.linewidth':1.5,
 'legend.fontsize':7,'legend.frameon':True,'legend.facecolor':'white','legend.edgecolor':'white',
 'legend.framealpha':.95,'axes.titlesize':8.5,'xtick.labelsize':7,'ytick.labelsize':7})
common=pd.read_csv(ROOT/'results/summary.csv')
affine=pd.read_csv(ROOT/'results_condition_affine/summary.csv')
fig,axes=plt.subplots(2,2,figsize=(7.5,6.1))
fig.subplots_adjust(left=.10,right=.98,bottom=.10,top=.92,wspace=.37,hspace=.55)
for letter,ax in zip('ABCD',axes.flat):
 ax.text(-.18,1.13,letter,transform=ax.transAxes,fontweight='bold',fontsize=12)
a=axes[0,0]
for data,color,label,shift in [(common,ORANGE,'Shared background/gain',-.08),(affine,BLUE,'Condition-specific fit',.08)]:
 q=data[(data.stage=='evaluation')&(data.alpha==.17)]
 for j,k in enumerate([7.8,11.7,15.6]):
  sub=q[q.kg==k]
  jitter=np.linspace(-.055,.055,len(sub))
  a.scatter(np.ones(len(sub))*j+shift+jitter,sub.rejection_rate,s=16,color=color,
            alpha=.8,label=label if j==0 else None)
a.axhline(.05,color=GRAY,ls='--',lw=.9)
a.set(xticks=[0,1,2],xticklabels=['7.8','11.7','15.6'],ylim=(-.003,.065),
      xlabel=r'$K_G$ (mM)',ylabel='Null rejection rate',title='Independent audit of null rejection')
a.legend(loc='upper left')
a=axes[0,1]
for kg,color in [(7.8,BLUE),(11.7,ORANGE),(15.6,GREEN)]:
 d=affine[(affine.stage=='evaluation')&(affine.kg==kg)]
 s=d.groupby('alpha').rejection_rate.agg(['min','max','mean'])
 a.fill_between(s.index,s['min'],s['max'],color=color,alpha=.14)
 a.plot(s.index,s['mean'],'o-',color=color,ms=3,label=fr'$K_G={kg:g}$ mM')
a.axhline(.05,color=GRAY,ls='--',lw=.9)
a.set(xlabel=r'Regulatory coefficient $\alpha_1$',ylabel='Reference rejection rate',
      ylim=(-.03,1.04),title='Power depends on nuisance values',xticks=[.17,.19,.25,.34])
a.legend(loc='lower right')
a=axes[1,0]
for data,color,label in [(common,ORANGE,'Shared fit'),(affine,BLUE,'Condition-specific fit')]:
 d=data[(data.stage=='signal_bias')&(data.kg==15.6)&(data.differential_gain==0)].sort_values('differential_background')
 a.plot(d.differential_background,d.rejection_rate,'o-',color=color,ms=3,label=label)
 a.fill_between(d.differential_background,d.ci_low,d.ci_high,color=color,alpha=.15)
a.axhline(.05,color=GRAY,ls='--',lw=.9)
a.set(xlabel='Background shift at nonzero glutamine',ylabel='Rejection with unchanged mechanism',
      ylim=(-.03,1.04),title='Fit condition-specific signal offsets',xticks=[-.02,-.01,0,.01,.02])
a.legend(loc='upper right')
a=axes[1,1]
for data,color,label in [(common,ORANGE,'Shared fit'),(affine,BLUE,'Condition-specific fit')]:
 d=data[(data.stage=='mixture')&(data.kg==15.6)].sort_values('mixture')
 a.plot(d.mixture,d.rejection_rate,'o-',color=color,ms=3,label=label)
 a.fill_between(d.mixture,d.ci_low,d.ci_high,color=color,alpha=.15)
a.axhline(.05,color=GRAY,ls='--',lw=.9)
a.set(xlabel='Hidden affinity mixture strength',ylabel='Single-class reference rejection rate',
      ylim=(-.03,1.04),title='Mixtures still alter the interpretation')
a.legend(loc='lower right')
fig.savefig(OUT/'Fig_bootstrap_controls.pdf')
fig.savefig(OUT/'Fig_bootstrap_controls.png',dpi=220)
plt.close(fig)
print(OUT)
