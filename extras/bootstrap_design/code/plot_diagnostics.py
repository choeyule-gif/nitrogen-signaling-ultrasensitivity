"""Retained raw-assay and GlnD-profile diagnostics for S8 Fig."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bootstrap_binding import occupancy

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'figures'
BLUE='#0690F6';ORANGE='#FF7B02';GREEN='#2BD302'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,
 'axes.spines.right':False,'axes.linewidth':.65,'pdf.fonttype':42,'lines.linewidth':1.5,
 'legend.fontsize':7,'legend.frameon':True,'legend.facecolor':'white','legend.edgecolor':'white',
 'legend.framealpha':.95,'axes.titlesize':8.5,'xtick.labelsize':7,'ytick.labelsize':7})
fig,axs=plt.subplots(1,2,figsize=(7.5,3.3))
fig.subplots_adjust(left=.09,right=.90,bottom=.20,top=.86,wspace=.58)
for l,a in zip('AB',axs):a.text(-.18,1.12,l,transform=a.transAxes,fontweight='bold',fontsize=12)
ex=json.loads((ROOT/'results_condition_affine/example.json').read_text())
p=np.array(ex['doses_for_fit']);obs=np.array(ex['observations']);z=np.array(ex['fitted_parameters'])
x=np.r_[0,np.geomspace(.001,5,180)]
for j,color in enumerate([BLUE,ORANGE,GREEN]):
 axs[0].scatter(np.repeat(p[j],2),obs[j].ravel(),s=10,facecolors='none',edgecolors=color,lw=.6)
 axs[0].plot(x,z[3+j]+z[6+j]*occupancy(x,np.exp(z[j]),ex['enzyme_measured'])[0],color=color,label=f'{[0,15.6,31.2][j]:g} mM Gln')
axs[0].set(xscale='symlog',xlabel=r'Total PII ($\mu$M)',ylabel='Normalized binding signal',title='Raw observations to halfpoints')
axs[0].set_xscale('symlog',linthresh=.005)
axs[0].set_xticks([0,.01,.1,1,5],['0','.01','.1','1','5']);axs[0].legend(loc='lower right')
pr=pd.read_csv(ROOT.parent/'observation_design/results/glnd_constraint_profile.csv')
a=axs[1];a.plot(pr.K1_mM,pr.RMSE,'o-',color=BLUE,ms=3,label='Titration RMSE')
a.set(xscale='log',xlabel=r'Fixed GlnD $K_1$ (mM)',ylabel='RMSE (UMP/trimer)',title='Fit quality and occupancy contrast')
a.yaxis.label.set_color(BLUE);a.tick_params(axis='y',colors=BLUE)
b=a.twinx();b.spines['right'].set_visible(True)
b.plot(pr.K1_mM,pr.max_p1_separation,'s--',color=ORANGE,ms=3,label='Occupancy contrast')
b.set(ylabel='Singly liganded contrast',ylim=(0,.5));b.yaxis.label.set_color(ORANGE);b.tick_params(axis='y',colors=ORANGE)
a.legend(loc='upper right',bbox_to_anchor=(1,.97));b.legend(loc='upper right',bbox_to_anchor=(1,.81))
fig.savefig(OUT/'Fig_raw_and_profile.pdf');fig.savefig(OUT/'Fig_raw_and_profile.png',dpi=220)
plt.close(fig)
