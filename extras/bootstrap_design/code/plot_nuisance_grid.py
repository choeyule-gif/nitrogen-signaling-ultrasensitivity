"""Display every audited null and small-effect nuisance cell (S9 Fig)."""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
R=Path(__file__).resolve().parents[1];O=R/'figures';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':8,'font.family':'DejaVu Sans','pdf.fonttype':42,'axes.titlesize':9})
fig,axs=plt.subplots(2,2,figsize=(7.5,6.3));fig.subplots_adjust(left=.12,right=.90,top=.92,bottom=.10,hspace=.48,wspace=.54)
settings=[('results',.17,'Shared-fit null',.05),('results_condition_affine',.17,'Condition-specific null',.05),('results_condition_affine',.19,r'Condition-specific: $\alpha_1=0.19$',1),('results_condition_affine',.25,r'Condition-specific: $\alpha_1=0.25$',1)]
pairs=[(c,s) for c in [0,.5,1] for s in [0,1,2]]
cmap=LinearSegmentedColormap.from_list('blue',['white','#0690F6'])
for letter,ax,(folder,alpha,title,upper) in zip('ABCD',axs.flat,settings):
 df=pd.read_csv(R/folder/'summary.csv');df=df[(df.stage=='evaluation')&(df.alpha==alpha)]
 values=np.array([[df[(df.kg==kg)&(df.calibration==c)&(df.signal==s)].rejection_rate.item() for kg in [7.8,11.7,15.6]] for c,s in pairs])
 im=ax.imshow(values,vmin=0,vmax=upper,cmap=cmap,aspect='auto')
 for (i,j),value in np.ndenumerate(values):ax.text(j,i,f'{value:.3f}',ha='center',va='center',fontsize=6.5,color='black')
 ax.set_xticks(range(3),['7.8','11.7','15.6']);ax.set_yticks(range(9),[f'{c:g}, {s:g}' for c,s in pairs])
 ax.set(xlabel=r'$K_G$ (mM)',ylabel='Error scales: concentration, signal',title=title)
 ax.text(-.30,1.10,letter,transform=ax.transAxes,fontsize=12,fontweight='bold')
 bar=fig.colorbar(im,ax=ax,fraction=.045,pad=.035);bar.set_label('Rejection probability',fontsize=7);bar.ax.tick_params(labelsize=7)
fig.savefig(O/'Fig_nuisance_grid.pdf');fig.savefig(O/'Fig_nuisance_grid.png',dpi=220)
