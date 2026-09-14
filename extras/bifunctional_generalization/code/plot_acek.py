"""Review figure: source-summary calibration, finite kinetics and enclosure."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from acek_kinetic_ambiguity import law,network,stationary

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
BLUE='#0690F6';ORANGE='#FF7B02';GREEN='#2BD302';GRAY='#797979'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,
 'axes.spines.right':False,'axes.linewidth':.65,'pdf.fonttype':42,'lines.linewidth':1.5,
 'legend.fontsize':7,'legend.frameon':True,'legend.facecolor':'white','legend.edgecolor':'white',
 'legend.framealpha':.95,'axes.titlesize':8.5,'xtick.labelsize':7,'ytick.labelsize':7})
fig,axes=plt.subplots(2,2,figsize=(7.5,6.1))
fig.subplots_adjust(left=.10,right=.98,bottom=.10,top=.92,wspace=.39,hspace=.55)
for letter,ax in zip('ABCD',axes.flat):
 ax.text(-.18,1.13,letter,transform=ax.transAxes,fontweight='bold',fontsize=12)
a=axes[0,0]
source=pd.read_csv(ROOT/'data/miller1996_table1.csv')
row=source[(source.protein=='WT')&(source.activity=='phosphatase')].iloc[0].to_dict()
b=np.geomspace(.01,100,200)*row['Km_IDH_form_uM'];atp=row['Km_ATP_uM']
reference=np.array([law(atp,s,row) for s in b])
random=np.array([stationary(network(atp,s,row,'random',1000))[-1] for s in b])
a.plot(b,reference,color=GRAY,lw=2.5,label='Published rate-law summary')
a.plot(b,reference,'--',color=BLUE,label='Ordered finite mechanism')
a.plot(b,random,':',color=ORANGE,lw=2,label='Random finite mechanism')
a.set(xscale='log',xlabel=r'Phospho-IDH ($\mu$M)',ylabel=r'Initial rate / $V_{\max}$',
      title='Distinct mechanisms, similar initial rates')
a.legend(loc='upper left')
a=axes[0,1];limits=pd.read_csv(ROOT/'results/kinetic_limit.csv')
for protein,color in [('WT',BLUE),('AceK3',ORANGE),('AceK4',GREEN)]:
 for activity,style in [('kinase','-'),('phosphatase','--')]:
  d=limits[(limits.protein==protein)&(limits.activity==activity)]
  a.plot(d.binding_to_catalysis_ratio,d.max_random_relative_error,style,color=color,
         marker='o',ms=3,label=protein if activity=='kinase' else None)
a.axhline(.01,color=GRAY,ls=':',lw=.9)
a.text(.04,.12,'Solid: kinase; dashed: phosphatase',transform=a.transAxes,fontsize=6.8)
a.set(xscale='log',yscale='log',xlabel='Binding dissociation / catalytic rate',
      ylabel='Largest relative rate-law deviation',title='Finite kinetics quantify the approximation')
a.legend(loc='upper right')
a=axes[1,0];design=pd.read_csv(ROOT/'results/binding_design.csv')
for protein,color in [('WT',BLUE),('AceK3',ORANGE),('AceK4',GREEN)]:
 for activity,style in [('kinase','-'),('phosphatase','--')]:
  d=design[(design.protein==protein)&(design.activity==activity)&(design.protein_over_apparent_Kd==1)&(design.ATP_over_Km>0)]
  a.plot(d.ATP_over_Km,d.occupancy_contrast,style,color=color,marker='o',ms=3,
         label=protein if activity=='kinase' else None)
a.set(xscale='log',xlabel=r'ATP / $K_{m,\mathrm{ATP}}$',ylabel='Protein-bound fraction contrast',
      ylim=(0,.55),title='Binding reveals the ambiguity')
a.legend(loc='lower left')
a=axes[1,1]
bounds=pd.DataFrame(json.loads((ROOT/'results/verification.json').read_text())['observation_gate_checks'])
for key,color,label in [('generic_width',ORANGE,'General enclosure'),('conditioned_width',BLUE,'With phosphate-count gate')]:
 a.plot(bounds.AceK_uM,bounds[key],'o-',color=color,ms=3,label=label)
a.annotate('0.0570',(.05,float(bounds[bounds.AceK_uM==.05].generic_width.iloc[0])),
           xytext=(-8,10),textcoords='offset points',fontsize=7,color=ORANGE)
a.annotate('0.0363',(.05,float(bounds[bounds.AceK_uM==.05].conditioned_width.iloc[0])),
           xytext=(4,-13),textcoords='offset points',fontsize=7,color=BLUE)
a.axhline(.05,color=GRAY,ls='--',lw=.9)
a.set(xscale='log',yscale='log',xlabel=r'AceK total ($\mu$M)',ylabel='Class interval width',
      title='An observation tightens the criterion',ylim=(.016,.19))
a.set_xticks([.02,.05,.1],['0.02','0.05','0.10'])
a.xaxis.set_minor_formatter(NullFormatter())
a.legend(loc='upper left')
fig.savefig(OUT/'Fig_AceK_generalization.pdf')
fig.savefig(OUT/'Fig_AceK_generalization.png',dpi=220)
plt.close(fig)
print(OUT)
