"""Main Figures 3 and 4. All input values are read from disclosed results."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,LogNorm
from observation_certificate import certificate
from binding_observation import occupancy
R=Path(__file__).resolve().parents[1];F=R/'figures';F.mkdir(exist_ok=True)
# This fallback is only for the author's staging folder; installed code uses sibling layer.
TOPO=R.parent/'closed_cascade/results/productive_topology.csv'
if not TOPO.exists():TOPO=R/'baseline/S1_Data/extras/closed_cascade/results/productive_topology.csv'
BLUE='#0690F6';ORANGE='#FF7B02';GREEN='#2BD302';RED='#FF1F1F';NAVY='#0E56FF';GRAY='#797979'
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'font.size':8.5,'axes.labelsize':8.5,'axes.titlesize':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'xtick.direction':'out','ytick.direction':'out','xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'ps.fonttype':42,'lines.linewidth':1.5,'legend.fontsize':7,'legend.frameon':True,'legend.facecolor':'white','legend.edgecolor':'white','legend.framealpha':.95})
def canvas():
 fig,axs=plt.subplots(2,2,figsize=(7.5,6.1));fig.subplots_adjust(left=.09,right=.91,bottom=.09,top=.93,wspace=.42,hspace=.56)
 for lab,a in zip('ABCD',axs.flat):a.text(-.19,1.14,lab,transform=a.transAxes,fontsize=12,fontweight='bold')
 return fig,axs
def save(fig,n):
 fig.savefig(F/f'Fig{n}.pdf');fig.savefig(F/f'Fig{n}.png',dpi=180);fig.savefig(F/f'Fig{n}.tif',dpi=600,pil_kwargs={'compression':'tiff_lzw'});plt.close(fig)
def readcsv(n):
 return list(csv.DictReader((R/'results'/n).open()))
q=json.loads((R/'results/observation_certificate.json').read_text());cfg=q['parameters'];fig,axs=canvas();a=axs[0,0]
es=np.geomspace(.002,.15,220);res=[certificate(**(cfg|{'Et':e}))for e in es]
a.plot(es,[x['generic_bound']for x in res],color=ORANGE,label='General bound');a.plot(es,[x['bound']for x in res],color=BLUE,label='With PII count');a.plot(es,es/cfg['St'],':',color=GRAY,label=r'$E_T/S_T$');a.axhline(.05,color='black',ls='--',lw=.9,label='0.05 tolerance')
for e in [.02,.1]:
 val=certificate(**(cfg|{'Et':e}))['bound'];a.plot(e,val,'o',color=BLUE,ms=4)
 a.annotate(f'{val:.3f}',(e,val),xytext=(2,7),textcoords='offset points',fontsize=7,color=NAVY)
a.set(xscale='log',yscale='log',xlabel=r'GlnE total $E_T$ ($\mu$M)',ylabel='Total GS difference: upper bound',ylim=(.003,3.),title='A count gate certifies a lower enzyme total');a.legend(loc='upper left',ncol=2,fontsize=6.5,columnspacing=.8,handlelength=1.6)
a=axs[0,1];m=np.linspace(.05,1.05,151);kp=np.geomspace(.08,250,141);z=np.array([[certificate(**(cfg|{'m_upper':mi,'KP':k}))['bound']for mi in m]for k in kp])
cmap=LinearSegmentedColormap.from_list('blue_orange',['#e7f5ff',BLUE,NAVY,ORANGE,RED])
im=a.pcolormesh(m,kp,z,cmap=cmap,norm=LogNorm(.04,.5),shading='auto',rasterized=True);cs=a.contour(m,kp,z,levels=[.05],colors=['black'],linewidths=1);a.clabel(cs,fmt={.05:'0.05'},fontsize=7)
for k,marker,label in [(.18,'o','Main text'),(1,'s','Sensitivity'),(180,'^','SI literal')]:
 a.scatter([.67],[k],marker=marker,s=26,facecolor='white',edgecolor='black',linewidth=.7,zorder=5);a.text(.70,k,label,fontsize=6.8,va='center',bbox={'fc':'white','ec':'none','alpha':.75,'pad':.6})
a.set(yscale='log',xlabel=r'PII count upper limit $M$ (UMP/trimer)',ylabel=r'PII coefficient $K_P$ ($\mu$M)',title='The observation and units matter');cb=fig.colorbar(im,ax=a,fraction=.045,pad=.03);cb.set_ticks([.04,.05,.1,.2,.4]);cb.set_ticklabels(['.04','.05','.10','.20','.40']);cb.ax.tick_params(labelsize=7);cb.ax.minorticks_off()
a=axs[1,0];rows=json.loads((R/'results/certificate_full_reaction_checks.json').read_text())['rows'];x=np.array([r['actual_gap']for r in rows]);y=np.array([r['observation_bound']for r in rows]);g=np.array([r['generic_bound']for r in rows])
a.scatter(x,g,s=12,marker='x',lw=.6,color=ORANGE,alpha=.65,label='General');a.scatter(x,y,s=13,marker='o',facecolors='none',edgecolors=BLUE,lw=.7,label='Count-conditioned');lo=max(1e-16,x[x>0].min()/2);a.set(xscale='log',yscale='log',xlabel='Actual pairwise total GS difference',ylabel='Upper bound',xlim=(lo,x.max()*3),ylim=(.002,1.1),title='Bounds for 120 paired realizations');a.legend(loc='lower right')
a=axs[1,1];t=np.genfromtxt(TOPO,delimiter=',',names=True);eta=np.unique(t['eta'])
for key,col,sty,label in [('fixed_free_GS_gap',NAVY,'-','Matched free pools'),('closed_GS_gap',ORANGE,'--','Conserved totals')]:a.plot(eta,[t[key][t['eta']==e].max()for e in eta],sty,color=col,label=label)
a.set_xscale('symlog',linthresh=.0001);a.set(xlabel=r'Added productive EPU coefficient $\eta$',ylabel='Maximum GS fraction difference',title='Productive EPU changes equivalence');a.set_xticks([0,.001,.01,.1,1],['0',r'$10^{-3}$',r'$10^{-2}$',r'$10^{-1}$','1']);a.legend(loc='upper left');a.text(.05,.54,r'$\alpha_2=2.17$ versus $8.85$',transform=a.transAxes,fontsize=7.5)
save(fig,3)
fig,axs=canvas();a=axs[0,0];ex=json.loads((R/'results/binding_observation_example.json').read_text());P=np.array(ex['P_nominal_uM']);obs=np.array(ex['observations']);h=np.array(ex['H_estimated_uM']);e=ex['E_measured_uM'];v=occupancy(P[None,:],h[:,None],e);gain,bg=np.linalg.lstsq(np.column_stack([np.repeat(v.ravel(),2),np.ones(obs.size)]),obs.ravel(),rcond=None)[0];xx=np.r_[0,np.geomspace(.001,5,180)]
for j,(color,label)in enumerate(zip([BLUE,ORANGE,GREEN],['0 mM Gln','15.6 mM','31.2 mM'])):
 a.scatter(np.repeat(P,2),obs[j].ravel(),s=9,edgecolors=color,facecolors='none',lw=.6);a.plot(xx,bg+gain*occupancy(xx,h[j],e),color=color,label=label)
a.set_xscale('symlog',linthresh=.005);a.set(xlabel=r'Total unmodified PII ($\mu$M)',ylabel='Binding signal (normalized units)',title='Dose points to fitted halfpoints',ylim=(-.08,1.1));a.set_xticks([0,.01,.1,1,5],['0','.01','.1','1','5']);a.legend(loc='lower right')
a=axs[0,1];main=readcsv('binding_observation.csv');extra=readcsv('binding_local_effects.csv')
for kg,col,marker in [(7.8,BLUE,'o'),(11.7,ORANGE,'s'),(15.6,GREEN,'^')]:
 d=[r for r in main if r['scenario']=='dose_bias_corrected' and float(r['KG_mM'])==kg and float(r['alpha'])<=.25]+[r for r in extra if float(r['KG_mM'])==kg];d.sort(key=lambda r:float(r['alpha']));xx=np.array([float(r['alpha'])for r in d]);yy=np.array([float(r['rejection_rate'])for r in d]);lo=np.array([float(r['interval_low'])for r in d]);hi=np.array([float(r['interval_high'])for r in d]);a.errorbar(xx,yy,yerr=[np.maximum(0,yy-lo),np.maximum(0,hi-yy)],color=col,marker=marker,ms=3,capsize=1.5,lw=1,label=fr'$K_G={kg:g}$ mM')
a.axhline(.05,color=GRAY,ls='--',lw=.8);a.set(xlabel=r'True GlnE regulatory coefficient $\alpha_1$',ylabel='Reference rejection rate',title='Small effects depend on nuisance values',ylim=(-.04,1.04));a.set_xticks([.17,.19,.21,.23,.25]);a.legend(loc='lower right')
a=axs[1,0];names=['ideal','calibration_noise','dose_bias','dose_bias_corrected','hidden_conformation'];labels=['Ideal\ncalibration','Stock +\nenzyme','Dose\nbias','Dose\ncorrected','Hidden\nmixture'];cols=[GRAY,BLUE,ORANGE,GREEN,RED]
for j,(name,col) in enumerate(zip(names,cols)):
 r=next(r for r in main if r['scenario']==name and float(r['alpha'])==.17 and float(r['KG_mM'])==15.6);v=float(r['rejection_rate']);l=float(r['interval_low']);u=float(r['interval_high']);a.errorbar(j,v,yerr=[[max(0,v-l)],[max(0,u-v)]],fmt='o',color=col,ms=4,capsize=3,lw=1.2)
a.axhline(.05,color=GRAY,ls='--',lw=.8);a.set_xticks(range(5),labels,fontsize=7);a.set(ylim=(-.04,1.05),xlim=(-.5,4.5),ylabel=r'Rejection at reference $\alpha_1=0.17$',title='Controls change the inference');a.set_yticks([0,.25,.5,.75,1]);
a=axs[1,1];pr=readcsv('glnd_constraint_profile.csv');x=np.array([float(r['K1_mM'])for r in pr]);y=np.array([float(r['RMSE'])for r in pr]);z=np.array([float(r['max_p1_separation'])for r in pr]);a.plot(x,y,'o-',color=BLUE,ms=3,label='Titration RMSE');a.set(xscale='log',xlabel=r'Fixed GlnD binding scale $K_1$ (mM)',ylabel='RMSE (UMP/trimer)',title='Better fits can reduce assay contrast');a.tick_params(axis='y',colors=BLUE);a.yaxis.label.set_color(BLUE);b=a.twinx();b.spines['right'].set_visible(True);b.plot(x,z,'s--',color=ORANGE,ms=3,label='Occupancy contrast');b.set(ylabel='Maximum singly liganded contrast',ylim=(0,.5));b.tick_params(axis='y',colors=ORANGE);b.yaxis.label.set_color(ORANGE);a.legend(loc='upper right',bbox_to_anchor=(1.,.95));b.legend(loc='upper right',bbox_to_anchor=(1.,.81));save(fig,4)
print(F)
