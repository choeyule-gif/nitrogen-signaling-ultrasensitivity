from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Helvetica','Arial','DejaVu Sans'],'font.size':11,'axes.labelsize':11,'axes.titlesize':12,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#52514e','axes.linewidth':.8,'xtick.direction':'out','ytick.direction':'out','pdf.fonttype':42,'ps.fonttype':42,'lines.linewidth':2})
c=['#4a3aa7','#3987e5','#008300','#eda100','#d55181']
r=np.genfromtxt(ROOT/'results/cascade_sweep.csv',delimiter=',',names=True)
t=np.genfromtxt(ROOT/'results/productive_topology.csv',delimiter=',',names=True)
def envelope(d,J,col):
 z=r[(r['direct']==d)&(r['J']==J)];x=np.unique(z['epsilon']);return x,np.array([z[col][z['epsilon']==e].max() for e in x])
fig,axs=plt.subplots(2,2,figsize=(10.4,8.2));fig.subplots_adjust(left=.09,right=.98,bottom=.09,top=.93,wspace=.32,hspace=.5)
a=axs[0,0];x,y=envelope(1,0,'PII_TV');a.loglog(x,y,color=c[0],label='Total PII states');a.loglog(x,x,'--',color='#777777',label=r'Conservation bound $\epsilon_P$');a.set(xlabel=r'Sequestration capacity $\epsilon_P$',ylabel='Maximum total-variation distance',title='PII equivalence after protein conservation');a.legend(frameon=True,facecolor='white',edgecolor='white',framealpha=.95,fontsize=9,loc='upper left');a.text(.97,.05,'Free PII state distributions coincide',transform=a.transAxes,fontsize=9,color=c[0],ha='right')
a=axs[0,1]
for d,col,label in [(1,c[1],'Direct Gln route present'),(0,c[2],'Direct Gln route absent')]:
 x,y=envelope(d,0,'GS_total_gap');a.loglog(x,np.maximum(y,1e-16),linestyle='-' if d else '--',color=col,label=label)
x,y=envelope(1,0,'GS_bound');a.loglog(x,y,'--',color='#777777',label=r'Independent-site upper bound');a.set(xlabel=r'Sequestration capacity $\epsilon_P$',ylabel='Maximum total GS fraction difference',title='Total GS drift remains small in this example');a.legend(frameon=True,facecolor='white',edgecolor='white',framealpha=.95,fontsize=8.5,loc='upper left')
a=axs[1,0]
for J,col,label in [(0,c[1],'Independent GS'),(.4,c[4],'Interacting GS (J = 0.4)')]:
 x,y=envelope(1,J,'GS_free_gap');a.loglog(x,np.maximum(y,1e-16),linestyle='-' if d else '--',color=col,label=label)
x,y=envelope(1,.4,'GS_bound');a.loglog(x,np.tanh(3*(-np.log1p(-x))),'--',color='#777777',label=r'All positive 12-site ladders: bound');a.set(xlabel=r'Sequestration capacity $\epsilon_P$',ylabel='Maximum free GS fraction difference',ylim=(1e-10,2),title='GS interactions amplify free-pool drift');a.legend(frameon=True,facecolor='white',edgecolor='white',framealpha=.95,fontsize=8.5,loc='upper left')
a=axs[1,1];x=np.unique(t['eta'])
for key,col,style,label in [('fixed_free_GS_gap',c[0],'-','Matched free protein pools'),('closed_GS_gap',c[3],'--','Conserved protein totals')]:
 y=np.array([t[key][t['eta']==e].max() for e in x]);a.plot(x,y,style,color=col,label=label)
a.set_xscale('symlog',linthresh=.0001);a.set(xlabel=r'Added productive EPU coefficient $\eta$',ylabel='Maximum GS fraction difference',title='A productive route breaks silent-state equivalence');a.legend(frameon=True,facecolor='white',edgecolor='white',framealpha=.95,fontsize=9,loc='upper left');a.text(.05,.47,r'$\alpha_2=2.17$ versus $8.85$',transform=a.transAxes,fontsize=10);a.set_xticks([0,.001,.01,.1,1],['0',r'$10^{-3}$',r'$10^{-2}$',r'$10^{-1}$','1'])
for lab,a in zip('ABCD',axs.flat):a.text(-.17,1.12,lab,transform=a.transAxes,fontweight='bold',fontsize=15)
fig.savefig(OUT/'Fig11.pdf');fig.savefig(OUT/'Fig11.png',dpi=220);fig.savefig(OUT/'Fig11.tif',dpi=400,pil_kwargs={'compression':'tiff_lzw'})
