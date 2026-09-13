"""Rendering-only style derived from the user's zeolite zstyle.m/zax.m.
Data arrays, limits, normalizations and uncertainty definitions are untouched.
"""
import hashlib,json
from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.text import Text
from matplotlib.colors import LinearSegmentedColormap,to_rgba
from matplotlib.font_manager import fontManager
B=Path(__file__).resolve().parents[1]
INK='#0b0b0b';AX='#52514e';MUTED='#898781'
BLUE='#0690f6';RED='#ff7b02';GREEN='#2bd302';AMBER='#ff1f1f';PURPLE='#ff7b02';PINK='#0e56ff'
seq=LinearSegmentedColormap.from_list('figure1_blue',['#fcfcfb','#d0e8fd','#7abff9',BLUE,'#054778'])
div=LinearSegmentedColormap.from_list('figure1_div',[BLUE,'#fcfcfb',RED])
try:fontManager.addfont('/System/Library/Fonts/Helvetica.ttc')
except Exception:pass
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Helvetica','Arial','DejaVu Sans'],'mathtext.fontset':'dejavusans','pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':False,'axes.edgecolor':AX,'text.color':INK,'axes.labelcolor':INK,'xtick.color':AX,'ytick.color':AX,'axes.linewidth':.9,'xtick.direction':'out','ytick.direction':'out','legend.frameon':False})
colormap={}
for target,old in [(BLUE,['#4d88cf','#427ba8','#477c9b','#3f7193','#4d83bd','#2f5268','#4f82bd']), (RED,['#d87978','#b77a75','#b87973','#d37a73']), (GREEN,['#33a99e','#42978b','#7b9778','#75947a']), (AMBER,['#d69b41','#d28b50','#a59062','#aa9768']), (PURPLE,['#8764a1']), (PINK,['#887eaa']), (MUTED,['#cd cfd2'.replace(' ',''),'#cd cfd2'.replace(' ','')])]:
 for c in old:colormap[c]=target

def recolor(c):
 try:
  r=to_rgba(c);hx=mpl.colors.to_hex(r).lower()
  if hx in colormap:return to_rgba(colormap[hx],r[3])
  return c
 except Exception:return c

def fingerprint(fig):
 h=hashlib.sha256()
 for ax in fig.axes:
  for line in ax.lines:
   for dat in [line.get_xdata(),line.get_ydata()]:
    a=np.asarray(dat);h.update(str(a.shape).encode());h.update(repr(a.tolist()).encode())
  for im in ax.images:h.update(np.asarray(im.get_array()).tobytes())
  for co in ax.collections:
   try:h.update(np.asarray(co.get_offsets()).tobytes())
   except Exception:pass
   if co.get_array() is not None:h.update(np.asarray(co.get_array()).tobytes())
 return h.hexdigest()

def restyle(fig):
 if getattr(fig,'_journal_done',False):return
 fig._journal_done=True;before=fingerprint(fig)
 fig.set_facecolor('white')
 data_axes=[a for a in fig.axes if getattr(a,'_colorbar',None) is None]
 n=len(data_axes)
 if n==4:fig.set_size_inches(10.4,8.8)
 elif n==2:fig.set_size_inches(10.4,4.5)
 elif n==3:fig.set_size_inches(14.8,4.6)
 elif n==1:
  fig.set_size_inches(6.2,5.5) if data_axes[0].name=='3d' else fig.set_size_inches(4.2,3.8)

 for ax in fig.axes:
  ax.set_facecolor('white');ax.grid(False);ax.tick_params(which='both',direction='out',colors=AX,labelsize=10,width=.8)
  for key,sp in ax.spines.items():sp.set_color(AX);sp.set_linewidth(.8)
  if ax.name!='3d':
   ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
  else:
   for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.fill=False;axis.pane.set_edgecolor('#e1e0d9');axis.label.set_fontsize(11)
  ax.xaxis.label.set_fontsize(11);ax.yaxis.label.set_fontsize(11)
  for title in [ax.title,ax._left_title,ax._right_title]:title.set_fontsize(11);title.set_fontweight('normal')
  for line in ax.lines:
   line.set_color(recolor(line.get_color()));line.set_markeredgecolor(recolor(line.get_markeredgecolor()));line.set_markerfacecolor(recolor(line.get_markerfacecolor()))
   if line.get_linewidth()>=1.2:line.set_linewidth(2)
  for co in ax.collections:
   if getattr(ax,'_colorbar',None) is not None:continue
   if co.get_array() is not None:
    name=co.get_cmap().name;co.set_cmap(div if name in ['BrBG','coolwarm','RdBu_r','RdBu'] else seq)
    if co.__class__.__name__ in ['QuadMesh','QuadContourSet','Poly3DCollection']:co.set_rasterized(True)
   else:
    for kind in ['facecolor','edgecolor']:
     try:
      cs=getattr(co,'get_'+kind)();getattr(co,'set_'+kind)([recolor(c) for c in cs])
     except Exception:pass
  for im in ax.images:im.set_cmap(div if im.get_cmap().name in ['BrBG','coolwarm','RdBu','RdBu_r'] else seq)
  for pa in ax.patches:pa.set_facecolor(recolor(pa.get_facecolor()));pa.set_edgecolor(recolor(pa.get_edgecolor()))
  le=ax.get_legend()
  if le:
   le.set_frame_on(False)
   for t in le.get_texts():t.set_fontsize(9)
   for h in le.legend_handles:
    if hasattr(h,'get_color'):h.set_color(recolor(h.get_color()))
    if hasattr(h,'get_facecolor'):h.set_facecolor(recolor(h.get_facecolor()))
    if hasattr(h,'get_edgecolor'):h.set_edgecolor(recolor(h.get_edgecolor()))
 for t in fig.findobj(Text):
  t.set_fontfamily('Helvetica');t.set_color(recolor(t.get_color()))
  if t.get_text() in list('ABCD') and t.get_fontweight() in ['bold',700]:t.set_fontsize(15)
 if fig.get_layout_engine() is None:
  try:fig.set_layout_engine('constrained')
  except Exception:pass
 after=fingerprint(fig);assert before==after,'Rendering modified scientific data'
 fig._journal_hash=before

original=Figure.savefig

def savefig(self,fname,*args,**kwargs):
 restyle(self)
 kwargs['dpi']=600
 if str(fname).endswith('.pdf'):
  log=B/'figure_data_integrity.json';d=json.loads(log.read_text()) if log.exists() else {}
  d[str(fname)]={'before_equals_after':True,'artist_data_sha256':self._journal_hash};log.write_text(json.dumps(d,indent=2))
 return original(self,fname,*args,**kwargs)
Figure.savefig=savefig
