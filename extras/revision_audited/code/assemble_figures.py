"""Assemble supplied MATLAB panels and promoted checked calculations.
Run from any directory: python code/assemble_figures.py
Dependencies: pypdf, reportlab, pdfplumber.
"""
from pathlib import Path
import io
from pypdf import PdfReader, PdfWriter, Transformation
from reportlab.pdfgen.canvas import Canvas
import pdfplumber
R=Path(__file__).resolve().parents[1]
M=R/'rendered_matlab'
names=['Fig1_architecture','Fig2_sensitivity','Fig3_titration','Fig4_cascade','Fig5_experiments']
P=R/'figures/promoted'
panels={
 'Fig1_architecture':[M/'Fig1A.pdf',M/'Fig1B.pdf',P/'Fig1C_total_heatmap.pdf',P/'Fig1D_total_slices.pdf'],
 'Fig2_sensitivity':[M/'Fig2A.pdf',M/'Fig2B.pdf',M/'Fig2C.pdf',P/'Fig2D_ladder_constraints.pdf'],
 'Fig3_titration':[M/f'Fig3{x}.pdf' for x in 'ABCD'],
 'Fig4_cascade':[M/'Fig4A.pdf',M/'Fig4B.pdf',M/'Fig4C.pdf',P/'Fig4D_direct_route.pdf'],
 'Fig5_experiments':[M/f'Fig5{x}.pdf' for x in 'ABCD'],
}
for name in names:
 panel_paths=panels[name];count=len(panel_paths)
 cols=2; rows=2; pw,ph=320,290; W,H=cols*pw,rows*ph
 writer=PdfWriter(); page=writer.add_blank_page(W,H)
 for i in range(count):
  panel=PdfReader(panel_paths[i]).pages[0]
  box=panel.cropbox; w,h=float(box.width),float(box.height)
  scale=min((pw-12)/w,(ph-24)/h); x=(i%2)*pw+6; y=H-(i//2+1)*ph+6
  panel.add_transformation(Transformation().translate(-float(box.left),-float(box.bottom)))
  page.merge_transformed_page(panel,Transformation().scale(scale).translate(x,y))
 b=io.BytesIO(); c=Canvas(b,pagesize=(W,H)); c.setFont('Helvetica-Bold',15)
 for i in range(count):c.drawString((i%2)*pw+7,H-(i//2)*ph-14,chr(65+i))
 c.save();b.seek(0);page.merge_page(PdfReader(b).pages[0]);writer.write(R/'figures'/f'{name}.pdf')
 with pdfplumber.open(R/'figures'/f'{name}.pdf') as pdf:pdf.pages[0].to_image(resolution=130).save(R/'figures'/f'{name}.png')
