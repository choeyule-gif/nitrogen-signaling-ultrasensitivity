"""Assemble submission files from compiled, reviewed manuscript sources.

No scientific simulations are run. Requires pypdf, Pillow, and pdftoppm.
Build LaTeX PDFs twice before running. Post-publication archive metadata can be
supplied through --zenodo-record without changing the archived manuscript.
"""
from pathlib import Path
import argparse,hashlib,json,re,shutil,subprocess,zipfile
from PIL import Image
from pypdf import PdfReader
R=Path(__file__).resolve().parents[1];M=R/'manuscript'
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=R/'output/release_v1.4.0');p.add_argument('--pdftoppm',default=shutil.which('pdftoppm'));p.add_argument('--zenodo-record',type=Path);p.add_argument('--prior-figures',type=Path,default=R/'output/release_v1.3.4/PLOS_CB_final_submission_v1.3.4');a=p.parse_args()
O=a.output.resolve();O.mkdir(parents=True,exist_ok=True);F=O/'PLOS_CB_submission_v1.4.0';F.mkdir(exist_ok=True)
items={'Manuscript.pdf':M/'manuscript_submission.pdf','Manuscript_with_figures.pdf':M/'manuscript_reader.pdf','Cover_letter.pdf':M/'cover_letter.pdf'}
items.update({f'S{i}_Appendix.pdf':M/f'S{i}_Appendix.pdf' for i in range(1,4)})
items.update({f'S{i}_Fig.pdf':M/f'S{i}_Fig.pdf' for i in range(1,10)})
for name,source in items.items():shutil.copy2(source,F/name)
# Keep the author's original figure and the two unchanged result figures.
# Fallback renders their existing vector PDFs if prior TIFFs are unavailable.
for i in range(1,6):
 prior=a.prior_figures/f'Fig{i}.tif'
 if i<4 and prior.exists():shutil.copy2(prior,F/prior.name)
 else:
  if not a.pdftoppm:raise RuntimeError('pdftoppm is required to render figure TIFFs')
  prefix=O/f'render_Fig{i}'
  subprocess.run([a.pdftoppm,'-r','600','-singlefile','-png',str(M/'main_figures'/f'Fig{i}.pdf'),str(prefix)],check=True)
  with Image.open(prefix.with_suffix('.png')) as im:im.convert('RGB').save(F/f'Fig{i}.tif',compression='tiff_lzw',dpi=(600,600))
  prefix.with_suffix('.png').unlink()
subprocess.run([__import__('sys').executable,str(R/'scripts/package_research_data.py'),'--output',str(F/'S1_Data.zip')],check=True)
commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip()
record=json.loads(a.zenodo_record.read_text()) if a.zenodo_record else None
archive=None
if record:
 assert record['metadata']['version']=='1.4.0'
 archive={'doi':record['doi'],'record_id':record['id'],'title':record['metadata']['title'],'files':record['files']}
files={q.name:{'bytes':q.stat().st_size,'sha256':hashlib.sha256(q.read_bytes()).hexdigest()} for q in sorted(F.iterdir()) if q.suffix in ['.pdf','.tif','.zip']}
manifest={'version':'1.4.0','source_commit':commit,'concept_doi':'10.5281/zenodo.22731768','archive':archive,'files':files}
(F/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
guide='''PLOS Computational Biology / source v1.4.0

No journal submission has been made.
For the journal's initial single-PDF option, use Manuscript_with_figures.pdf.
Alternatively upload Manuscript.pdf with Fig1.tif--Fig5.tif separately.
Upload one manuscript version, not both.
Upload S1--S3 Appendices, S1--S9 Figs, and S1_Data.zip as supporting information.
Upload Cover_letter.pdf in the cover-letter category.
MANIFEST.json and this guide are author records, not journal supplements.

The manuscript is double-spaced, including main captions, with continuous
body line numbering. The user's Fig. 1 is preserved. Main figures 4 and 5
are supplied at 600 dpi; all main figures also have vector/PDF sources.
Both concentration examples require prospective matched observation gates.
The old v1.3.4 archive does not contain this new research extension.

Source: https://github.com/choeyule-gif/nitrogen-signaling-ultrasensitivity/releases/tag/v1.4.0
Concept DOI: https://doi.org/10.5281/zenodo.22731768
'''
if archive:guide+='Verified version DOI: https://doi.org/'+archive['doi']+'\n'
(F/'UPLOAD_GUIDE.txt').write_text(guide)
submission=O/'PLOS_CB_submission_v1.4.0.zip'
with zipfile.ZipFile(submission,'w',zipfile.ZIP_DEFLATED,9) as z:
 for q in sorted(F.iterdir()):z.write(q,F.name+'/'+q.name)
# Portable document source; include every referenced figure PDF and local style.
source_files=list(M.glob('*.tex'))+list(M.glob('*.sty'))+list((M/'main_figures').glob('*.pdf'))
for q in M.glob('*.tex'):
 text='\n'.join(line for line in q.read_text().splitlines() if not line.lstrip().startswith('%'))
 for name in re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}',text):
  f=M/name
  if f not in source_files:source_files.append(f)
assert all(q.exists() for q in source_files)
latex=O/'PLOS_CB_LaTeX_v1.4.0.zip'
with zipfile.ZipFile(latex,'w',zipfile.ZIP_DEFLATED,9) as z:
 for q in sorted(source_files):z.write(q,'LaTeX/'+str(q.relative_to(M)))
 z.writestr('LaTeX/README.txt','Compile each .tex file twice with pdfLaTeX, from this directory. manuscript_reader.tex embeds the main figures using the same manuscript source. No BibTeX step is needed: bibliography entries are embedded.\n')
for name,source in [('PLOS_CB_reader_v1.4.0.pdf',M/'manuscript_reader.pdf'),('Cover_letter_v1.4.0.pdf',M/'cover_letter.pdf')]:shutil.copy2(source,O/name)
for zpath in [submission,latex,F/'S1_Data.zip']:
 with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
with zipfile.ZipFile(F/'S1_Data.zip') as z:
 inner=json.loads(z.read('S1_Data/MANIFEST.json'))
 for name,h in inner.items():assert hashlib.sha256(z.read('S1_Data/'+name)).hexdigest()==h
for name,source in items.items():assert (F/name).read_bytes()==source.read_bytes()
for i in [4,5]:
 with Image.open(F/f'Fig{i}.tif') as im:assert im.width==4500 and im.info['dpi']==(600,600)
print(json.dumps({'submission_zip':str(submission),'latex_zip':str(latex),'reader':str(O/'PLOS_CB_reader_v1.4.0.pdf'),'files':len(files),'manuscript_pages':len(PdfReader(F/'Manuscript.pdf').pages),'reader_pages':len(PdfReader(F/'Manuscript_with_figures.pdf').pages),'archive_doi':archive['doi'] if archive else None},indent=2))
