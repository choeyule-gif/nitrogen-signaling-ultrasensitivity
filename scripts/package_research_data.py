"""Package the current research code/data with synchronized documentation.
Run after numerical regeneration. Unlike the legacy baseline S1_Data folder,
this ZIP includes every manuscript research extension.
"""
from pathlib import Path
import argparse,zipfile,json,hashlib
R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=R/'S1_Data.zip');a=p.parse_args()
roots=['esbm','scripts','data','results','figdata','extras','tests','validation','matlab']
files=[]
for base in roots:
 for f in (R/base).rglob('*'):
  if not f.is_file():continue
  if any(x in f.parts for x in ['__pycache__','.pytest_cache','rendered_matlab','panels']):continue
  if f.suffix in ['.pyc','.log','.png','.pdf','.tif','.buildlog']:continue
  if f.name in ['gosztolai2017.xml']:continue # read-only literature copy, not an analysis input
  files.append(f)
for name in ['README.md','REPRODUCIBILITY.md','LICENSE','CITATION.cff','pyproject.toml','requirements.txt','run_all.py','Makefile']:
 if (R/name).exists():files.append(R/name)
manifest={str(f.relative_to(R)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(files)}
a.output.parent.mkdir(parents=True,exist_ok=True)
with zipfile.ZipFile(a.output,'w',zipfile.ZIP_DEFLATED,9) as z:
 for f in sorted(files):z.write(f,Path('S1_Data')/f.relative_to(R))
 z.writestr('S1_Data/MANIFEST.json',json.dumps(manifest,indent=2))
with zipfile.ZipFile(a.output) as z:
 assert z.testzip() is None
 assert z.read('S1_Data/README.md')==(R/'README.md').read_bytes()
print(a.output,len(files),'files',a.output.stat().st_size,'bytes')
