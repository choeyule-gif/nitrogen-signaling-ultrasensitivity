"""Independent algebraic checks of the sufficient bounds and hidden-state example."""
import numpy as np,json
from pathlib import Path
from model import Cascade
rng=np.random.default_rng(20260915)
err=0.;violation=0.
for _ in range(10000):
 w=np.exp(rng.normal(size=4));a=np.exp(rng.normal(size=4));b=np.exp(rng.normal(size=4));h=np.exp(rng.normal(size=4));f=rng.uniform(.001,.999);j=1
 ww=np.r_[w[:j],f*w[j],(1-f)*w[j],w[j+1:]]
 def split(v):return np.r_[v[:j],v[j]/f,0,v[j+1:]]
 err=max(err,abs(ww.sum()/w.sum()-1),abs(ww@split(a)/(w@a)-1),abs(ww@split(b)/(w@b)-1),abs(ww@split(h)/(w@h)-1))
 # Arbitrary positive GS ladder coefficients, not just the simulated J family.
 coeff=rng.normal(0,4,size=13);L=rng.uniform(.001,3);t=rng.normal(0,3)
 def mean(v):
  z=coeff+np.arange(13)*v;z=np.exp(z-z.max());return z@np.arange(13)/(12*z.sum())
 violation=max(violation,abs(mean(t+L)-mean(t))-np.tanh(3*L))
assert err<1e-12 and violation<1e-12
out=Path(__file__).resolve().parents[1]/'results'
(out/'algebra_verification.json').write_text(json.dumps(dict(draws=10000,hidden_state_relative_error=err,arbitrary_ladder_bound_violation=violation),indent=2))
for name,m in [('reference',Cascade()),('productive',Cascade(productive=1))]:
 data=dict(species=m.names,conserved_proteins=['PII','GS','GlnD','GlnE'],conservation_matrix=m.cons.tolist(),stoichiometric_matrix=m.N.tolist(),reactant_indices=[[int(i) for i in r] for r in m.rs],rate_constants=m.ks.tolist(),free_glutamine_uM=m.G,units='micromolar and arbitrary time units; buffered glutamine enters pseudo-first-order rates')
 (out/(name+'_network.json')).write_text(json.dumps(data))
print('Algebra and network export verified',err,violation)
