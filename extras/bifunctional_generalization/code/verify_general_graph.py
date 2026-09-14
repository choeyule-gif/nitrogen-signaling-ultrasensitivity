"""Independent tree-polynomial derivation and enclosure checks beyond AceK."""
from pathlib import Path
import json
import numpy as np
import sympy as sp
from graph_enclosure import tree_coefficients,bernstein_vertices,observation_interval
from acek_certificate import EDGES,source_rates,coefficients

ROOT=Path(__file__).resolve().parents[1]


def main():
    rng=np.random.default_rng(4161927)
    graph_reports=[];max_error=0.;max_violation=0.;conditioned_checks=0
    for n in [3,4,5,6]:
        edges=[]
        for i in range(n):
            edges.append((i,(i+1)%n,[int(rng.integers(1,6)),int(rng.integers(1,6))]))
            edges.append(((i+1)%n,i,[int(rng.integers(1,6)),int(rng.integers(1,6))]))
        coef,exact,q,e=tree_coefficients(n,edges)
        function=sp.lambdify(e,q,'numpy')
        lo,hi=.01,3.
        vertices=bernstein_vertices(coef,lo,hi)
        for x in np.geomspace(lo,hi,101):
            mat=np.array(function(x),float);mat[-1]=1.
            b=np.zeros(n);b[-1]=1.
            direct=np.linalg.solve(mat,b)
            polynomial=coef@x**np.arange(coef.shape[1]);polynomial/=polynomial.sum()
            max_error=max(max_error,float(max(abs(direct-polynomial))))
            # A genuinely conditioned mixture, with a known compatible point.
            reporter=rng.random(n);measurement=np.arange(n)/(n-1)
            bound=rng.dirichlet(np.ones(n));fraction=rng.uniform(0,.2)
            total=(1-fraction)*direct+fraction*bound
            observed=float(measurement@total)
            interval=observation_interval(vertices,.2,reporter,[(measurement,observed-.01,observed+.01)])
            assert interval is not None
            value=float(reporter@total)
            max_violation=max(max_violation,interval[0]-value,value-interval[1])
            conditioned_checks+=1
        graph_reports.append(dict(states=n,edges=len(edges),polynomial_degree=coef.shape[1]-1,coefficient_rank=exact.rank()))
    assert max_error<1e-10 and max_violation<1e-9
    # Independent exact matrix-tree derivation for the published AceK topology.
    k=source_rates()
    edges=[(i,j,[0,k[idx]] if binding else [k[idx]]) for idx,(i,j,binding) in enumerate(EDGES,1)]
    tree,exact,q,e=tree_coefficients(9,edges,[0,1,2])
    hand=coefficients(k)
    difference=float(np.max(abs(tree/tree[0,0]-hand/hand[0,0])))
    assert tree.shape==(3,3) and exact.rank()==3 and difference<1e-10
    # Rank-one graph: equal enzyme factors multiply both directions.
    cancels,rankone,_,_=tree_coefficients(2,[(0,1,[0,2]),(1,0,[0,3])])
    assert rankone.rank()==1
    assert np.max(np.ptp(bernstein_vertices(cancels,0,4),axis=1))<1e-12
    # A contradictory observation must return empty, never a zero-width success.
    assert observation_interval(np.eye(3),.1,[1,0,0],[([1,1,1],1.1,1.2)]) is None
    # A recurrent free subset must not hide a transient excluded state.
    try:
        tree_coefficients(3,[(0,1,[1]),(1,0,[1]),(2,1,[1])],[0,1])
    except ValueError:
        pass
    else:
        raise AssertionError('Reducible full graph accepted')
    try:
        observation_interval(np.array([[.2],[.3]]),.1,[1,0])
    except ValueError:
        pass
    else:
        raise AssertionError('Unnormalized probability vertex accepted')
    result=dict(passed=True,generic_graphs=graph_reports,conditioned_mixture_checks=conditioned_checks,
                max_stationary_polynomial_error=max_error,max_enclosure_violation=max_violation,
                independent_AceK_tree_max_difference=difference,AceK_exact_coefficient_rank=3,
                rank_one_control_passed=True,infeasible_observation_control_passed=True,
                excluded_transient_state_control_passed=True,invalid_vertex_control_passed=True,
                scope='Algorithmic verification for stated positive-polynomial irreducible graphs; no claim that all biochemical networks satisfy these assumptions.')
    (ROOT/'results/general_graph_verification.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
