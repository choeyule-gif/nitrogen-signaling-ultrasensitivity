"""Reusable positive-polynomial stationary enclosure for finite state graphs.

Edges use (source, destination, coefficient_sequence), with rates expressed
in increasing powers of the free-enzyme concentration. Each edge coefficient
must be nonnegative. Irreducibility for positive input is required.
"""
from math import comb
import numpy as np
import sympy as sp
from scipy.optimize import linprog


def tree_coefficients(states,edges,free_states=None):
    if not isinstance(states,int) or states<1:raise ValueError('Positive integer state count required')
    edges=list(edges)
    forward=[set() for _ in range(states)];reverse=[set() for _ in range(states)]
    e=sp.Symbol('e',positive=True)
    q=sp.zeros(states)
    for source,target,coefficients in edges:
        if not (0<=source<states and 0<=target<states):raise ValueError('State index outside graph')
        cc=[sp.Rational(str(c)) for c in coefficients]
        if any(c<0 for c in cc):raise ValueError('Negative edge coefficient')
        if any(c>0 for c in cc):
            forward[source].add(target);reverse[target].add(source)
        rate=sum(c*e**i for i,c in enumerate(cc))
        q[target,source]+=rate;q[source,source]-=rate
    # Positive edge polynomials have the same support for every e>0.
    # Check all states, including transient states excluded from free_states.
    for adjacency in [forward,reverse]:
        seen={0};pending=[0]
        while pending:
            for target in adjacency[pending.pop()]-seen:
                seen.add(target);pending.append(target)
        if len(seen)!=states:raise ValueError('Graph must be strongly connected for positive input')
    chosen=list(range(states)) if free_states is None else list(free_states)
    if not chosen or len(set(chosen))!=len(chosen) or any(i<0 or i>=states for i in chosen):
        raise ValueError('Invalid free-state subset')
    polys=[sp.Poly((-q).minor_submatrix(i,i).det(method='domain-ge'),e) for i in chosen]
    if any(p.is_zero for p in polys):raise ValueError('Graph is not irreducible on the chosen domain')
    common=polys[0]
    for p in polys[1:]:common=sp.gcd(common,p)
    reduced=[sp.exquo(p,common) for p in polys]
    # Raw tree weights always have nonnegative coefficients. Use a reduced
    # representation only when it preserves that property.
    if all(all(c>=0 for c in p.all_coeffs()) for p in reduced):polys=reduced
    degree=max(p.degree() for p in polys)
    exact=sp.Matrix([[p.nth(j) for j in range(degree+1)] for p in polys])
    if any(c<0 for c in exact):raise ArithmeticError('Tree polynomial positivity violated')
    # A common scalar normalization improves floating-point conditioning.
    exact=exact/max(exact)
    return np.array(exact,dtype=float),exact,q,e


def bernstein_vertices(coef,lower,upper):
    if not (0<=lower<=upper and np.all(coef>=0)):raise ValueError('Invalid positive polynomial interval')
    n=coef.shape[1]-1
    shifted=np.zeros_like(coef,dtype=float)
    for j in range(n+1):
        for l in range(j,n+1):
            shifted[:,j]+=coef[:,l]*comb(l,j)*lower**(l-j)*(upper-lower)**j
    bb=np.zeros_like(shifted)
    for k in range(n+1):
        for j in range(k+1):bb[:,k]+=shifted[:,j]*comb(k,j)/comb(n,j)
    good=bb.sum(axis=0)>0
    return bb[:,good]/bb[:,good].sum(axis=0)


def observation_interval(vertices,epsilon,reporter,observations=()):
    if not (0<=epsilon<=1):raise ValueError('Bound fraction must lie in [0,1]')
    vertices=np.asarray(vertices,dtype=float)
    if vertices.ndim!=2 or not all(vertices.shape) or not np.all(np.isfinite(vertices)) or np.any(vertices<0) or not np.allclose(vertices.sum(axis=0),1,rtol=0,atol=1e-12):
        raise ValueError('Vertices must be nonnegative probability columns')
    n,m=vertices.shape;mapping=np.c_[vertices,np.eye(n)]
    c=np.asarray(reporter)@mapping
    au=[np.r_[np.zeros(m),np.ones(n)]];bu=[epsilon]
    for weights,low,high in observations:
        row=np.asarray(weights)@mapping
        au.extend([row,-row]);bu.extend([high,-low])
    kw=dict(A_ub=au,b_ub=bu,A_eq=np.ones((1,n+m)),b_eq=[1.],bounds=(0,None),method='highs')
    fits=[linprog(s*c,**kw) for s in [1.,-1.]]
    if any(f.status==2 for f in fits):return None
    if not all(f.success for f in fits):raise RuntimeError('Observation LP did not solve')
    return [float(fits[0].fun),float(-fits[1].fun)]
