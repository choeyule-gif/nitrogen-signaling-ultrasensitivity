"""Author-derived Markov-chain helpers; see audit report for interpretation."""
import numpy as np

def chain(n, A=1.0, Bh=0.4, Gh=0.4, Kh=1.0, l=1.0):
    """Auxiliary chain at the identical-independent steady state.  Order: S_0..S_n, B_*, G_*."""
    m = 3*n+1
    b2 = b3 = g2 = g3 = 1.0
    Bhi = [(n-i)*Bh for i in range(n)]
    Ghi = [(i+1)*Gh for i in range(n)]
    b1 = [Bhi[i]*(b2+b3) for i in range(n)]
    g1 = [Ghi[i]*(g2+g3) for i in range(n)]
    e = 1.0; eps = A*e*l
    S = lambda k: k; B = lambda k: (n+1)+k; G = lambda k: (n+1)+n+k
    L = np.zeros((m, m))                       # row generator: L[from, to]
    for i in range(n):
        L[S(i),   B(i)] += b1[i]*e
        L[B(i),   S(i)] += b2
        L[B(i), S(i+1)] += b3
        L[S(i+1), G(i)] += g1[i]*eps
        L[G(i), S(i+1)] += g2
        L[G(i),   S(i)] += g3
    np.fill_diagonal(L, -L.sum(1))
    # stationary distribution
    w, V = np.linalg.eig(L.T)
    pi = np.real(V[:, np.argmin(np.abs(w))]); pi = np.abs(pi); pi /= pi.sum()
    q = -np.diag(L).copy()
    c = float(pi @ q); L, q = L/c, q/c            # mean exit rate 1
    phi = np.array([1.0]*(n+1) + [-1.0]*(2*n))
    return L, q, pi, phi, m

def ReF(L, q, pi, phi, m, om):
    P = np.diag(1.0/q) @ (L + np.diag(q))
    nu = pi*q; Sm = (np.eye(m)+P)/2
    T = np.diag(om/q); d = float(pi@phi); r = d*phi - np.ones(m)
    R2 = 1.0 - d*d
    if om == 0: return R2
    Pp = np.eye(m) - np.outer(phi, nu*phi)
    Am = 2*Sm + 1j*(Pp @ T @ (np.eye(m) - np.outer(phi, pi*phi)))
    g = Pp @ (T @ r)
    Dn = np.diag(np.sqrt(nu)); Di = np.diag(1.0/np.sqrt(nu))
    pt = Dn@phi; pt = pt/np.linalg.norm(pt)
    Zt = np.linalg.qr(np.column_stack([pt] + list(np.eye(m).T)))[0][:, 1:m]
    h = Di @ (Zt @ np.linalg.solve(Zt.conj().T@(Dn@Am@Di)@Zt, Zt.conj().T@(Dn@g)))
    return R2 + float(np.imag(np.sum(pi*np.conj(r)*h)))

def circ_threshold(L, q, pi, m):
    Jm = np.outer(pi, np.ones(m))*L - (np.outer(pi, np.ones(m))*L).T
    Pi = np.diag(1.0/np.sqrt(pi))
    return 0.5*np.linalg.norm(Pi@Jm@Pi, 2)
