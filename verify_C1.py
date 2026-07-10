"""Proposition C verification: fresh-ring increment positivity + bond factorization.
Structure: ring i fresh, neighbours marginalized. Per-bond block:
  N+ = 1 + prod_k(1+r_k b_k),  N- = 1 + prod_k(1-r_k b_k),  sigma_bond = ln(N+/N-)
  b_k = sin(4 g_k dt) in (0,1) fresh shared; true branch all +1 -> r_k ~ (1+r_k b_k)/2.
Private records: independent KL terms, positive trivially (paper-1 Th2 structure).
"""
import numpy as np
from itertools import product
import sympy as sp

print("=== LEMMA: bond factorization of the local marginal (exact check vs engine) ===")
try:
    import lattice_p3 as eng          # prototype (session file)
except ImportError:
    try:
        import lattice as eng          # production module (src/lattice/)
    except ImportError:
        raise SystemExit("Neither lattice_p3.py nor lattice.py importable from cwd. "
                         "Copy this script into src/lattice/ (or alongside lattice_p3.py). "
                         "Required API: build_chain, all_branches, beta_matrix, "
                         "sample_records, local_marginal_sigma.")
for fn in ("build_chain","all_branches","beta_matrix","sample_records","local_marginal_sigma"):
    assert hasattr(eng, fn), f"engine module lacks {fn}; adapt the import section (names may differ in production lattice.py)"
# chain L=3, ring i=1 has two bonds (deg=2): engine local sigma vs per-bond product formula
G, Tc, owner = eng.build_chain(L=3, Np=2, Ns=2, seed=7)
t = 0.8; S_all = eng.all_branches(3); beta = eng.beta_matrix(G, Tc, t, S_all)
rng = np.random.default_rng(7)
R = eng.sample_records(beta[0], 200, rng)   # true branch all +1 (index 0)
sig_engine = eng.local_marginal_sigma(R, G, Tc, t, owner, S_all, 0, 1, nbrs=[0,1])
# independent per-bond + private reconstruction
priv = [n for n,o in enumerate(owner) if o==('priv',1)]
sig_manual = np.zeros(R.shape[0])
for n in priv:
    eps = np.sin(2*G[1,n]*t)
    sig_manual += np.log((1+R[:,n]*eps)/(1-R[:,n]*eps))
for bond in (0,1):   # link i in {0,1} touches ring 1
    W = [n for n,o in enumerate(owner) if o==('link',bond)]
    b = np.sin(4*G[1,W]*t)                     # equal couplings both rings
    Xp = np.prod(1+R[:,W]*b[None,:], axis=1); Xm = np.prod(1-R[:,W]*b[None,:], axis=1)
    sig_manual += np.log((1+Xp)/(1+Xm))
print(f"max |engine - per-bond-factorized| = {np.abs(sig_engine-sig_manual).max():.2e}")

print()
print("=== C1 symbolic: per-bond increment K=0->1 and K=1->2 ===")
b1, b2 = sp.symbols('b1 b2', positive=True)
def sigma_bond(rs, bs):
    Xp = sp.prod([1+r*b for r,b in zip(rs,bs)])
    Xm = sp.prod([1-r*b for r,b in zip(rs,bs)])
    return sp.log((1+Xp)/(1+Xm))
def Q(rs, bs):
    return sp.prod([(1+r*b)/2 for r,b in zip(rs,bs)])
# K=0->1
D1 = sum(Q([r],[b1]) * sigma_bond([r],[b1]) for r in (1,-1))
D1 = sp.simplify(D1)
print("Delta(0->1) =", D1)
# K=1->2 increment: E_Q[ sigma(r1,r2) - sigma(r1) ]
D2 = sum(Q([r1,r2],[b1,b2]) * (sigma_bond([r1,r2],[b1,b2]) - sigma_bond([r1],[b1]))
         for r1 in (1,-1) for r2 in (1,-1))
D2 = sp.simplify(D2)
print("Delta(1->2) simplified; scanning positivity on (0,1)^2 ...")
f2 = sp.lambdify((b1,b2), D2, 'numpy')
g = np.linspace(1e-4, 1-1e-4, 400)
B1, B2 = np.meshgrid(g, g)
vals = f2(B1, B2)
print(f"min Delta(1->2) on 400x400 grid = {vals.min():.3e}  (at b1={B1.flat[vals.argmin()]:.3f}, b2={B2.flat[vals.argmin()]:.3f})")

print()
print("=== C1 numeric: per-bond increments K->K+1, K up to 6, random b draws ===")
rng = np.random.default_rng(0)
worst = np.inf; worst_cfg = None
for trial in range(20000):
    K = rng.integers(1, 7)
    bs = rng.uniform(1e-3, 1-1e-3, size=K+1)
    rs = np.array(list(product([1,-1], repeat=K+1)))
    q = np.prod((1+rs*bs[None,:])/2, axis=1)
    Xp = np.prod(1+rs*bs[None,:], axis=1); Xm = np.prod(1-rs*bs[None,:], axis=1)
    s_new = np.log((1+Xp)/(1+Xm))
    XpK = np.prod(1+rs[:,:K]*bs[None,:K], axis=1); XmK = np.prod(1-rs[:,:K]*bs[None,:K], axis=1)
    s_old = np.log((1+XpK)/(1+XmK))
    D = (q*(s_new-s_old)).sum()
    if D < worst: worst, worst_cfg = D, (K, bs.copy())
print(f"min increment over 20000 random (K<=6) configs = {worst:.3e}")
print(f"worst config: K={worst_cfg[0]}")

print()
print("=== Full monotonicity of <sigma_i^loc> (private + 2 bonds, random, sequential) ===")
worst_full = np.inf
for trial in range(2000):
    Ka, Kb, Pv = rng.integers(1,4), rng.integers(1,4), rng.integers(1,4)
    seq = []
    for _ in range(Pv): seq.append(('p', rng.uniform(0.01,0.99)))
    for _ in range(Ka): seq.append(('a', rng.uniform(0.01,0.99)))
    for _ in range(Kb): seq.append(('b', rng.uniform(0.01,0.99)))
    rng.shuffle(seq)
    # accumulate expected sigma after each record; increments must all be >= 0
    # private: adds eps*ln((1+eps)/(1-eps)); bonds: per-bond block increments
    got_a, got_b, cur = [], [], 0.0
    for kind, val in seq:
        if kind == 'p':
            inc = val*np.log((1+val)/(1-val))
        else:
            blk = got_a if kind=='a' else got_b
            K = len(blk)
            bs = np.array(blk+[val])
            rs = np.array(list(product([1,-1], repeat=K+1)))
            q = np.prod((1+rs*bs[None,:])/2, axis=1)
            Xp = np.prod(1+rs*bs[None,:], axis=1); Xm = np.prod(1-rs*bs[None,:], axis=1)
            s_new = np.log((1+Xp)/(1+Xm))
            if K:
                XpK = np.prod(1+rs[:,:K]*bs[None,:K],axis=1); XmK = np.prod(1-rs[:,:K]*bs[None,:K],axis=1)
                s_old = np.log((1+XpK)/(1+XmK))
            else:
                s_old = 0.0
            inc = (q*(s_new-s_old)).sum()
            blk.append(val)
        worst_full = min(worst_full, inc)
print(f"min per-record increment over 2000 random sequences = {worst_full:.3e}")
