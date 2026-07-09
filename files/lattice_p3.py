"""
PAPER 3 CORE ENGINE -- lattice of chiral rings with shared recording baths.
Validation-first: gates G1-G4 must pass before any physics run.
Branch structure: joint branches S in {-1,+1}^L, record bias
    beta_n(S) = sin(2 * a_n(S) * (t - Tc_n)),  a_n(S) = sum_i S_i G[i,n]
Fresh witnesses: Tc=0. Charged(mirror) witnesses: Tc=T_charge.
Contrarian ring = its PRIVATE bath charged; shared link baths fresh by default.
"""
import numpy as np
from itertools import product

# ---------------------------------------------------------------- topology
def build_chain(L, Np, Ns, g_scale=0.25, seed=0, contrarians=(), T_charge=1.0,
                shared_strength=1.0):
    """Returns G (L x N couplings), Tc (N,), witness ownership lists."""
    rng = np.random.default_rng(seed)
    cols, Tc, owner = [], [], []
    for i in range(L):                               # private baths
        for _ in range(Np):
            c = np.zeros(L); c[i] = rng.uniform(0.5, 1.5) * g_scale
            cols.append(c); Tc.append(T_charge if i in contrarians else 0.0)
            owner.append(('priv', i))
    for i in range(L - 1):                           # shared link baths
        for _ in range(Ns):
            g = rng.uniform(0.5, 1.5) * g_scale * shared_strength
            c = np.zeros(L); c[i] = g; c[i+1] = g
            cols.append(c); Tc.append(0.0); owner.append(('link', i))
    return np.array(cols).T, np.array(Tc), owner     # G: (L, N)

# ------------------------------------------------------------- record laws
def all_branches(L):
    return np.array(list(product([1,-1], repeat=L)))          # (2^L, L)

def beta_matrix(G, Tc, t, S_all):
    """(n_branches, N): bias of each witness in each joint branch."""
    a = S_all @ G                                             # (B, N)
    return np.sin(2.0 * a * (t - Tc)[None, :])

def sample_records(beta_true, n_traj, rng):
    """records r in {-1,+1}: P(r=+1) = (1+beta)/2 per witness."""
    return np.where(rng.random((n_traj, beta_true.size)) < (1+beta_true)/2, 1, -1)

def loglik(R, beta):
    """(n_traj, B): sum_n log[(1 + r_n * beta_{S,n})/2]."""
    return np.log((1 + R[:, None, :] * beta[None, :, :]) / 2).sum(axis=2)

# --------------------------------------------------------- entropy outputs
def productions(R, beta, S_all, true_idx):
    """Per-trajectory sigma_S and Option-A global marginals sigma_i."""
    LL = loglik(R, beta)                                      # (T, B)
    B, L = S_all.shape[0], S_all.shape[1]
    flip = np.array([np.where((S_all == -S_all[b]).all(axis=1))[0][0]
                     for b in range(B)])
    sig_S = LL[:, true_idx] - LL[:, flip[true_idx]]
    # marginal likelihoods: logsumexp over branches with s_i fixed
    sig_i = np.empty((R.shape[0], L))
    s_true = S_all[true_idx]
    for i in range(L):
        for sgn in (+1, -1):
            mask = S_all[:, i] == sgn
            m = LL[:, mask]; mx = m.max(axis=1, keepdims=True)
            lp = mx.squeeze(1) + np.log(np.exp(m - mx).mean(axis=1))
            if sgn == +1: lp_p = lp
            else: lp_m = lp
        sig_i[:, i] = np.where(s_true[i] == 1, lp_p - lp_m, lp_m - lp_p)
    return sig_S, sig_i

def local_marginal_sigma(R, G, Tc, t, owner, S_all, true_idx, i, nbrs):
    """sigma_i using ONLY ring i's accessible witnesses (private + its links),
    marginalizing uniformly over the neighbour chiralities."""
    W = [n for n, o in enumerate(owner)
         if (o == ('priv', i)) or (o[0] == 'link' and o[1] in nbrs)]
    vars_ = [i] + sorted({j for n in W for j in np.nonzero(G[:, n])[0] if j != i})
    combos = np.array(list(product([1,-1], repeat=len(vars_))))
    a = np.zeros((len(combos), len(W)))
    for k, v in enumerate(vars_):
        a += combos[:, k:k+1] * G[v, W][None, :]
    beta = np.sin(2 * a * (t - Tc[W])[None, :])
    LL = np.log((1 + R[:, None, :][:, :, W] * beta[None, :, :]) / 2).sum(axis=2)
    s_true = S_all[true_idx, i]
    out = np.empty(R.shape[0])
    for sgn_case, store in ((+1, 'p'), (-1, 'm')):
        mask = combos[:, 0] == sgn_case
        m = LL[:, mask]; mx = m.max(axis=1, keepdims=True)
        lp = mx.squeeze(1) + np.log(np.exp(m - mx).mean(axis=1))
        if store == 'p': lp_p = lp
        else: lp_m = lp
    return np.where(s_true == 1, lp_p - lp_m, lp_m - lp_p)

# ======================================================================
# GATES
# ======================================================================
def gate_G2_paper1(tol=1e-12):
    """f=0 (no sharing): <sigma_i> must equal D_i = sum eps*ln((1+e)/(1-e)) EXACTLY
    (exhaustive over ring-i private records)."""
    G, Tc, owner = build_chain(L=2, Np=4, Ns=0, seed=1)
    t = 0.8; S_all = all_branches(2); beta = beta_matrix(G, Tc, t, S_all)
    worst = 0.0
    for i in (0, 1):
        Wi = [n for n, o in enumerate(owner) if o == ('priv', i)]
        eps = np.abs(np.sin(2 * G[i, Wi] * t))
        D = np.sum(eps * np.log((1+eps)/(1-eps)))
        # exhaustive expectation of sigma_i in branch s_i=+1
        E = 0.0
        for r in product([1,-1], repeat=len(Wi)):
            r = np.array(r)
            p = np.prod((1 + r*eps)/2)
            E += p * np.sum(np.log((1+r*eps)/(1-r*eps)))
        worst = max(worst, abs(E - D))
    return worst < tol, worst

def gate_G1_decomposition(tol=1e-12):
    """L=2 shared: Delta = I + C_bwd identity on the LATTICE-generated law,
    exhaustively at small N (paper-2 Theorem, hypotheses: uniform + antipodal)."""
    G, Tc, owner = build_chain(L=2, Np=1, Ns=2, seed=2)
    t = 0.7; S_all = all_branches(2); beta = beta_matrix(G, Tc, t, S_all)
    N = beta.shape[1]
    xs = np.array(list(product([1,-1], repeat=N)))
    P = np.prod((1 + xs[:, None, :] * beta[None, :, :]) / 2, axis=2)   # (X, B)
    Pxi = P.mean(axis=1)
    P1 = {s: P[:, S_all[:,0]==s].mean(axis=1) for s in (1,-1)}
    P2 = {s: P[:, S_all[:,1]==s].mean(axis=1) for s in (1,-1)}
    flip = [np.where((S_all == -S_all[b]).all(axis=1))[0][0] for b in range(4)]
    Delta = I = Cb = 0.0
    for b in range(4):
        s1, s2 = S_all[b]
        w = P[:, b] / 4
        sS = np.log(P[:, b] / P[:, flip[b]])
        si = np.log(P1[s1]/P1[-s1]) + np.log(P2[s2]/P2[-s2])
        J = np.log(P[:, b] * Pxi / (P1[s1] * P2[s2]))
        Delta += (w * (sS - si)).sum()
        I     += (w * J).sum()
        Cb    += -(P[:, flip[b]]/4 * J).sum()
    err = abs(Delta - (I + Cb))
    return err < tol and Delta >= -tol, (err, Delta, I, Cb)

def gate_G3_bruteforce(tol=1e-10):
    """L=3 chain, chirality-qubit representation vs full state-vector expm."""
    from scipy.linalg import expm
    G, Tc, owner = build_chain(L=3, Np=1, Ns=1, seed=3)
    L, N = G.shape; t = 0.6
    sz = np.diag([1.,-1.]).astype(complex); I2 = np.eye(2, dtype=complex)
    def kron(ops):
        out = np.array([[1.]], dtype=complex)
        for o in ops: out = np.kron(out, o)
        return out
    Zi = [kron([sz if k==i else I2 for k in range(L)]) for i in range(L)]
    H = np.zeros((2**L * 2**N,)*2, dtype=complex)
    for n in range(N):
        zn = kron([sz if k==n else I2 for k in range(N)])
        A = sum(G[i, n] * Zi[i] for i in range(L))
        H += np.kron(A, zn)
    plus = np.ones(2)/np.sqrt(2)
    psi0 = kron([plus.reshape(2,1)]*L).ravel()
    psi0 = np.kron(psi0, kron([plus.reshape(2,1)]*N).ravel())
    psi_t = expm(-1j*H*t) @ psi0
    S_all = all_branches(L)
    psi_ana = np.zeros_like(psi_t)
    for b, S in enumerate(S_all):
        a = S @ G
        bath = kron([ (expm(-1j*a[n]*sz*t) @ plus).reshape(2,1) for n in range(N)]).ravel()
        ring = np.zeros(2**L, dtype=complex)
        idx = int(''.join('0' if s==1 else '1' for s in S), 2)
        ring[idx] = 2**(-L/2)
        psi_ana += np.kron(ring, bath)
    err = np.abs(psi_t - psi_ana).max()
    return err < tol, err

def gate_G4_local_vs_global(L=4, n_traj=4000, seed=5):
    """Chain L=4: sign agreement and correlation of local vs global marginals."""
    G, Tc, owner = build_chain(L=L, Np=3, Ns=2, seed=seed)
    t = 0.8; S_all = all_branches(L); beta = beta_matrix(G, Tc, t, S_all)
    rng = np.random.default_rng(seed)
    true_idx = 0                                  # all +1 branch
    R = sample_records(beta[true_idx], n_traj, rng)
    _, sig_g = productions(R, beta, S_all, true_idx)
    agree, corr = [], []
    for i in range(L):
        nbrs = [j for j in (i-1, i) if 0 <= j < L-1]     # links adjacent to i
        sl = local_marginal_sigma(R, G, Tc, t, owner, S_all, true_idx, i, nbrs)
        m_g, m_l = sig_g[:, i].mean(), sl.mean()
        agree.append(np.sign(m_g) == np.sign(m_l))
        corr.append(np.corrcoef(sig_g[:, i], sl)[0,1])
    return all(agree), (np.mean(corr), [f"{c:.3f}" for c in corr])

# ======================================================================
if __name__ == "__main__":
    print("=== PAPER 3 GATES ===")
    ok, v = gate_G2_paper1();        print(f"[{'PASS' if ok else 'FAIL'}] G2 f=0 -> paper-1 Th2 exact: {v:.2e}")
    ok, v = gate_G1_decomposition(); print(f"[{'PASS' if ok else 'FAIL'}] G1 L=2 Delta=I+Cbwd on lattice law: err={v[0]:.2e}, Delta={v[1]:.4f}, I={v[2]:.4f}, Cbwd={v[3]:.4f}")
    ok, v = gate_G3_bruteforce();    print(f"[{'PASS' if ok else 'FAIL'}] G3 L=3 brute force vs branches: {v:.2e}")
    ok, v = gate_G4_local_vs_global(); print(f"[{'PASS' if ok else 'FAIL'}] G4 local vs global marginal: sign-agree={ok}, mean corr={v[0]:.3f} {v[1]}")
