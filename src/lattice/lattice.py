"""Paper 3 — lattice of chiral rings with shared recording baths.

Production engine, ported from the validated prototype (files/lattice_p3.py,
4 gates passed) and extended for Paper 3:
  - build_square: 2D square lattice, 4-neighbour shared-link baths;
  - collision protocol (fixed window dt, permanent monotone records, as in
    Paper-1 src/): the production protocol, kills the recurrences seen in the
    instantaneous preview;
  - local_marginal_sigma made topology-agnostic (accessible witnesses read
    straight from the coupling matrix), so it works on chain and square alike.

Branch structure: joint branches S in {-1,+1}^L, record bias
    beta_n(S) = sin(2 * a_n(S) * (t - Tc_n)),   a_n(S) = sum_i S_i G[i,n].
Fresh witnesses Tc=0; charged (mirror) witnesses Tc=T_charge. A contrarian ring
has its PRIVATE bath charged; shared link baths are fresh by default.

Arrow of a ring (Paper 3 definition): tau_i = sign(d<sigma_i^local>/dt), the
direction of ACCUMULATION of ring i's accessible records (see
notes/paper3/arrow_definition.md). Not the inference-ratio sign of Paper-2 C5.
"""
import numpy as np
from itertools import product


# ---------------------------------------------------------------- topology
def build_chain(L, Np, Ns, g_scale=0.25, seed=0, contrarians=(), T_charge=1.0,
                shared_strength=1.0):
    """1D chain of L rings. Np private witnesses/ring, Ns shared per link.
    Returns G (L x N couplings), Tc (N,), owner list. (Unchanged prototype;
    link owner is ('link', i) meaning the bond (i, i+1).)"""
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
            c = np.zeros(L); c[i] = g; c[i + 1] = g
            cols.append(c); Tc.append(0.0); owner.append(('link', i))
    return np.array(cols).T, np.array(Tc), owner     # G: (L, N)


def build_square(Lx, Ly, Np, Ns, g_scale=0.25, seed=0, contrarians=(),
                 T_charge=1.0, shared_strength=1.0):
    """2D square lattice (Lx*Ly rings, row-major index i = x + Lx*y).
    4-neighbour bonds (right and up); Ns shared witnesses per bond. Link owner
    is ('link', (i, j)) with i < j the coupled ring pair."""
    rng = np.random.default_rng(seed)
    L = Lx * Ly
    cols, Tc, owner = [], [], []
    for i in range(L):                               # private baths
        for _ in range(Np):
            c = np.zeros(L); c[i] = rng.uniform(0.5, 1.5) * g_scale
            cols.append(c); Tc.append(T_charge if i in contrarians else 0.0)
            owner.append(('priv', i))
    bonds = []
    for y in range(Ly):
        for x in range(Lx):
            i = x + Lx * y
            if x + 1 < Lx:
                bonds.append((i, i + 1))             # horizontal
            if y + 1 < Ly:
                bonds.append((i, i + Lx))            # vertical
    for (i, j) in bonds:                             # shared link baths
        for _ in range(Ns):
            g = rng.uniform(0.5, 1.5) * g_scale * shared_strength
            c = np.zeros(L); c[i] = g; c[j] = g
            cols.append(c); Tc.append(0.0); owner.append(('link', (i, j)))
    return np.array(cols).T, np.array(Tc), owner


# ------------------------------------------------------------- record laws
def all_branches(L):
    return np.array(list(product([1, -1], repeat=L)))            # (2^L, L)


def beta_matrix(G, Tc, t, S_all):
    """Instantaneous protocol: (n_branches, N) bias at snapshot time t."""
    a = S_all @ G                                               # (B, N)
    return np.sin(2.0 * a * (t - Tc)[None, :])


def collision_beta(G, Tc, dt, S_all):
    """Collision protocol: FIXED window dt (permanent, monotone records; Paper-1
    src/ convention). Fresh witness bias sin(2 a dt); charged witness carries the
    same magnitude but its records are consumed (drawn backward) -- the draw sign
    lives in draw_signs(), the scoring bias is returned here."""
    a = S_all @ G                                               # (B, N)
    # charged witnesses: residual phase (T_charge - dt) i.e. sin(2 a (dt - Tc))
    return np.sin(2.0 * a * (dt - Tc)[None, :])


def draw_signs(Tc):
    """+1 for fresh witnesses (forward draw), -1 for charged (backward draw)."""
    return np.where(np.asarray(Tc) > 0, -1.0, 1.0)


def sample_records(beta_true, n_traj, rng, draw_sign=None):
    """records r in {-1,+1}: fresh P(r=+1)=(1+beta)/2; charged (backward draw)
    P(r=+1)=(1-beta)/2. beta_true is (N,)."""
    if draw_sign is None:
        draw_sign = np.ones_like(beta_true)
    p_plus = (1.0 + draw_sign * beta_true) / 2.0
    return np.where(rng.random((n_traj, beta_true.size)) < p_plus, 1, -1)


def loglik(R, beta):
    """(n_traj, B): sum_n log[(1 + r_n * beta_{S,n})/2] (forward scoring)."""
    return np.log((1 + R[:, None, :] * beta[None, :, :]) / 2).sum(axis=2)


# --------------------------------------------------------- entropy outputs
def _flip_index(S_all):
    B = S_all.shape[0]
    return np.array([np.where((S_all == -S_all[b]).all(axis=1))[0][0]
                     for b in range(B)])


def productions(R, beta, S_all, true_idx):
    """Per-trajectory sigma_S and Option-A GLOBAL marginals sigma_i."""
    LL = loglik(R, beta)                                        # (T, B)
    L = S_all.shape[1]
    flip = _flip_index(S_all)
    sig_S = LL[:, true_idx] - LL[:, flip[true_idx]]
    sig_i = np.empty((R.shape[0], L))
    s_true = S_all[true_idx]
    for i in range(L):
        for sgn in (+1, -1):
            m = LL[:, S_all[:, i] == sgn]
            mx = m.max(axis=1, keepdims=True)
            lp = mx.squeeze(1) + np.log(np.exp(m - mx).mean(axis=1))
            if sgn == +1:
                lp_p = lp
            else:
                lp_m = lp
        sig_i[:, i] = np.where(s_true[i] == 1, lp_p - lp_m, lp_m - lp_p)
    return sig_S, sig_i


def local_marginal_sigma(R, G, Tc, param, S_true, i):
    """sigma_i using ONLY ring i's accessible witnesses (its private bath + every
    shared link touching it), marginalizing uniformly over the neighbour
    chiralities. SCALABLE: cost is O(2^(deg+1)) per ring, independent of L (it
    never enumerates the 2^L global branches). `S_true` is the true chirality
    vector (only S_true[i] is used); `param` is the snapshot t or collision dt.

    Topology-agnostic: the accessible set is read from the coupling matrix, so
    chain and square use the same code."""
    W = np.nonzero(G[i])[0]                                     # witnesses on i
    vars_ = [i] + sorted({j for n in W for j in np.nonzero(G[:, n])[0] if j != i})
    combos = np.array(list(product([1, -1], repeat=len(vars_))))  # local branches
    a = combos @ G[np.ix_(vars_, W)]                           # (n_combos, |W|)
    beta = np.sin(2 * a * (param - Tc[W])[None, :])
    RW = R[:, W]                                                # (n_traj, |W|)
    LL = np.log((1 + RW[:, None, :] * beta[None, :, :]) / 2).sum(axis=2)
    lp = {}
    for sgn in (+1, -1):                                        # marginalize vars>0
        m = LL[:, combos[:, 0] == sgn]
        mx = m.max(axis=1, keepdims=True)
        lp[sgn] = mx.squeeze(1) + np.log(np.exp(m - mx).mean(axis=1))
    return lp[+1] - lp[-1] if S_true[i] == 1 else lp[-1] - lp[+1]


# --------------------------------------------------- closed-form distinguishab.
def local_D(G, Tc, param, i, protocol="collision"):
    """Closed-form D_i over ring i's PRIVATE witnesses (single-ring, disjoint
    part): sum eps*ln((1+eps)/(1-eps)), eps = |sin(2 G[i,n] (param - Tc_n))|.
    For collision, `param` is the window dt; instantaneous, the snapshot t."""
    Wi = [n for n in np.nonzero(G[i])[0]
          if np.count_nonzero(G[:, n]) == 1]                    # private only
    eps = np.abs(np.sin(2 * G[i, Wi] * (param - Tc[Wi])))
    eps = np.clip(eps, 0, 1 - 1e-15)
    return float(np.sum(eps * np.log((1 + eps) / (1 - eps))))
