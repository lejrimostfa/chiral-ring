"""Closed-form baseline quantities (pure dephasing, spec §2).

Global state: (1/sqrt2)[|+> prod_n e^{-i g_n t sz}|+x>  +  |-> prod_n e^{+i g_n t sz}|+x>]
All bipartite reduced states are rank <= 2; entropies are binary entropies of
(1 +- |product of cosines|)/2.
"""
import numpy as np

def H2(x):
    """Binary entropy (bits) of eigenvalue pair (x, 1-x)."""
    x = np.clip(x, 1e-15, 1 - 1e-15)
    return -(x * np.log2(x) + (1 - x) * np.log2(1 - x))

def ent_rank2(c):
    """von Neumann entropy (bits) of a rank-2 state with eigenvalues (1+-|c|)/2."""
    return H2((1 + np.abs(c)) / 2)

def mutual_information(g, t, frag_mask):
    """I(S:F) for fragment given by boolean mask over the N bath qubits."""
    cos_all = np.cos(2 * g * t)
    cF = np.prod(cos_all[frag_mask])
    cFbar = np.prod(cos_all[~frag_mask])
    S_S = ent_rank2(cF * cFbar)
    S_F = ent_rank2(cF)
    S_SF = ent_rank2(cFbar)
    return S_S + S_F - S_SF

def partial_info_plot(g, t, n_samples=200, rng=None):
    """Average I(S:F) over random fragments for each size m = 0..N."""
    rng = rng or np.random.default_rng(0)
    N = len(g)
    out = np.zeros(N + 1)
    for m in range(1, N + 1):
        vals = []
        for _ in range(n_samples):
            idx = rng.choice(N, size=m, replace=False)
            mask = np.zeros(N, bool); mask[idx] = True
            vals.append(mutual_information(g, t, mask))
        out[m] = np.mean(vals)
    return out

def redundancy(g, t, delta=0.1, n_samples=200, rng=None):
    """R_delta = N / m_delta, with m_delta smallest m s.t. avg I >= (1-delta) S_S."""
    N = len(g)
    cos_all = np.cos(2 * g * t)
    S_S = ent_rank2(np.prod(cos_all))
    if S_S < 1e-9:
        return 1.0  # nothing to record yet
    pip = partial_info_plot(g, t, n_samples, rng)
    target = (1 - delta) * S_S
    hits = np.where(pip[1:] >= target)[0]
    if len(hits) == 0:
        return 1.0
    m_delta = hits[0] + 1
    return N / m_delta

def sample_polarization(g, t, n_traj=4000, rng=None):
    """Branch-structure sampling: P(t) = E|2q-1| and the <I>_xi/I0 sample set.

    Snapshot protocol: all bath qubits measured in sigma_y at time t.
    P(r|s) = (1 + r s eps_n)/2 with eps_n = sin(2 g_n t).
    Posterior q = P(s=+1 | r_1..r_N).
    """
    rng = rng or np.random.default_rng(1)
    eps = np.sin(2 * g * t)  # signed; formula handles sign consistently
    s = rng.choice([1, -1], size=n_traj)
    u = rng.random((n_traj, len(g)))
    p_plus = (1 + s[:, None] * eps[None, :]) / 2  # P(r=+1|s)
    r = np.where(u < p_plus, 1, -1)
    # log-odds of s=+1 vs s=-1 given record
    lo = np.sum(np.log1p(r * eps[None, :]) - np.log1p(-r * eps[None, :]), axis=1)
    q = 1 / (1 + np.exp(-lo))
    pol = 2 * q - 1  # = <I>_xi / I0
    return np.mean(np.abs(pol)), pol

def record_distinguishability(g, t):
    """D(t) = sum_n D_KL(p(r|s) || p(r|-s)) in nats (note Eq. 7, snapshot version)."""
    eps = np.abs(np.sin(2 * g * t))
    eps = np.clip(eps, 0, 1 - 1e-12)
    return np.sum(eps * np.log((1 + eps) / (1 - eps)))
