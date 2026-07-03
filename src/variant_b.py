"""Variant B (spec §7): H_E = sum_n omega_n sigma_x^(n). Records degrade.

Chirality is still conserved => two-branch structure survives:
|Psi(t)> = (1/sqrt2)[ |+> prod|phi_+^n(t)>  +  |-> prod|phi_-^n(t)> ]
with |phi_s^n(t)> = exp(-i(s g_n sz + omega_n sx) t)|+x>.

All bipartite entropies remain rank-2 formulas in the complex overlaps
c_n(t) = <phi_-^n|phi_+^n>. Optimal (Helstrom) per-record quality:
eps_n(t) = sqrt(1 - |c_n(t)|^2).

Microreversibility note: the correct antiunitary is
Theta = K_ring (x) (sigma_x K)^{(x)N}  (flips sz, preserves sx, fixes |+x>).
"""
import numpy as np
from exact import ent_rank2

SX = np.array([[0, 1], [1, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)
PLUS_X = np.array([1, 1], complex) / np.sqrt(2)

def conditional_states(g_n, w_n, t):
    """|phi_+> and |phi_-> for one qubit at time t (analytic Rabi rotation)."""
    out = []
    for s in (+1, -1):
        a = np.array([w_n, 0.0, s * g_n])
        na = np.linalg.norm(a)
        if na < 1e-15:
            out.append(PLUS_X.copy()); continue
        ah = a / na
        U = np.cos(na * t) * np.eye(2) - 1j * np.sin(na * t) * (
            ah[0] * SX + ah[2] * SZ)
        out.append(U @ PLUS_X)
    return out[0], out[1]

def overlaps(g, w, t):
    """Complex overlaps c_n(t) for all N qubits."""
    return np.array([np.vdot(*conditional_states(g[n], w[n], t)[::-1])
                     for n in range(len(g))])
    # note: c_n = <phi_-|phi_+>

def eps_from_c(c):
    return np.sqrt(np.clip(1 - np.abs(c) ** 2, 0, 1))

def S_system(c):
    return ent_rank2(np.prod(c))

def mutual_information_B(c, frag_mask):
    cF = np.prod(c[frag_mask]); cFb = np.prod(c[~frag_mask])
    return ent_rank2(cF * cFb) + ent_rank2(cF) - ent_rank2(cFb)

def redundancy_B(c, delta=0.1, n_samples=80, rng=None):
    rng = rng or np.random.default_rng(0)
    N = len(c)
    S_S = ent_rank2(np.prod(c))
    if S_S < 1e-9:
        return 1.0
    target = (1 - delta) * S_S
    for m in range(1, N + 1):
        vals = []
        for _ in range(n_samples):
            idx = rng.choice(N, size=m, replace=False)
            mask = np.zeros(N, bool); mask[idx] = True
            vals.append(mutual_information_B(c, mask))
        if np.mean(vals) >= target:
            return N / m
    return 1.0

def polarization_D(eps, n_traj=4000, rng=None):
    """P = E|2q-1| and D = sum KL, for given per-record qualities eps_n."""
    rng = rng or np.random.default_rng(1)
    eps = np.clip(eps, 0, 1 - 1e-12)
    s = rng.choice([1, -1], size=n_traj)
    u = rng.random((n_traj, len(eps)))
    p_plus = (1 + s[:, None] * eps[None, :]) / 2
    r = np.where(u < p_plus, 1, -1)
    lo = np.sum(np.log1p(r * eps[None, :]) - np.log1p(-r * eps[None, :]), axis=1)
    q = 1 / (1 + np.exp(-lo))
    P = float(np.mean(np.abs(2 * q - 1)))
    D = float(np.sum(eps * np.log((1 + eps) / (1 - eps))))
    return P, D
