"""Collision protocol (note §2): qubit n interacts during [t_{n-1}, t_n], then
is measured. Records are permanent; all quantities are functions of the
collision index k (time t = k*dt).

eps_n = sin(2 g_n dt): record quality of collision n.
c_n   = cos(2 g_n dt): per-collision decoherence factor.
"""
import numpy as np
from exact import ent_rank2, H2

def collision_S_entropy(g, dt, k):
    """S(rho_S) after k collisions."""
    return ent_rank2(np.prod(np.cos(2 * g[:k] * dt)))

def collision_D(g, dt, k):
    """Accumulated record distinguishability = <sigma_s> (note Eq. 7), nats."""
    eps = np.clip(np.abs(np.sin(2 * g[:k] * dt)), 0, 1 - 1e-12)
    return np.sum(eps * np.log((1 + eps) / (1 - eps)))

def collision_sample_sigma_pol(g, dt, k, n_traj, rng):
    """Sample forward trajectories of the first k records within a branch.

    Returns (sigma samples, polarization samples, IFT estimator <e^-sigma>).
    sigma_s[xi] = sum_n ln[(1 + r_n s eps_n)/(1 - r_n s eps_n)]   (note Eq. 6)
    """
    eps = np.sin(2 * g[:k] * dt)
    s = rng.choice([1, -1], size=n_traj)
    u = rng.random((n_traj, k))
    p_plus = (1 + s[:, None] * eps[None, :]) / 2
    r = np.where(u < p_plus, 1, -1)
    x = r * s[:, None] * eps[None, :]
    sigma = np.sum(np.log1p(x) - np.log1p(-x), axis=1)
    # posterior polarization from the same records
    lo = np.sum(np.log1p(r * eps[None, :]) - np.log1p(-r * eps[None, :]), axis=1)
    q = 1 / (1 + np.exp(-lo))
    pol = np.abs(2 * q - 1)
    ift = float(np.mean(np.exp(-sigma)))
    return sigma, pol, ift

def collision_redundancy(g, dt, k, delta=0.1, n_samples=80, rng=None):
    """R_delta over the k collided qubits (uncollided qubits are uncorrelated)."""
    rng = rng or np.random.default_rng(0)
    cos_k = np.cos(2 * g[:k] * dt)
    S_S = ent_rank2(np.prod(cos_k))
    if S_S < 1e-9 or k == 0:
        return 1.0
    target = (1 - delta) * S_S
    for m in range(1, k + 1):
        vals = []
        for _ in range(n_samples):
            idx = rng.choice(k, size=m, replace=False)
            mask = np.zeros(k, bool); mask[idx] = True
            cF = np.prod(cos_k[mask]); cFb = np.prod(cos_k[~mask])
            vals.append(ent_rank2(cF * cFb) + ent_rank2(cF) - ent_rank2(cFb))
        if np.mean(vals) >= target:
            return k / m
    return 1.0
