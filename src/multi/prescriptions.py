"""Shared-qubit prescriptions for the C5 frustration scenario (consolidation T1).

Definitions: notes/multi_definitions.md §6bis. model_multi.py (validated) is NOT
modified; this module generalises its frustration sampler to per-qubit charge
durations T_c:

    fresh  (T_c=0):  mu_n(S,t) = sin(2 a_n(S) t),         draw (1+r mu)/2
    charged (T_c>0): mu_n(S,t) = sin(2 a_n(S) (T_c - t)), draw (1-r mu)/2

Prescriptions (shared qubits only; exclusives fixed: ring-1-private fresh,
ring-2-private charged with T_c=T):
    F: shared fresh (T_c=0)      — reference, must reproduce
       model_multi.sample_frustration to machine identity (same rng).
    M: shared charged, T_c = T
    H: shared charged, T_c = T/2

Run self-test:  python3 prescriptions.py   (from src/multi)
"""
import os, sys
import numpy as np
from scipy.special import logsumexp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import model_multi as mm


def masks_from_couplings(g_mat):
    """(shared, ring1_private, ring2_private) boolean masks from g_mat (2,N)."""
    g_mat = np.asarray(g_mat, float)
    c0, c1 = g_mat[0] > 0, g_mat[1] > 0
    return c0 & c1, c0 & ~c1, ~c0 & c1


def charge_profile(prescription, g_mat, T):
    """Per-qubit Theta-charge duration T_c (N,): 0 = fresh."""
    shared, r1p, r2p = masks_from_couplings(g_mat)
    cT = np.zeros(g_mat.shape[1])
    cT[r2p] = T                                  # ring-2-private: always mirror
    if prescription == "F":
        pass                                     # shared fresh
    elif prescription == "M":
        cT[shared] = T
    elif prescription == "H":
        cT[shared] = T / 2.0
    else:
        raise ValueError(prescription)
    return cT


def sample_frustration_rx(g_mat, t, charge_T, n_traj=20_000, rng=None):
    """Frustration sampler with per-qubit charge durations (§6bis conventions).

    Structurally identical to model_multi.sample_frustration (same rng call
    order) so that charge_T = T*mask reproduces it bit-for-bit; the only
    generalisation is mu on charged qubits: sin(2 a (T_c - t)) per qubit.
    Returns {"sigma": (M,), "P_i": (M,), "sigma_S": float}.
    """
    rng = rng or np.random.default_rng(0)
    g_mat = np.asarray(g_mat, float)
    M, N = g_mat.shape
    charge_T = np.asarray(charge_T, float)
    charged = charge_T > 0
    S_all = mm.all_branches(M)
    B = S_all.shape[0]
    a = S_all @ g_mat                                    # (B,N)
    mu = np.empty((B, N))
    mu[:, ~charged] = np.sin(2.0 * a[:, ~charged] * t)   # fresh: forward
    mu[:, charged] = np.sin(2.0 * a[:, charged]
                            * (charge_T[charged][None, :] - t))  # residual
    mu = np.clip(mu, -1 + 1e-12, 1 - 1e-12)

    bidx = rng.integers(0, B, size=n_traj)
    mu_true = mu[bidx]
    sgn = np.where(charged, -1.0, 1.0)[None, :]          # charged: backward draw
    u = rng.random((n_traj, N))
    r = np.where(u < (1 + sgn * mu_true) / 2, 1, -1).astype(float)

    logL = np.log1p(r[:, None, :] * mu[None, :, :]).sum(axis=2)
    post = np.exp(logL - logL.max(axis=1, keepdims=True))
    post /= post.sum(axis=1, keepdims=True)
    m_i = post @ S_all.astype(float)
    P_i = np.mean(np.abs(m_i), axis=0)

    sig_S = np.sum(np.log1p(r * mu_true) - np.log1p(-r * mu_true), axis=1)
    sigma = np.zeros(M)
    for i in range(M):
        plus = S_all[:, i] > 0
        Lp = logsumexp(logL[:, plus], axis=1)
        Lm = logsumexp(logL[:, ~plus], axis=1)
        sigma[i] = float(np.mean(S_all[bidx, i].astype(float) * (Lp - Lm)))
    return {"sigma": sigma, "P_i": P_i, "sigma_S": float(np.mean(sig_S))}


def _selftest():
    """F must reproduce model_multi.sample_frustration bit-for-bit; M/H limits."""
    import run_scenarios as rs
    rng = np.random.default_rng(77)
    N, T, t = 32, 1.5, 0.6
    g, mirror_mask = rs.couplings_shared(N, 0.4, np.random.default_rng(3))
    # F == existing sampler (same seed -> identical draws -> identical floats)
    ref = mm.sample_frustration(g, t, T, mirror_mask, 5000,
                                np.random.default_rng(9))
    new = sample_frustration_rx(g, t, charge_profile("F", g, T), 5000,
                                np.random.default_rng(9))
    ok = (np.array_equal(ref["sigma"], new["sigma"])
          and np.array_equal(ref["P_i"], new["P_i"])
          and ref["sigma_S"] == new["sigma_S"])
    print(f"  F == model_multi.sample_frustration (bit-for-bit) ... "
          f"{'PASS' if ok else 'FAIL'}")

    # M at f=0 has no shared qubits -> identical to F at f=0
    g0, _ = rs.couplings_shared(N, 0.0, np.random.default_rng(3))
    a = sample_frustration_rx(g0, t, charge_profile("F", g0, T), 5000,
                              np.random.default_rng(9))
    b = sample_frustration_rx(g0, t, charge_profile("M", g0, T), 5000,
                              np.random.default_rng(9))
    ok2 = np.array_equal(a["sigma"], b["sigma"])
    print(f"  M == F at f=0 (no shared qubits) ... {'PASS' if ok2 else 'FAIL'}")

    # M at f=1: all qubits charged -> both rings consume: sigma_1 ~ -D_1 < 0
    g1, _ = rs.couplings_shared(N, 1.0, np.random.default_rng(3))
    c = sample_frustration_rx(g1, t, charge_profile("M", g1, T), 100_000,
                              np.random.default_rng(9))
    ok3 = c["sigma"][0] < 0 and c["sigma"][1] < 0
    print(f"  M at f=1: both arrows negative (both consume) "
          f"sigma={np.round(c['sigma'],2)} ... {'PASS' if ok3 else 'FAIL'}")
    return ok and ok2 and ok3


if __name__ == "__main__":
    print("=== prescriptions.py self-test ===")
    sys.exit(0 if _selftest() else 1)
