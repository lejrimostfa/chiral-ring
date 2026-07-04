"""Multi-ring chiral model — closed forms (Paper 2).

Foundation: Paper 1 modules in ../ (exact.py, variant_b.py) are NOT modified.
All quantities here are closed forms derived from the single-qubit conditional
state |phi_n(S,t)> = e^{-i a_n(S) t sigma_z}|+x>, a_n(S) = sum_i s_i g_{i,n}.
Only posterior *sampling* (records + Bayes) is Monte-Carlo.

Definitions and equation numbers: notes/multi_definitions.md.

Run tests:  python3 model_multi.py    (from src/multi)
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import exact  # Paper-1 closed forms, used as the M=1 ground truth


# ----------------------------------------------------------------------
# Branch bookkeeping
# ----------------------------------------------------------------------
def all_branches(M):
    """All 2^M chirality tuples S in {+1,-1}^M, shape (2^M, M).

    Row 0 is all +1; the last is all -1 (antipode of row 0). Ordering is fixed
    and used throughout (branch index b <-> row)."""
    bits = ((np.arange(2 ** M)[:, None] >> np.arange(M)[::-1][None, :]) & 1)
    return 1 - 2 * bits  # 0 -> +1, 1 -> -1


def branch_rates(g_mat, S_all):
    """a_n(S) = sum_i s_i g_{i,n}  (Eq. 1). g_mat (M,N), S_all (B,M) -> (B,N)."""
    return np.asarray(S_all) @ np.asarray(g_mat)


# ----------------------------------------------------------------------
# Coherences: Gamma, rho_joint, DFS  (Eqs. 2-4)
# ----------------------------------------------------------------------
def coherence(g_mat, S, Sp, t):
    """Gamma_{S,S'}(t) = prod_n cos((a_n(S)-a_n(S')) t)  (Eqs. 2-3)."""
    g_mat = np.asarray(g_mat, float)
    aS = np.asarray(S, float) @ g_mat
    aSp = np.asarray(Sp, float) @ g_mat
    return float(np.prod(np.cos((aS - aSp) * t)))


def rho_joint(g_mat, t):
    """Joint ring reduced state (2^M x 2^M) in the chiral basis (Eq. 3)."""
    M = np.asarray(g_mat).shape[0]
    S_all = all_branches(M)
    a = branch_rates(g_mat, S_all)  # (B,N)
    B = a.shape[0]
    rho = np.empty((B, B), float)
    for i in range(B):
        for j in range(B):
            rho[i, j] = np.prod(np.cos((a[i] - a[j]) * t))
    return rho / B


def is_dfs_pair(g_mat, S, Sp):
    """True iff branches S,S' are mutually decoherence-free: a_n(S)=a_n(S') all n
    (Eq. 4). Exact, t-independent."""
    g_mat = np.asarray(g_mat, float)
    return bool(np.allclose((np.asarray(S, float) - np.asarray(Sp, float)) @ g_mat,
                            0.0, atol=1e-12))


# ----------------------------------------------------------------------
# Records: bias, distinguishability  (Eqs. 5-8)
# ----------------------------------------------------------------------
def record_bias(g_mat, S, t):
    """mu_n(S,t) = sin(2 a_n(S) t)  (Eq. 5), signed, shape (N,)."""
    a = np.asarray(S, float) @ np.asarray(g_mat, float)
    return np.sin(2.0 * a * t)


def distinguishability(g_mat, S, t):
    """D(S,t) = sum_n eps_n ln((1+eps_n)/(1-eps_n)), eps_n=|mu_n(S,t)|  (Eq. 8).

    Branch-antipode KL distinguishability = <sigma_S>. For M=1 this equals
    exact.record_distinguishability."""
    eps = np.clip(np.abs(record_bias(g_mat, S, t)), 0, 1 - 1e-12)
    return float(np.sum(eps * np.log((1 + eps) / (1 - eps))))


def relative_distinguishability(g_mat, S, Sp, t):
    """D_rel(t) between two branches S,S' (Eq. 20): KL of their record laws.
    eps_n = |sin((a_n(S)-a_n(S')) t)| (half-angle of the collision bias diff)."""
    g_mat = np.asarray(g_mat, float)
    d = (np.asarray(S, float) - np.asarray(Sp, float)) @ g_mat
    eps = np.clip(np.abs(np.sin(d * t)), 0, 1 - 1e-12)
    return float(np.sum(eps * np.log((1 + eps) / (1 - eps))))


# ----------------------------------------------------------------------
# Concurrence for the inter-ring entangled state (Eq. 18-19), M=2
# ----------------------------------------------------------------------
def concurrence_pm(g_mat, t):
    """Concurrence C(t) of rho_12 for initial |Psi>=(|+->+|-+>)/sqrt2 (M=2).

    This state lives in span{(+,-),(-,+)}; under pure dephasing rho_12 stays an
    X-state whose only dynamical off-diagonal is Gamma between those two
    branches, so C(t) = |Gamma_{(+,-),(-,+)}(t)| (Eq. 19)."""
    S_pm, S_mp = np.array([1, -1]), np.array([-1, 1])
    return abs(coherence(g_mat, S_pm, S_mp, t))


# ----------------------------------------------------------------------
# Posterior sampling (the ONLY Monte-Carlo piece)  (Eqs. 9-13)
# ----------------------------------------------------------------------
def sample_joint(g_mat, t, n_traj=20_000, rng=None):
    """Forward collision sampling of records + Bayes posterior over 2^M branches.

    Draws S ~ uniform, records r_n ~ P(r|S)=(1+r mu_n(S))/2, forms the posterior
    P(S|xi) (Eq. 9) and returns per-ring/joint arrows and correlations.

    Returns dict:
      P_i        : polarization per ring, E_xi|E[s_i|xi]|            (Eq. 10)
      chi        : <E[s_1 s_2|xi]> ring-ring record correlation (M=2, Eq. 11)
      sigma_S    : <sigma_S[xi]> joint branch-relative EP           (Eq. 7)
      sigma_i    : <sigma_i[xi]> marginal EP per ring, Option A     (Eq. 13)
      deficit    : sigma_S - sum_i sigma_i                          (Eq. 14)
      D          : closed-form sum_S P(S) D(S,t) (matches sigma_S)  (Eq. 8)
      pol_total  : polarization of the TOTAL chirality sum_i s_i (for DFS/C1)
    """
    rng = rng or np.random.default_rng(0)
    g_mat = np.asarray(g_mat, float)
    M, N = g_mat.shape
    S_all = all_branches(M)                 # (B,M)
    B = S_all.shape[0]
    mu = np.sin(2.0 * (S_all @ g_mat) * t)  # (B,N) signed bias per branch
    mu = np.clip(mu, -1 + 1e-12, 1 - 1e-12)

    bidx = rng.integers(0, B, size=n_traj)  # true branch per trajectory
    mu_true = mu[bidx]                       # (n_traj,N)
    u = rng.random((n_traj, N))
    r = np.where(u < (1 + mu_true) / 2, 1, -1).astype(float)  # records (n_traj,N)

    # log-likelihood of each branch: sum_n log(1 + r mu_b)  (prior uniform)
    logL = np.log1p(r[:, None, :] * mu[None, :, :]).sum(axis=2)  # (n_traj,B)
    logL -= logL.max(axis=1, keepdims=True)
    post = np.exp(logL)
    post /= post.sum(axis=1, keepdims=True)  # P(S|xi)  (n_traj,B)

    s_i = S_all.astype(float)                 # (B,M)
    m_i = post @ s_i                          # E[s_i|xi]  (n_traj,M)
    P_i = np.mean(np.abs(m_i), axis=0)        # (M,)

    # total-chirality polarization (objectivity of the sum, C1)
    tot = S_all.sum(axis=1).astype(float)     # (B,)
    m_tot = post @ tot
    denom = abs(tot).max() if abs(tot).max() > 0 else 1.0
    pol_total = float(np.mean(np.abs(m_tot)) / denom)

    out = {"P_i": P_i, "pol_total": pol_total}
    if M == 2:
        s1s2 = (S_all[:, 0] * S_all[:, 1]).astype(float)
        chi_traj = post @ s1s2                       # E[s1 s2|xi] per trajectory
        out["chi"] = float(np.mean(chi_traj))        # ~0 always (tower property)
        # record-induced posterior covariance: 0 for disjoint baths, >0 shared
        out["cov"] = float(np.mean(chi_traj - m_i[:, 0] * m_i[:, 1]))

    # joint branch-relative EP: sigma_S = sum_n log((1+r mu_true)/(1-r mu_true))
    sig_S = np.sum(np.log1p(r * mu_true) - np.log1p(-r * mu_true), axis=1)
    out["sigma_S"] = float(np.mean(sig_S))
    out["sigma_S_traj"] = sig_S          # kept for C4 histograms
    out["D"] = float(np.mean([distinguishability(g_mat, S_all[b], t)
                              for b in bidx]))

    # marginal EP (Option A / Eq. 13): sigma_i = s_i * log(L_i^+ / L_i^-)
    from scipy.special import logsumexp
    sigma_i = np.zeros(M)
    logL_raw = np.log1p(r[:, None, :] * mu[None, :, :]).sum(axis=2)  # unnormalised
    for i in range(M):
        plus = S_all[:, i] > 0
        Lp = logsumexp(logL_raw[:, plus], axis=1)
        Lm = logsumexp(logL_raw[:, ~plus], axis=1)
        s_i_true = S_all[bidx, i].astype(float)
        sigma_i[i] = float(np.mean(s_i_true * (Lp - Lm)))
    out["sigma_i"] = sigma_i
    out["deficit"] = out["sigma_S"] - float(sigma_i.sum())
    return out


# ----------------------------------------------------------------------
# Mirror preparation (Eqs. 16-17)
# ----------------------------------------------------------------------
def mirror_eps(g_row, T, t):
    """eps_n^mir(t) = |sin(2 g_n (T-t))| (Eq. 16), single ring, shape (N,)."""
    return np.abs(np.sin(2.0 * np.asarray(g_row, float) * (T - t)))


def mirror_distinguishability(g_row, T, t):
    """D^mir(t) = sum_n eps^mir ln((1+eps^mir)/(1-eps^mir)) (magnitude of Eq.17)."""
    eps = np.clip(mirror_eps(g_row, T, t), 0, 1 - 1e-12)
    return float(np.sum(eps * np.log((1 + eps) / (1 - eps))))


def sample_mirror_sigma(g_row, T, t, n_traj=20_000, rng=None, s=+1):
    """Draw records from the BACKWARD law and score sigma with the FORWARD law
    (Eq. 17). <sigma> -> -D^mir(t). Single ring, branch s (default +1).

    Returns (mean, sem): sem is the standard error of the mean, which grows like
    sqrt(N) since sigma is a sum over qubits — the correct scale for comparison
    against the closed form -D^mir."""
    rng = rng or np.random.default_rng(0)
    b = s * np.sin(2.0 * np.asarray(g_row, float) * (T - t))  # signed forward bias
    b = np.clip(b, -1 + 1e-12, 1 - 1e-12)
    # backward law: P_bwd(r|s) = (1 - r b)/2
    u = rng.random((n_traj, len(b)))
    r = np.where(u < (1 - b) / 2, 1, -1).astype(float)
    sig = np.sum(np.log1p(r * b) - np.log1p(-r * b), axis=1)  # forward score
    return float(np.mean(sig)), float(np.std(sig) / np.sqrt(n_traj))


# ----------------------------------------------------------------------
# Unit tests (work-order step 2: a, b, c)
# ----------------------------------------------------------------------
def _test_a_reduce_to_M1():
    """(a) M=1 reproduces exact.py to machine precision: Gamma, mu, D."""
    rng = np.random.default_rng(1)
    g = rng.uniform(0.05, 0.4, size=10)
    g_mat = g[None, :]              # M=1
    Splus, Sminus = np.array([1]), np.array([-1])
    ok = True
    for t in [0.3, 1.1, 2.7]:
        # Gamma_{+,-} vs prod cos(2 g t)
        G = coherence(g_mat, Splus, Sminus, t)
        G_ref = float(np.prod(np.cos(2 * g * t)))
        # bias vs sin(2 g t)
        mu = record_bias(g_mat, Splus, t)
        mu_ref = np.sin(2 * g * t)
        # D vs exact.record_distinguishability
        D = distinguishability(g_mat, Splus, t)
        D_ref = exact.record_distinguishability(g, t)
        ok &= abs(G - G_ref) < 1e-13
        ok &= np.max(np.abs(mu - mu_ref)) < 1e-13
        ok &= abs(D - D_ref) < 1e-12
    print(f"  (a) M=1 == exact.py (Gamma, mu, D) ... {'PASS' if ok else 'FAIL'}")
    return ok


def _test_b_dfs():
    """(b) Identical couplings => (+,-)<->(-,+) protected, (+,+)<->(-,-) decays."""
    rng = np.random.default_rng(2)
    g = rng.uniform(0.05, 0.4, size=12)
    g_mat = np.vstack([g, g])       # M=2 identical couplings
    S_pm, S_mp = np.array([1, -1]), np.array([-1, 1])
    S_pp, S_mm = np.array([1, 1]), np.array([-1, -1])
    ok = is_dfs_pair(g_mat, S_pm, S_mp)                       # exact DFS
    ok &= abs(relative_distinguishability(g_mat, S_pm, S_mp, 1.7)) < 1e-13
    ok &= abs(concurrence_pm(g_mat, 3.3) - 1.0) < 1e-13       # C(t)=1 protected
    # (+,+)<->(-,-) must decay: pick t where prod cos(4 g t) clearly < 1
    G_dec = coherence(g_mat, S_pp, S_mm, 1.7)
    ok &= (not is_dfs_pair(g_mat, S_pp, S_mm)) and abs(G_dec) < 0.99
    print(f"  (b) DFS: (+,-)~(-,+) protected, (+,+)~(-,-) decays ... "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def _test_c_mirror():
    """(c) Mirror at f=0: sampled <sigma_2>(t) = -D_2^mir(t) (forward correspond.)."""
    rng = np.random.default_rng(3)
    g = rng.uniform(0.05, 0.4, size=40)
    T = 2.0
    ok = True
    for t in [0.3, 0.9, 1.5]:
        s_samp, sem = sample_mirror_sigma(g, T, t, n_traj=200_000,
                                          rng=np.random.default_rng(100 + int(10 * t)))
        D_mir = mirror_distinguishability(g, T, t)
        ok &= abs(s_samp - (-D_mir)) < 5 * sem  # within 5 standard errors
    print(f"  (c) mirror f=0: <sigma>= -D^mir(t) (within 5 SEM) ... "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def run_unit_tests():
    print("=== model_multi.py unit tests ===")
    a = _test_a_reduce_to_M1()
    b = _test_b_dfs()
    c = _test_c_mirror()
    allok = a and b and c
    print(f"ALL: {'PASS' if allok else 'FAIL'}")
    return allok


if __name__ == "__main__":
    ok = run_unit_tests()
    sys.exit(0 if ok else 1)
