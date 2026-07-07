"""C2 re-analysis: what IS the additivity deficit? (exact enumeration)

The first pass compared Delta = <sigma_S> - sum_i <sigma_i> against the
POSTERIOR total correlation (bounded by M ln2) and rejected the identity as a
category error. That comparison was too weak: "multi-information of the
records" admits EXTENSIVE readings. This script tests them exactly.

Exact setup (no sampling): M=2, N small enough to enumerate all 2^N records.
For each record xi and branch S, P(xi|S) = prod_n (1 + xi_n mu_n(S))/2 with
mu_n(S) = sin(2 a_n(S) t). Everything below is an exact sum over xi.

Quantities (all in nats):
  Delta   = D_joint - D_1 - D_2,
            D_joint = E_S KL(P(.|S) || P(.|-S))          <- <sigma_S>
            D_i     = E_{s_i} KL(P(.|s_i) || P(.|-s_i))  <- <sigma_i>, Option A
                      (identity <sigma_i> = marginal KL: proven by
                       marginalisation over the uniform independent s_j)
  S_syn   = sum_n [I(r_n;S) - I(r_n;s_1) - I(r_n;s_2)]   per-qubit cross-ring
            record synergy: extensive, =0 for disjoint baths (closed form)
  TC_rec  = sum_n H(r_n) - H(xi) = N ln2 - H(xi)         full record
            multi-information: extensive but NOT 0 for disjoint baths
            (intra-ring correlations through each s_i)
  I_cond  = E_xi [H(s1|xi)+H(s2|xi)-H(S|xi)]             posterior TC,
            bounded by M ln2 (the first-pass comparator, kept for reference)

Run:  python3 analyze_c2.py    (from src/multi)
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import model_multi as mm

FIG = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "paper2", "")
LN2 = np.log(2.0)


def h_bias(b):
    """Binary entropy (nats) of a +-1 variable with bias b: p(+1)=(1+b)/2."""
    p = np.clip((1.0 + np.asarray(b, float)) / 2.0, 1e-15, 1 - 1e-15)
    return -(p * np.log(p) + (1 - p) * np.log(1 - p))


def exact_quantities(g_mat, t):
    """All C2 candidates by exact enumeration. g_mat (2,N), N <= ~14."""
    M, N = g_mat.shape
    S_all = mm.all_branches(M)                    # rows: (++),(+-),(-+),(--)
    B = S_all.shape[0]
    mu = np.clip(np.sin(2.0 * (S_all @ g_mat) * t), -1 + 1e-12, 1 - 1e-12)

    # all 2^N records as +-1 patterns
    bits = ((np.arange(2 ** N)[:, None] >> np.arange(N)[::-1][None, :]) & 1)
    xi = (1 - 2 * bits).astype(float)             # (2^N, N)

    # P(xi|S): (B, 2^N)
    P = np.prod((1.0 + mu[:, None, :] * xi[None, :, :]) / 2.0, axis=2)
    P = np.clip(P, 1e-300, None)
    anti = B - 1 - np.arange(B)                   # antipode: bit complement

    # D_joint = E_S sum_xi P(xi|S) ln P(xi|S)/P(xi|-S)
    D_joint = float(np.mean(np.sum(P * (np.log(P) - np.log(P[anti])), axis=1)))

    # marginal channels P(xi|s_i)
    D_i = []
    for i in range(M):
        plus = S_all[:, i] > 0
        Pp = P[plus].mean(axis=0)                 # P(xi|s_i=+1)
        Pm = P[~plus].mean(axis=0)
        D_i.append(0.5 * float(np.sum(Pp * (np.log(Pp) - np.log(Pm))))
                   + 0.5 * float(np.sum(Pm * (np.log(Pm) - np.log(Pp)))))
    Delta = D_joint - sum(D_i)

    # per-qubit cross-ring synergy, closed form (no enumeration needed)
    Ebar1 = np.array([mu[S_all[:, 0] == s].mean(axis=0) for s in (1, -1)])
    Ebar2 = np.array([mu[S_all[:, 1] == s].mean(axis=0) for s in (1, -1)])
    I_rS = LN2 - h_bias(mu).mean(axis=0)          # (N,)  I(r_n;S)
    I_r1 = LN2 - h_bias(Ebar1).mean(axis=0)       # I(r_n;s_1)
    I_r2 = LN2 - h_bias(Ebar2).mean(axis=0)
    S_syn = float(np.sum(I_rS - I_r1 - I_r2))

    # full record multi-information
    Pxi = P.mean(axis=0)
    H_xi = float(-np.sum(Pxi * np.log(Pxi)))
    TC_rec = N * LN2 - H_xi

    # posterior TC (bounded comparator from the first pass)
    post = (P / B) / Pxi[None, :]                 # (B, 2^N) P(S|xi)
    post = np.clip(post, 1e-15, 1.0)
    H_S_xi = -np.sum(post * np.log(post), axis=0)
    Hm = np.zeros(2 ** N)
    for i in range(M):
        pi = np.clip(post[S_all[:, i] > 0].sum(axis=0), 1e-15, 1 - 1e-15)
        Hm += -(pi * np.log(pi) + (1 - pi) * np.log(1 - pi))
    I_cond = float(np.sum(Pxi * (Hm - H_S_xi)))

    # --- exact decomposition:  Delta = I_cond + C_bwd ------------------
    # Derivation (antipodal symmetry P(xi|-S) = P(-xi|S), uniform prior):
    #   Delta = E_fwd[i(xi;s1;s2)] - E_bwd[i(xi;s1;s2)]
    # with the pointwise co-information
    #   i(xi;s1;s2) = ln [ P(xi|S) P(xi) / (P(xi|s1) P(xi|s2)) ],
    # E_fwd over xi ~ P(.|S) (true branch), E_bwd over xi ~ P(.|-S) (antipode),
    # both scored on the true branch S; the -ln P(xi) pieces cancel between the
    # two ensembles because the unconditional mixture is antipode-invariant.
    # Standard identity (s1 independent of s2):
    #   E_fwd[i] = I(xi;S) - I(xi;s1) - I(xi;s2) = I(s1;s2|xi)  (bounded, = I_cond)
    # so the extensive part of Delta is entirely the backward cross-synergy
    #   C_bwd = -E_bwd[i]  (a cross-entropy-type term, NOT a bounded Shannon MI).
    logi = np.log(P) + np.log(Pxi)[None, :]       # ln P(xi|S) + ln P(xi)
    for i in range(M):
        plus = S_all[:, i] > 0
        Pmarg_i = np.empty_like(P)
        Pmarg_i[plus] = P[plus].mean(axis=0)[None, :]     # P(xi|s_i) per branch
        Pmarg_i[~plus] = P[~plus].mean(axis=0)[None, :]
        logi -= np.log(Pmarg_i)
    C_bwd = float(-np.mean(np.sum(P[anti] * logi, axis=1)))
    decomp_err = abs(Delta - (I_cond + C_bwd))

    return {"Delta": Delta, "S_syn": S_syn, "TC_rec": TC_rec,
            "I_cond": I_cond, "C_bwd": C_bwd, "decomp_err": decomp_err,
            "D_joint": D_joint, "D_i": D_i}


def couplings_f(N, f, rng):
    """Same layout as run_scenarios.couplings_shared (shared block first)."""
    n_sh = int(round(f * N)); n_pr = N - n_sh; n1 = n_pr // 2
    g = np.zeros((2, N))
    g[0, :n_sh] = rng.uniform(0.1, 0.4, n_sh)
    g[1, :n_sh] = rng.uniform(0.1, 0.4, n_sh)
    g[0, n_sh:n_sh + n1] = rng.uniform(0.1, 0.4, n1)
    g[1, n_sh + n1:] = rng.uniform(0.1, 0.4, n_pr - n1)
    return g


def main():
    N = 10
    ts = np.geomspace(0.02, 2.0, 25)
    print("=== C2 exact re-analysis (enumeration, N=%d) ===" % N)

    # --- disjoint control: Delta must be 0, S_syn must be 0, TC_rec is NOT 0
    g0 = couplings_f(N, 0.0, np.random.default_rng(11))
    q = exact_quantities(g0, 0.8)
    print(f"  f=0.0 t=0.8: Delta={q['Delta']:.2e}  S_syn={q['S_syn']:.2e}  "
          f"TC_rec={q['TC_rec']:.3f}  I_cond={q['I_cond']:.2e}")
    print("    -> TC_rec (full record multi-info) is nonzero for DISJOINT baths:")
    print("       intra-ring record correlations; cannot equal Delta. Excluded.")

    # --- fully shared: candidate comparison + exact decomposition check
    for seed in (11, 12, 13):
        g1 = couplings_f(N, 1.0, np.random.default_rng(seed))
        rows = [exact_quantities(g1, t) for t in ts]
        Delta = np.array([r["Delta"] for r in rows])
        Ssyn = np.array([r["S_syn"] for r in rows])
        derr = max(r["decomp_err"] for r in rows)
        ratio = Delta / np.where(Ssyn > 1e-12, Ssyn, np.nan)
        i_small = ts < 0.2
        print(f"  f=1.0 seed={seed}: Delta/S_syn small-t "
              f"{np.nanmean(ratio[i_small]):.3f} +/- {np.nanstd(ratio[i_small]):.3f}"
              f" | mid-t (t~0.5) {ratio[np.argmin(abs(ts-0.5))]:.3f}"
              f" | large-t (t={ts[-1]:.1f}) {ratio[-1]:.3f}"
              f" | decomp max|err|={derr:.2e}")

    # --- the reference curve figure (seed 11)
    g1 = couplings_f(N, 1.0, np.random.default_rng(11))
    rows = [exact_quantities(g1, t) for t in ts]
    Delta = np.array([r["Delta"] for r in rows])
    Ssyn = np.array([r["S_syn"] for r in rows])
    TCr = np.array([r["TC_rec"] for r in rows])
    Icond = np.array([r["I_cond"] for r in rows])
    Cbwd = np.array([r["C_bwd"] for r in rows])

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].loglog(ts, Delta, "o-", label=r"deficit $\Delta$")
    ax[0].loglog(ts, np.maximum(Icond, 1e-12), "v-",
                 label=r"$I(s_1;s_2|\xi)$ (bounded part)")
    ax[0].loglog(ts, np.maximum(Cbwd, 1e-12), "d-",
                 label=r"$C_{\rm bwd}$ (extensive part)")
    ax[0].loglog(ts, Ssyn, "s-", label=r"per-qubit synergy $S_{\rm syn}$")
    ax[0].loglog(ts, TCr, "^-", alpha=0.5, label=r"$TC_{\rm rec}$ (excluded)")
    ax[0].axhline(LN2, color="k", ls=":", lw=0.8, label=r"$\ln 2$")
    ax[0].set_xlabel("t"); ax[0].set_ylabel("nats"); ax[0].legend(fontsize=7)
    ax[0].set_title(f"exact decomposition, f=1, N={N}")
    ax[1].semilogx(ts, Icond / Delta, "v-", label=r"$I(s_1;s_2|\xi)/\Delta$")
    ax[1].semilogx(ts, Cbwd / Delta, "d-", label=r"$C_{\rm bwd}/\Delta$")
    ax[1].axhline(1.0, color="k", ls=":", lw=0.8)
    ax[1].set_xlabel("t"); ax[1].set_ylabel("fraction of $\\Delta$")
    ax[1].legend(fontsize=8)
    ax[1].set_title(r"$\Delta = I(s_1;s_2|\xi) + C_{\rm bwd}$ (exact)")
    fig.suptitle("(C2 re-analysis) the deficit decomposed — exact enumeration")
    fig.tight_layout(); fig.savefig(FIG + "fig_multi_c2_exact.png", dpi=150)
    plt.close(fig)

    # --- sweep f at fixed t: Delta vs S_syn turn on together?
    print("  sweep f at t=0.5 (exact):")
    for f in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
        g = couplings_f(N, f, np.random.default_rng(11))
        q = exact_quantities(g, 0.5)
        r = q["Delta"] / q["S_syn"] if q["S_syn"] > 1e-12 else np.nan
        print(f"    f={f:.1f}: Delta={q['Delta']:8.4f}  S_syn={q['S_syn']:8.4f}"
              f"  ratio={r:6.3f}")

    # --- N scaling at fixed t (both should be extensive with same slope)
    print("  N scaling at t=0.5, f=1 (exact):")
    for Ns in [6, 8, 10, 12]:
        g = couplings_f(Ns, 1.0, np.random.default_rng(11))
        q = exact_quantities(g, 0.5)
        print(f"    N={Ns:2d}: Delta={q['Delta']:8.4f}  S_syn={q['S_syn']:8.4f}"
              f"  ratio={q['Delta']/q['S_syn']:6.3f}")


if __name__ == "__main__":
    main()
