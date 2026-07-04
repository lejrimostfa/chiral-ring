"""Deficit scaling in the number of shared records (consolidation T3, C2).

Reinterpretation under test: the additivity deficit Delta = sigma_S - sum_i
sigma_i is a SYNERGISTIC irreversibility, extensive in the number k of shared
qubits.

Parts:
  A. N_tot=64 fixed, k in {0,8,...,64} shared, t in {0.5, 1.0}, R=10,
     n_traj=20k: fit Delta = alpha*k + beta -> alpha, beta, R^2.
  B. N in {16,32,64,128} fully shared (k=N), t=1.0: same fit in N.
  C. Per-record decomposition test (EXACT, enumeration N=10):
        delta_n = KL_n^joint - sum_i KL_n^(i),
        KL_n^joint = E_S[ mu_n(S) L(mu_n(S)) ],        L(x)=ln((1+x)/(1-x))
        KL_n^(i)   = E_{s_i}[ Ebar_n(s_i) L(Ebar_n(s_i)) ],
                     Ebar_n(s_i) = E_{s_j}[ mu_n(S) ]
     delta_n = 0 exactly on private qubits (closed form), so sum_n delta_n is
     manifestly extensive in the shared records. Exactness of
     Delta = sum_n delta_n is tested by enumeration; the residual, if any, is
     the intra-marginal record-correlation term sum_i [sum_n KL_n^(i) - D_i].

Output: fig_multi_deficit_scaling.png + printed fits.
Run:  python3 run_deficit_scaling.py    (from src/multi)
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import model_multi as mm
import run_scenarios as rs
import analyze_c2 as ac

FIG = rs.FIG
NTRAJ = 20_000
R = 10


def sum_delta_n(g_mat, t):
    """Closed-form per-record synergistic deficit sum_n delta_n (nats)."""
    g_mat = np.asarray(g_mat, float)
    M, N = g_mat.shape
    S_all = mm.all_branches(M)
    mu = np.clip(np.sin(2.0 * (S_all @ g_mat) * t), -1 + 1e-12, 1 - 1e-12)

    def xLx(x):
        x = np.clip(x, -1 + 1e-12, 1 - 1e-12)
        return x * np.log((1 + x) / (1 - x))

    kl_joint = xLx(mu).mean(axis=0)                    # (N,) E_S per record
    out = kl_joint.copy()
    for i in range(M):
        Ebar = np.array([mu[S_all[:, i] == s].mean(axis=0) for s in (1, -1)])
        out -= xLx(Ebar).mean(axis=0)                  # E_{s_i}
    return float(out.sum()), out                       # total, per-record


def linfit(x, y):
    """Least-squares y = a x + b with R^2."""
    A = np.vstack([x, np.ones_like(x)]).T
    (a, b), res, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = a * x + b
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return float(a), float(b), 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0


def main():
    print("=== T3 deficit scaling: Delta vs number of shared records ===")

    # ---------- C. exact per-record decomposition (enumeration, N=10) -----
    print("  (C) exact decomposition test, N=10 enumeration:")
    worst = 0.0
    for f in (0.3, 0.6, 1.0):
        for t in (0.2, 0.5, 1.0):
            g = ac.couplings_f(10, f, np.random.default_rng(11))
            q = ac.exact_quantities(g, t)
            sd, _ = sum_delta_n(g, t)
            rel = abs(q["Delta"] - sd) / max(q["Delta"], 1e-12)
            worst = max(worst, rel)
            print(f"    f={f:.1f} t={t:.1f}: Delta={q['Delta']:8.4f}  "
                  f"sum delta_n={sd:8.4f}  rel.err={rel:.3f}")
    exact_holds = worst < 1e-10
    print(f"    -> per-record decomposition exact: {exact_holds} "
          f"(worst rel. residual {worst:.3f})")

    # ---------- A. k sweep at N=64 ----------------------------------------
    N = 64
    ks = np.arange(0, N + 1, 8)
    seeds = [4000 + r for r in range(R)]
    panelA = {}
    for t in (0.5, 1.0):
        dm, ds, sd_cf = [], [], []
        for k in ks:
            vals, cf = [], []
            for s in seeds:
                g, _ = rs.couplings_shared(N, k / N, np.random.default_rng(s))
                d = mm.sample_joint(g, t, NTRAJ, np.random.default_rng(77))
                vals.append(d["deficit"])
                cf.append(sum_delta_n(g, t)[0])
            dm.append(np.mean(vals)); ds.append(np.std(vals))
            sd_cf.append(np.mean(cf))
        dm, ds, sd_cf = map(np.array, (dm, ds, sd_cf))
        a, b, r2 = linfit(ks.astype(float), dm)
        panelA[t] = (dm, ds, sd_cf, a, b, r2)
        print(f"  (A) N=64 t={t}: Delta = {a:.4f}*k + {b:+.4f}   R^2={r2:.5f}")

    # ---------- B. N sweep at f=1 ------------------------------------------
    Ns = np.array([16, 32, 64, 128])
    t_B = 1.0
    dmB, dsB = [], []
    for Nv in Ns:
        vals = []
        for s in seeds:
            g = rs.couplings_random_shared(Nv, np.random.default_rng(s))
            d = mm.sample_joint(g, t_B, NTRAJ, np.random.default_rng(77))
            vals.append(d["deficit"])
        dmB.append(np.mean(vals)); dsB.append(np.std(vals))
    dmB, dsB = np.array(dmB), np.array(dsB)
    aB, bB, r2B = linfit(Ns.astype(float), dmB)
    print(f"  (B) f=1 t={t_B}: Delta = {aB:.4f}*N + {bB:+.4f}   R^2={r2B:.5f}")

    # ---------- figure ------------------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    for t, c in ((0.5, "C0"), (1.0, "C1")):
        dm, ds, sd_cf, a, b, r2 = panelA[t]
        ax[0].errorbar(ks, dm, yerr=ds, fmt="o", color=c, capsize=3,
                       label=f"$\\Delta$ sampled, t={t}")
        ax[0].plot(ks, a * ks + b, "-", color=c, lw=1,
                   label=f"fit {a:.3f}k{b:+.3f} ($R^2$={r2:.4f})")
        ax[0].plot(ks, sd_cf, "--", color=c, lw=1.4, alpha=0.8,
                   label=r"$\sum_n \delta_n$ closed form")
    ax[0].set_xlabel("k shared qubits (N=64)"); ax[0].set_ylabel("nats")
    ax[0].legend(fontsize=7); ax[0].set_title("deficit vs shared count")
    ax[1].errorbar(Ns, dmB, yerr=dsB, fmt="s", color="C3", capsize=3,
                   label=f"$\\Delta$ sampled (f=1, t={t_B})")
    ax[1].plot(Ns, aB * Ns + bB, "-", color="C3", lw=1,
               label=f"fit {aB:.3f}N{bB:+.3f} ($R^2$={r2B:.4f})")
    ax[1].set_xlabel("N (all shared)"); ax[1].set_ylabel("nats")
    ax[1].legend(fontsize=8); ax[1].set_title("deficit vs bath size")
    fig.suptitle("(T3) Synergistic irreversibility: deficit extensive in "
                 "shared records")
    fig.tight_layout()
    fig.savefig(FIG + "fig_multi_deficit_scaling.png", dpi=150)
    plt.close(fig)
    print("figure -> fig_multi_deficit_scaling.png")


if __name__ == "__main__":
    main()
