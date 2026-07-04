"""Statistical stability of the frustration threshold f* (consolidation T2).

Prescription F only. f* is extracted PER disorder realization (zero crossing of
sigma_2(f), linearly interpolated), for N in {32, 64, 128}, at t in {0.5, 1.0},
T=1.5, 20 realizations, n_traj=20k.

Criterion (reported for results_multi.md): if std(f*) > 0.1 the wording must be
"crossover", not "threshold".

Output: fig_multi_fstar_stability.png + printed table.

Run:  python3 run_fstar_stability.py    (from src/multi)
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import run_scenarios as rs
import prescriptions as px
from run_c5_robustness import zero_crossings

FIG = rs.FIG
T = 1.5
R = 20
NTRAJ = 20_000
FS = np.linspace(0.0, 1.0, 11)
NS = [32, 64, 128]
TS = [0.5, 1.0]


def fstar_one(N, t, seed):
    """f* for one disorder realization (nan if no crossing in [0,1])."""
    y = np.empty(len(FS))
    for j, f in enumerate(FS):
        g, _ = rs.couplings_shared(N, f, np.random.default_rng(seed))
        d = px.sample_frustration_rx(g, t, px.charge_profile("F", g, T),
                                     NTRAJ, np.random.default_rng(99))
        y[j] = d["sigma"][1]
    cr = zero_crossings(FS, y)
    return cr[0] if cr else np.nan


def main():
    print(f"=== f* stability (prescription F, T={T}, R={R}, n_traj={NTRAJ}) ===")
    stats = {}
    for t in TS:
        for N in NS:
            vals = np.array([fstar_one(N, t, 6000 + r) for r in range(R)])
            ok = ~np.isnan(vals)
            m, s = vals[ok].mean(), vals[ok].std()
            stats[(t, N)] = (m, s, int(ok.sum()))
            word = "crossover" if s > 0.1 else "threshold"
            print(f"  t={t:.1f} N={N:4d}: f* = {m:.3f} +/- {s:.3f}  "
                  f"({ok.sum()}/{R} realizations w/ crossing)  -> {word}")

    plt.figure(figsize=(6.5, 4.4))
    for t, mk in zip(TS, ("o-", "s--")):
        m = [stats[(t, N)][0] for N in NS]
        s = [stats[(t, N)][1] for N in NS]
        plt.errorbar(NS, m, yerr=s, fmt=mk, capsize=4, label=f"t={t}")
    plt.axhspan(0, 0, color="k")   # keep y range honest around data
    plt.xscale("log")
    plt.xticks(NS, [str(n) for n in NS])
    plt.xlabel("N (bath size)")
    plt.ylabel(r"$f^{*}$ (zero of $\sigma_2$, per-realization)")
    plt.title(f"Frustration crossing point: disorder mean ± std (R={R})")
    plt.legend(); plt.tight_layout()
    plt.savefig(FIG + "fig_multi_fstar_stability.png", dpi=150)
    plt.close()
    print("figure -> fig_multi_fstar_stability.png")


if __name__ == "__main__":
    main()
