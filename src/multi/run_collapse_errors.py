"""Collapse RMS with disorder error bars (consolidation T4).

Redo the step-5 collapse per disorder realization (R=10): for each realization
build the generic master curve (shared f in {0,0.25,0.5,0.75,1} + random) from
that realization's points only, compute RMS_generic and RMS_DFS against it.
Report mean ± std of both and the worst-case separation.

Output: fig_multi_collapse_errors.png + printed numbers.
Run:  python3 run_collapse_errors.py    (from src/multi)
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import model_multi as mm
import run_scenarios as rs
from collapse_multi import master_curve

FIG = rs.FIG
N = 48
NTRAJ = 20_000
R = 10
TS = np.linspace(0.1, 2.5, 14)
GEN_CONFIGS = [("shared", f) for f in (0.0, 0.25, 0.5, 0.75, 1.0)] + \
              [("random", None)]


def points_one(kind, f, rng_geom, t_list):
    """(D_i, P_i) points for one coupling draw across times and rings."""
    if kind == "shared":
        g, _ = rs.couplings_shared(N, f, rng_geom)
    elif kind == "random":
        g = rs.couplings_random_shared(N, rng_geom)
    else:                                   # identical -> DFS
        g = rs.couplings_identical(N, rng_geom)
    D, P = [], []
    for t in t_list:
        d = mm.sample_joint(g, t, NTRAJ, np.random.default_rng(1234))
        for i in range(2):
            D.append(d["sigma_i"][i]); P.append(d["P_i"][i])
    return np.array(D), np.array(P)


def main():
    print(f"=== T4 collapse RMS with disorder error bars (R={R}) ===")
    rms_gen, rms_dfs = [], []
    for r in range(R):
        Dg, Pg = [], []
        for kind, f in GEN_CONFIGS:
            D, P = points_one(kind, f, np.random.default_rng(9000 + r), TS)
            Dg.append(D); Pg.append(P)
        Dg, Pg = np.concatenate(Dg), np.concatenate(Pg)
        keep = Dg > 0
        Dg, Pg = Dg[keep], Pg[keep]
        mD, mP = master_curve(Dg, Pg, nb=16)
        rms_gen.append(float(np.sqrt(np.mean(
            (Pg - np.interp(Dg, mD, mP)) ** 2))))
        Dd, Pd = points_one("identical", None, np.random.default_rng(9000 + r),
                            TS)
        kd = Dd > 0
        rms_dfs.append(float(np.sqrt(np.mean(
            (Pd[kd] - np.interp(Dd[kd], mD, mP)) ** 2))))

    rms_gen, rms_dfs = np.array(rms_gen), np.array(rms_dfs)
    print(f"  generic RMS = {rms_gen.mean():.4f} +/- {rms_gen.std():.4f}  "
          f"(max over realizations {rms_gen.max():.4f})")
    print(f"  DFS     RMS = {rms_dfs.mean():.4f} +/- {rms_dfs.std():.4f}  "
          f"(min over realizations {rms_dfs.min():.4f})")
    sep = rms_dfs.min() / rms_gen.max()
    print(f"  worst-case separation: min(DFS)/max(generic) = {sep:.1f}x  "
          f"-> {'separation SURVIVES disorder' if sep > 3 else 'MARGINAL'}")

    plt.figure(figsize=(6, 4.4))
    x = np.arange(R)
    plt.plot(x, rms_gen, "o-", label="generic (shared+random)")
    plt.plot(x, rms_dfs, "s-", label="identical (DFS)")
    plt.axhline(rms_gen.mean(), color="C0", ls=":", lw=1)
    plt.axhline(rms_dfs.mean(), color="C1", ls=":", lw=1)
    plt.yscale("log")
    plt.xlabel("disorder realization"); plt.ylabel("collapse RMS residual")
    plt.title(f"Collapse RMS per realization: generic "
              f"{rms_gen.mean():.3f}±{rms_gen.std():.3f} vs DFS "
              f"{rms_dfs.mean():.3f}±{rms_dfs.std():.3f}")
    plt.legend(); plt.tight_layout()
    plt.savefig(FIG + "fig_multi_collapse_errors.png", dpi=150)
    plt.close()
    print("figure -> fig_multi_collapse_errors.png")


if __name__ == "__main__":
    main()
