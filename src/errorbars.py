"""Error-bar production run (paper checklist item 2).

Run from chiral-ring/src:   python3 errorbars.py
Expected runtime on M1 Max: ~20-40 min (all analytic; trajectory sampling
dominates). Outputs: printed table for the paper + three figures with bands.

Statistical design:
- DISORDER: R realizations of couplings g_n. Reported errors = std over
  realizations. This is THE error bar of the paper.
- SAMPLING: n_traj = 50k makes Monte-Carlo noise << disorder noise.
- Hysteresis floor test: residual vs n_traj; 1/sqrt(n) decay => sampling
  floor (report high-n value); plateau => physical deviation (report as is).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import exact, variant_b as vb

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures", "paper1", "")
R = 10                 # disorder realizations (raise to 20 for final)
NTRAJ = 50_000
GBAR = 0.2

# ---------- 1. t_P / t_S ratio with disorder error bars ----------
def timescales_one(N, seed):
    g = np.random.default_rng(seed).uniform(0.5*GBAR, 1.5*GBAR, size=N)
    ts = np.linspace(0, 8.0 * np.sqrt(8.0 / N), 241)   # adaptive window
    S = np.array([exact.ent_rank2(np.prod(np.cos(2*g*t))) for t in ts])
    P = np.array([exact.sample_polarization(g, t, n_traj=NTRAJ,
                  rng=np.random.default_rng(seed+1))[0] for t in ts])
    def cross(y, lvl):
        i = np.argmax(y >= lvl)
        return np.nan if i == 0 else ts[i-1] + (lvl-y[i-1])*(ts[i]-ts[i-1])/(y[i]-y[i-1])
    return cross(S, 0.5), cross(P, 0.5), cross(P, 0.9) - cross(P, 0.1)

print("=== t_P/t_S ratio, disorder-averaged ===")
Ns = [8, 16, 32, 64, 128, 256]
table = {}
for N in Ns:
    vals = np.array([timescales_one(N, 1000*r + N) for r in range(R)])
    tS, tP, w = vals[:, 0], vals[:, 1], vals[:, 2]
    ratio = tP / tS
    table[N] = (tS.mean(), tS.std(), ratio.mean(), ratio.std(),
                (w/tP).mean(), (w/tP).std())
    table[N] = table[N] + (ratio,)
    print(f"N={N:4d}: t_S={tS.mean():.3f}±{tS.std():.3f}  "
          f"t_P/t_S={ratio.mean():.3f}±{ratio.std():.3f}  "
          f"rel.width={(w/tP).mean():.2f}±{(w/tP).std():.2f}")
allr = np.concatenate([table[N][6] for N in Ns])
print(f"POOLED t_P/t_S = {allr.mean():.3f} ± {allr.std():.3f}  "
      f"<- number for the abstract")

plt.figure(figsize=(6, 4.2))
m = [table[N][2] for N in Ns]; e = [table[N][3] for N in Ns]
plt.errorbar(Ns, m, yerr=e, fmt="o-", capsize=4)
plt.xscale("log"); plt.axhline(1.0, color="k", lw=0.5, ls=":")
plt.xlabel("N"); plt.ylabel(r"$t_P/t_S$")
plt.title(f"Polarization-decoherence locking (R={R} realizations)")
plt.tight_layout(); plt.savefig(FIG + "fig5b_ratio_errorbars.png", dpi=150)

# ---------- 2. collapse residuals with bootstrap ----------
print("\n=== P(D) collapse residuals ===")
# master curve: pooled dense baseline over R realizations at N=64
Dp, Pp = [], []
for r in range(R):
    g = np.random.default_rng(2000+r).uniform(0.05, 0.6, size=64)  # wide g mix
    for t in np.linspace(0.05, 4.0, 200):
        d = exact.record_distinguishability(g, t)
        if d < 60:
            Dp.append(d)
            Pp.append(exact.sample_polarization(g, t, n_traj=20_000,
                      rng=np.random.default_rng(3000+r))[0])
Dp, Pp = np.array(Dp), np.array(Pp)
o = np.argsort(Dp)
# residuals: distance of each point to a running-median master curve
from numpy import interp
nb = max(8, min(40, len(Dp) // 8))          # quantile bins: equal occupancy
edges = np.quantile(Dp, np.linspace(0, 1, nb + 1))
idx = np.clip(np.digitize(Dp, edges[1:-1]), 0, nb - 1)
med_D = np.array([Dp[idx == i].mean() for i in range(nb)])
med_P = np.array([np.median(Pp[idx == i]) for i in range(nb)])
srt = np.argsort(med_D); med_D, med_P = med_D[srt], med_P[srt]
resid = Pp - interp(Dp, med_D, med_P)
boot = [np.abs(np.random.default_rng(b).choice(resid, size=len(resid))).max()
        for b in range(200)]
print(f"collapse: RMS residual = {np.sqrt((resid**2).mean()):.4f}, "
      f"max = {np.abs(resid).max():.4f} "
      f"(bootstrap 95% CI on max: [{np.percentile(boot, 2.5):.4f}, "
      f"{np.percentile(boot, 97.5):.4f}])")

# ---------- 3. hysteresis residual: sampling-floor test ----------
print("\n=== hysteresis residual vs n_traj (floor test) ===")
N = 16
g = np.random.default_rng(42).uniform(0.5*GBAR, 1.5*GBAR, size=N)
w = 2.0 * g
ts = np.linspace(0, 30, 121)
cs = [vb.overlaps(g, w, t) for t in ts]
for nt in [3_000, 30_000, 300_000]:
    P_B = np.array([vb.polarization_D(vb.eps_from_c(c), n_traj=nt,
                    rng=np.random.default_rng(77))[0] for c in cs])
    ep = [np.clip(vb.eps_from_c(c), 0, 1 - 1e-12) for c in cs]
    D_B = np.array([float(np.sum(e * np.log((1 + e) / (1 - e)))) for e in ep])
    ok = (D_B > med_D.min()) & (D_B < med_D.max())
    res = np.abs(P_B[ok] - interp(D_B[ok], med_D, med_P)).max()
    print(f"n_traj={nt:7d}: max residual = {res:.4f}")
print("interpretation: 1/sqrt(n) decay => sampling floor; plateau => physical")
print("\nDone. Paper numbers: pooled ratio above; collapse RMS; final residual.")
