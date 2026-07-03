"""Rerun of errorbars.py Sections 2 (collapse) + 3 (hysteresis floor test)
with the DENSIFIED master curve (linspace 0.05..4.0, 200 points).

Statistical definitions are copied verbatim from errorbars.py; only the
master-curve sampling density (50 -> 200) changes. Section 1 (t_P/t_S,
the expensive part) is intentionally skipped -- it is unaffected.
"""
import numpy as np
from numpy import interp
import exact, variant_b as vb

R = 10
GBAR = 0.2

# ---------- 2. collapse residuals with bootstrap (master curve @ 200) ----
print("\n=== P(D) collapse residuals (master curve = 200 pts) ===")
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
nb = max(8, min(40, len(Dp) // 8))          # quantile bins: equal occupancy
edges = np.quantile(Dp, np.linspace(0, 1, nb + 1))
idx = np.clip(np.digitize(Dp, edges[1:-1]), 0, nb - 1)
med_D = np.array([Dp[idx == i].mean() for i in range(nb)])
med_P = np.array([np.median(Pp[idx == i]) for i in range(nb)])
srt = np.argsort(med_D); med_D, med_P = med_D[srt], med_P[srt]
resid = Pp - interp(Dp, med_D, med_P)
boot = [np.abs(np.random.default_rng(b).choice(resid, size=len(resid))).max()
        for b in range(200)]
print(f"n_master_pts = {len(Dp)}, n_bins = {nb}")
print(f"collapse: RMS residual = {np.sqrt((resid**2).mean()):.4f}, "
      f"max = {np.abs(resid).max():.4f} "
      f"(bootstrap 95% CI on max: [{np.percentile(boot, 2.5):.4f}, "
      f"{np.percentile(boot, 97.5):.4f}])")

# ---------- 3. hysteresis residual: sampling-floor test -----------------
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
