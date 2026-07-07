"""Run 3 — Variant B: degrading records, hysteresis test, ratchet contrast."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import exact, variant_b as vb

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures", "")
rng = np.random.default_rng(42)

N, gbar = 16, 0.2
g = rng.uniform(0.5 * gbar, 1.5 * gbar, size=N)

# ---------- fig 7: rise and fall for several omega/g ----------
fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8), sharey=False)
ratios = [0.5, 1.0, 2.0]
ts = np.linspace(0, 30, 181)
curves = {}
for ax, rat in zip(axes, ratios):
    w = rat * g
    S, P, D, R = [], [], [], []
    for t in ts:
        c = vb.overlaps(g, w, t)
        S.append(vb.S_system(c))
        p, d = vb.polarization_D(vb.eps_from_c(c), n_traj=3000,
                                 rng=np.random.default_rng(17))
        P.append(p); D.append(d)
        R.append(vb.redundancy_B(c, rng=np.random.default_rng(18)))
    S, P, D, R = map(np.asarray, (S, P, D, R))
    curves[rat] = (S, P, D, R)
    ax.plot(ts, S, color="tab:gray", label=r"$S(\rho_S)$")
    ax.plot(ts, P, color="tab:red", label=r"$P$")
    ax.plot(ts, D / max(D.max(), 1e-9), color="tab:green", ls="--",
            label=r"$D$ (norm.)")
    ax.plot(ts, R / R.max(), color="tab:blue", label=r"$R_{0.1}$ (norm.)")
    ax.set_title(fr"$\omega/g$ = {rat}")
    ax.set_xlabel("t  [1/J]")
axes[0].set_ylabel("normalized observables")
axes[0].legend(fontsize=8)
plt.suptitle("Variant B: record degradation and revival", y=1.02)
plt.tight_layout()
plt.savefig(FIG + "fig7_variantB.png", dpi=150, bbox_inches="tight")

# ---------- fig 8: hysteresis test on the master curve ----------
# master curve from baseline (snapshot, no H_E): dense (D, P) samples
gg = np.random.default_rng(9).uniform(0.1, 0.3, size=64)
tsm = np.linspace(0, 6, 60)
Dm = np.array([exact.record_distinguishability(gg, t) for t in tsm])
Pm = np.array([exact.sample_polarization(gg, t, n_traj=3000,
               rng=np.random.default_rng(33))[0] for t in tsm])
im = np.argmax(Dm)
order = np.argsort(Dm[:im + 1])

plt.figure(figsize=(6.2, 4.3))
plt.plot(Dm[:im + 1][order], Pm[:im + 1][order], "k-", lw=2, alpha=0.6,
         label="master curve (baseline)")
cols = {0.5: "tab:orange", 1.0: "tab:red", 2.0: "tab:purple"}
for rat in ratios:
    S, P, D, R = curves[rat]
    plt.plot(D, P, ".", ms=4, color=cols[rat], alpha=0.6,
             label=fr"variant B, $\omega/g$={rat} (up & down)")
plt.xscale("log"); plt.xlim(3e-2, 50)
plt.xlabel("D  [nats]"); plt.ylabel("P")
plt.title("Hysteresis test: does P ride the master curve both ways?")
plt.legend(fontsize=8); plt.tight_layout()
plt.savefig(FIG + "fig8_hysteresis.png", dpi=150)

# quantify: residuals of variant-B points vs master curve (interpolated)
from numpy import interp
mask_mono = np.argsort(Dm[:im + 1])
Dmm, Pmm = Dm[:im + 1][mask_mono], Pm[:im + 1][mask_mono]
resid = []
for rat in ratios:
    S, P, D, R = curves[rat]
    ok = (D > Dmm.min()) & (D < Dmm.max())
    resid.append(np.max(np.abs(P[ok] - interp(D[ok], Dmm, Pmm))))
print("hysteresis residuals |P - P_master(D)| max:",
      {r: f"{v:.3f}" for r, v in zip(ratios, resid)})

# ---------- fig 9: ratchet contrast (snapshot vs collision records) ----------
# collision: qubit n interacts during [t_{n-1}, t_n] under variant-B conditional
# dynamics, measured at t_n -> eps frozen at eps_n(dt); snapshot: all evolve.
rat = 2.0
w = rat * g
dt = 1.0
Nc = 60
gC = np.tile(g, 4)[:Nc]; wC = np.tile(w, 4)[:Nc]
eps_frozen = np.array([vb.eps_from_c(np.array([vb.overlaps(
    np.array([gC[n]]), np.array([wC[n]]), dt)[0]]))[0] for n in range(Nc)])
tsr = np.arange(0, Nc + 1) * dt
P_ratchet, D_ratchet = [0.0], [0.0]
for k in range(1, Nc + 1):
    p, d = vb.polarization_D(eps_frozen[:k], n_traj=3000,
                             rng=np.random.default_rng(55))
    P_ratchet.append(p); D_ratchet.append(d)
P_snap, D_snap = [], []
for t in tsr:
    c = vb.overlaps(g, w, t)
    p, d = vb.polarization_D(vb.eps_from_c(c), n_traj=3000,
                             rng=np.random.default_rng(56))
    P_snap.append(p); D_snap.append(d)

plt.figure(figsize=(6.2, 4.3))
plt.plot(tsr, P_snap, color="tab:red",
         label="coherent records (snapshot): revocable")
plt.plot(tsr, P_ratchet, color="k",
         label="measured records (collision): ratchet")
plt.xlabel("t  [1/J]"); plt.ylabel("P")
plt.title(fr"The arrow ratchet ($\omega/g$ = {rat})")
plt.legend(fontsize=9); plt.tight_layout()
plt.savefig(FIG + "fig9_ratchet.png", dpi=150)
print("figures 7-9 written")
