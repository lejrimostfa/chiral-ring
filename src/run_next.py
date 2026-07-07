"""Run 2: collision protocol + sigma, finite-size scaling, data collapse."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import exact, collision

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures", "paper1", "")
rng = np.random.default_rng(42)

# ================= Part A: collision protocol, N = 120 =================
N, dt = 120, 0.6
g = rng.uniform(0.1, 0.3, size=N)
ks = np.arange(0, N + 1, 2)

S_S = np.array([collision.collision_S_entropy(g, dt, k) for k in ks])
D   = np.array([collision.collision_D(g, dt, k) for k in ks])
P, sig_mean, ift = np.zeros(len(ks)), np.zeros(len(ks)), np.zeros(len(ks))
for i, k in enumerate(ks):
    if k == 0:
        P[i], sig_mean[i], ift[i] = 0, 0, 1
        continue
    sg, pol, f = collision.collision_sample_sigma_pol(
        g, dt, k, n_traj=6000, rng=np.random.default_rng(100 + k))
    P[i], sig_mean[i], ift[i] = np.mean(pol), np.mean(sg), f
Rd = np.array([collision.collision_redundancy(
    g, dt, k, rng=np.random.default_rng(200 + k)) for k in ks])

id_err = np.max(np.abs(sig_mean - D) / np.maximum(D, 1e-9))
kW = (D > 0.5) & (D < 25)          # window where IFT estimator is reliable
ift_err = np.max(np.abs(ift[kW & (ks[:] <= 40)] - 1)) if np.any(kW) else 0.0
print(f"identity <sigma> = D : max relative deviation = {id_err:.3e}")
print(f"IFT <e^-sigma> = 1   : max |dev| over k<=40   = {ift_err:.3f} "
      f"(sampling-limited)")

t_ax = ks * dt
fig, ax1 = plt.subplots(figsize=(6.6, 4.3))
ax1.plot(t_ax, S_S, color="tab:gray", label=r"$S(\rho_S)$")
ax1.plot(t_ax, P, color="tab:red", label=r"$P$ polarization")
ax1.plot(t_ax, D / D.max(), color="tab:green", ls="--", label=r"$D$ (norm.)")
ax1.plot(t_ax, sig_mean / D.max(), color="k", ls=":", lw=2,
         label=r"$\langle\sigma\rangle$ sampled (norm.)")
ax1.set_xlabel("t = k dt  [1/J]"); ax1.set_ylabel("S, P, D (norm.)")
ax2 = ax1.twinx()
ax2.plot(t_ax, Rd, color="tab:blue", label=r"$R_{0.1}$")
ax2.set_ylabel(r"$R_\delta$", color="tab:blue")
h1, l1 = ax1.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2, fontsize=8, loc="center right")
plt.title(f"Collision protocol, N = {N}, dt = {dt} (monotone, no revivals)")
plt.tight_layout(); plt.savefig(FIG + "fig4_collision.png", dpi=150)

# ================= Part B: finite-size scaling (snapshot) =================
Ns = [8, 16, 32, 64, 128]
gbar = 0.2
res = {}
for Nn in Ns:
    gg = np.random.default_rng(5).uniform(0.5 * gbar, 1.5 * gbar, size=Nn)
    ts = np.linspace(0, 6.0, 241)
    Sv = np.array([exact.ent_rank2(np.prod(np.cos(2 * gg * t))) for t in ts])
    Pv = np.array([exact.sample_polarization(
        gg, t, n_traj=3000, rng=np.random.default_rng(21))[0] for t in ts])

    def cross(y, level):
        i = np.argmax(y >= level)
        if i == 0:
            return np.nan
        # linear interpolation
        return ts[i-1] + (level - y[i-1]) * (ts[i] - ts[i-1]) / (y[i] - y[i-1])

    # use the pre-revival rise only
    t_S50 = cross(Sv, 0.5 * 1.0)
    t_P50 = cross(Pv, 0.5)
    width = cross(Pv, 0.9) - cross(Pv, 0.1)
    res[Nn] = (t_S50, t_P50, width)
    print(f"N={Nn:4d}: t_S(50%)={t_S50:.3f}  t_P(50%)={t_P50:.3f}  "
          f"ratio={t_P50/t_S50:.3f}  width_P={width:.3f}")

Ns_a = np.array(Ns, float)
tS_a = np.array([res[n][0] for n in Ns])
tP_a = np.array([res[n][1] for n in Ns])
w_a = np.array([res[n][2] for n in Ns])
plt.figure(figsize=(6, 4.2))
plt.loglog(Ns_a, tS_a, "o-", label=r"$t_S$ (decoherence 50%)")
plt.loglog(Ns_a, tP_a, "s-", label=r"$t_P$ (polarization 50%)")
plt.loglog(Ns_a, w_a, "^-", label=r"width of $P$ rise (10–90%)")
guide = tS_a[0] * np.sqrt(Ns_a[0] / Ns_a)
plt.loglog(Ns_a, guide, "k:", label=r"$\propto N^{-1/2}$ guide")
plt.xlabel("N"); plt.ylabel("timescale  [1/J]")
plt.title("Finite-size scaling (snapshot protocol, g = 0.2 band)")
plt.legend(fontsize=8); plt.tight_layout()
plt.savefig(FIG + "fig5_scaling.png", dpi=150)

# ================= Part C: data collapse P vs D =================
plt.figure(figsize=(6, 4.2))
markers = ["o", "s", "^", "v", "D", "P"]
i = 0
for Nn in [16, 64]:
    for gb in [0.1, 0.2, 0.4]:
        gg = np.random.default_rng(9).uniform(0.5 * gb, 1.5 * gb, size=Nn)
        ts = np.linspace(0, 3.0 / gb, 40)
        Dv = np.array([exact.record_distinguishability(gg, t) for t in ts])
        Pv = np.array([exact.sample_polarization(
            gg, t, n_traj=3000, rng=np.random.default_rng(33))[0] for t in ts])
        # keep the first monotone rise of D (pre-revival)
        imax = np.argmax(Dv)
        plt.plot(Dv[:imax + 1], Pv[:imax + 1], markers[i % 6], ms=4, alpha=0.65,
                 label=f"N={Nn}, g={gb}")
        i += 1
# collision points too (different protocol entirely)
plt.plot(D[D <= 30], P[D <= 30], "k*", ms=6, alpha=0.7,
         label="collision, N=120")
plt.xscale("log"); plt.xlim(1e-2, 40)
plt.xlabel("accumulated record distinguishability D  [nats]")
plt.ylabel(r"branch polarization $P$")
plt.title("Data collapse test: P as a universal function of D")
plt.legend(fontsize=7); plt.tight_layout()
plt.savefig(FIG + "fig6_collapse.png", dpi=150)
print("figures 4-6 written")
