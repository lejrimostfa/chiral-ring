"""Baseline run: validation (N=8) + production curves R_delta(t), P(t) (N=16)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import run_unit_tests
import exact
import evolve
import os

rng = np.random.default_rng(42)
FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures", "paper1", "")

# ---------- 0. unit tests ----------
run_unit_tests()

# ---------- 1. validation at N = 8 ----------
N_v = 8
g_v = rng.uniform(0.1, 0.3, size=N_v)          # g = 0.2, spread [g/2, 3g/2]
ts_v = np.linspace(0.0, 12.0, 13)
H = evolve.build_hamiltonian(g_v)
psi0 = evolve.initial_state(N_v)
states = evolve.evolve(H, psi0, ts_v)

max_err = 0.0
frag_rng = np.random.default_rng(7)
for psi, t in zip(states, ts_v):
    for m in (1, 2, 4, 6):
        idx = frag_rng.choice(N_v, size=m, replace=False)
        mask = np.zeros(N_v, bool); mask[idx] = True
        Ia = exact.mutual_information(g_v, t, mask)
        In = evolve.mutual_information_num(psi, N_v, list(idx))
        max_err = max(max_err, abs(Ia - In))
print(f"validation I(S:F): max |analytic - numeric| = {max_err:.3e}")
assert max_err < 1e-9, "analytic/numeric mismatch"

# polarization: analytic sampling vs numeric sampling (statistical agreement)
t_chk = 6.0
psi_chk = evolve.evolve(H, psi0, [t_chk])[0]
Pa, _ = exact.sample_polarization(g_v, t_chk, n_traj=20000, rng=np.random.default_rng(3))
Pn, _ = evolve.sample_polarization_num(psi_chk, N_v, 20000, np.random.default_rng(4))
print(f"validation P(t={t_chk}): analytic {Pa:.4f} vs numeric {Pn:.4f} "
      f"(diff {abs(Pa-Pn):.4f}, sampling noise ~0.005)")
assert abs(Pa - Pn) < 0.02, "polarization sampling mismatch"

# ---------- 2. production run at N = 16 (analytic engine) ----------
N = 16
g = rng.uniform(0.1, 0.3, size=N)
ts = np.linspace(0.0, 25.0, 126)

S_S, Rd, P = [], [], []
for t in ts:
    cos_all = np.cos(2 * g * t)
    S_S.append(exact.ent_rank2(np.prod(cos_all)))
    Rd.append(exact.redundancy(g, t, delta=0.1, n_samples=150,
                               rng=np.random.default_rng(11)))
    p, _ = exact.sample_polarization(g, t, n_traj=4000,
                                     rng=np.random.default_rng(12))
    P.append(p)
S_S, Rd, P = map(np.asarray, (S_S, Rd, P))
D = np.array([exact.record_distinguishability(g, t) for t in ts])

# timescales (half-rise of asymptote over the pre-recurrence window)
w = ts <= 12.0
def half_rise(y, mask):
    ym = y[mask]; tgt = 0.5 * ym.max()
    return ts[mask][np.argmax(ym >= tgt)]
tP, tS = half_rise(P, w), half_rise(S_S, w)
tR = ts[w][np.argmax(Rd[w] >= 0.5 * Rd[w].max())]
print(f"timescales: t_S(decoherence) = {tS:.2f}, t_P(polarization) = {tP:.2f}, "
      f"t_R(redundancy half) = {tR:.2f}   [units 1/J]")

# ---------- 3. figures ----------
# fig 1: partial information plots
plt.figure(figsize=(6, 4.2))
for t, c in zip([1.0, 2.0, 4.0, 8.0], plt.cm.viridis(np.linspace(0.15, 0.85, 4))):
    pip = exact.partial_info_plot(g, t, n_samples=150, rng=np.random.default_rng(5))
    ssv = exact.ent_rank2(np.prod(np.cos(2 * g * t)))
    plt.plot(np.arange(N + 1), pip, "-o", ms=3, color=c, label=f"t = {t:g}")
    plt.axhline(ssv, color=c, lw=0.6, ls=":")
plt.xlabel("fragment size m"); plt.ylabel(r"$\bar I(S{:}F_m)$  [bits]")
plt.title("Partial information plots (N = 16)")
plt.legend(); plt.tight_layout(); plt.savefig(FIG + "fig1_pip.png", dpi=150)

# fig 2: co-transition curves
fig, ax1 = plt.subplots(figsize=(6.4, 4.2))
ax1.plot(ts, S_S, label=r"$S(\rho_S)$ (decoherence)", color="tab:gray")
ax1.plot(ts, P, label=r"$P(t)$ branch polarization", color="tab:red")
ax1.plot(ts, D / D[w].max(), label=r"$D(t)$ record disting. (norm.)",
         color="tab:green", ls="--")
ax1.set_xlabel("t  [1/J]"); ax1.set_ylabel("S, P, D (normalized)")
ax2 = ax1.twinx()
ax2.plot(ts, Rd, label=r"$R_{0.1}(t)$ redundancy", color="tab:blue")
ax2.set_ylabel(r"$R_{\delta}$", color="tab:blue")
h1, l1 = ax1.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2, loc="upper center",
           bbox_to_anchor=(0.5, -0.18), ncol=4, fontsize=8, frameon=False)
plt.title("Co-transition, baseline N = 16 (random couplings)")
plt.tight_layout(); plt.savefig(FIG + "fig2_cotransition_v2.png", dpi=150)

# fig 3: chirality histograms (unimodal -> bimodal)
fig, axes = plt.subplots(1, 4, figsize=(11, 2.8), sharey=True)
for ax, t in zip(axes, [0.1, 0.5, 1.5, 8.0]):
    _, pol = exact.sample_polarization(g, t, n_traj=20000,
                                       rng=np.random.default_rng(13))
    ax.hist(pol, bins=60, range=(-1, 1), density=True, color="tab:red", alpha=0.8)
    ax.set_title(f"t = {t:g}"); ax.set_xlabel(r"$\langle I\rangle_\xi / I_0$")
axes[0].set_ylabel("density")
plt.suptitle("Branch-level chirality: spontaneous T-symmetry breaking", y=1.02)
plt.tight_layout(); plt.savefig(FIG + "fig3_histograms_v2.png", dpi=150,
                                bbox_inches="tight")
print("figures written to", FIG)
