"""C5 robustness to the shared-qubit prescription (consolidation T1).

M=2, N=64, f in {0, 0.1, ..., 1.0}, t=0.5 fixed (same as the existing crossing
figure), T=1.5, R=10 disorder realizations (same seed family as scenario e),
n_traj=20k. Prescriptions F (fresh), M (mirror-charged), H (half-charged) on
the shared qubits; exclusives fixed.

Outputs: fig_multi_c5_prescriptions.png (sigma_2(f) panel a, sigma_1(f) panel b,
+-1 std disorder bands) and printed f* (interpolated zero crossing) per
prescription and per ring where it exists.

Run:  python3 run_c5_robustness.py    (from src/multi)
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import run_scenarios as rs
import prescriptions as px

FIG = rs.FIG
N = 64
T = 1.5
T_EVAL = 0.5
R = 10
NTRAJ = 20_000
FS = np.linspace(0.0, 1.0, 11)
PRESCRIPTIONS = ["F", "M", "H"]


def zero_crossings(f, y):
    """All sign-change crossings of y(f), linearly interpolated."""
    out = []
    for i in np.where(np.diff(np.sign(y)) != 0)[0]:
        out.append(f[i] + (0 - y[i]) * (f[i + 1] - f[i]) / (y[i + 1] - y[i]))
    return out


def main():
    print(f"=== C5 robustness: prescriptions F/M/H (N={N}, t={T_EVAL}, "
          f"T={T}, R={R}) ===")
    # curves[p] : (R, len(FS), 2) sigma_1, sigma_2
    curves = {p: np.zeros((R, len(FS), 2)) for p in PRESCRIPTIONS}
    for r in range(R):
        for j, f in enumerate(FS):
            g, _ = rs.couplings_shared(N, f, np.random.default_rng(6000 + r))
            for p in PRESCRIPTIONS:
                d = px.sample_frustration_rx(
                    g, T_EVAL, px.charge_profile(p, g, T),
                    NTRAJ, np.random.default_rng(99))
                curves[p][r, j] = d["sigma"]

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4), sharex=True)
    colors = {"F": "C0", "M": "C3", "H": "C2"}
    for p in PRESCRIPTIONS:
        m = curves[p].mean(axis=0)      # (len(FS), 2)
        s = curves[p].std(axis=0)
        # panel 0 = sigma_2 (mirror ring), panel 1 = sigma_1 (forward ring)
        for k, ring in ((0, 1), (1, 0)):
            ax[k].plot(FS, m[:, ring], "o-", color=colors[p], label=f"{p}")
            ax[k].fill_between(FS, m[:, ring] - s[:, ring],
                               m[:, ring] + s[:, ring],
                               color=colors[p], alpha=0.18)
        # f* per ring from the disorder-mean curve
        for ring, name in ((1, "sigma_2"), (0, "sigma_1")):
            cr = zero_crossings(FS, m[:, ring])
            txt = ", ".join(f"{c:.3f}" for c in cr) if cr else "none"
            print(f"  {p}: f*({name}) = {txt}")

    for k, ttl in ((0, r"$\sigma_2$ (ring 2, mirror side)"),
                   (1, r"$\sigma_1$ (ring 1, forward side)")):
        ax[k].axhline(0, color="k", lw=0.6)
        ax[k].set_xlabel("shared fraction f"); ax[k].set_ylabel("nats")
        ax[k].set_title(ttl); ax[k].legend(title="prescription")
    fig.suptitle(f"(C5 robustness) shared-qubit prescriptions F / M / H "
                 f"(t={T_EVAL}, T={T}, ±1σ disorder)")
    fig.tight_layout()
    fig.savefig(FIG + "fig_multi_c5_prescriptions.png", dpi=150)
    plt.close(fig)
    print("figure -> fig_multi_c5_prescriptions.png")


if __name__ == "__main__":
    main()
