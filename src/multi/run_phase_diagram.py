"""Phase diagram f*(t) of the arrow-frustration transition (final task A).

For each shared-qubit prescription F / M / H (multi_definitions.md §6bis):
sweep t in [0.05, 1.5] (30 points); at each t extract the sharing fraction f*
where an arrow flips (interpolated zero crossing of sigma_2(f) AND of
sigma_1(f); we record WHICH ring flips). M=2, N=64, R=10 disorder
realizations, n_traj=20k, per-realization f* -> mean ± std.

T_CHARGE = 1.0 (fixed, documented): differs deliberately from the T=1.5 of
run_c5_robustness.py. With T=1.0 the sweep crosses
  - the EQUAL-EFFORT point t = T/2 = 0.5, where forward and mirror record
    qualities are exactly comparable (|sin(2g t)| vs |sin(2g (T-t))| coincide),
    which answers whether the F-vs-M asymmetry (0.76 vs 0.10 measured at
    T=1.5, t=0.5 — an UNEQUAL-effort comparison: mirror reservoir stronger)
    persists at equal effort or was a head-start artefact;
  - the H-shared exhaustion point t = T/2 (H-shared qubits carry zero bias
    exactly at the equal-effort time: pure dilution);
  - the full-consumption point t = T (beyond it every charge has unwound and
    re-accumulates while the §6bis scoring keeps the backward reference).

Output: fig_multi_phase_diagram.png + printed f*(t->0), f*(t=T/2), f*(t->T).
Run:  python3 run_phase_diagram.py    (from src/multi)
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
N = 64
T_CHARGE = 1.0
R = 10
NTRAJ = 20_000
TS = np.linspace(0.05, 1.5, 30)
FS = np.linspace(0.0, 1.0, 11)
PRE = ["F", "M", "H"]


def fstars_one(p, t, seed):
    """(f*_sigma1, f*_sigma2) for one prescription/time/realization (nan if none)."""
    s1 = np.empty(len(FS)); s2 = np.empty(len(FS))
    for j, f in enumerate(FS):
        g, _ = rs.couplings_shared(N, f, np.random.default_rng(seed))
        d = px.sample_frustration_rx(g, t, px.charge_profile(p, g, T_CHARGE),
                                     NTRAJ, np.random.default_rng(99))
        s1[j], s2[j] = d["sigma"]
    c1 = zero_crossings(FS, s1)
    c2 = zero_crossings(FS, s2)
    return (c1[0] if c1 else np.nan), (c2[0] if c2 else np.nan)


def main():
    print(f"=== Phase diagram f*(t): F/M/H (N={N}, T_charge={T_CHARGE}, "
          f"R={R}, {len(TS)} t-points) ===")
    # res[p]: (R, len(TS), 2) -> f*_1, f*_2
    res = {p: np.full((R, len(TS), 2), np.nan) for p in PRE}
    for p in PRE:
        for r in range(R):
            for k, t in enumerate(TS):
                res[p][r, k] = fstars_one(p, t, 6000 + r)
        print(f"  [{p}] done")

    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    colors = {"F": "C0", "M": "C3", "H": "C2"}
    summary = {}
    for p in PRE:
        for ring, ls, mk in ((1, "-", "o"), (0, "--", "s")):   # idx1=sigma_2
            arr = res[p][:, :, ring]                            # (R, T)
            nok = np.sum(~np.isnan(arr), axis=0)
            m = np.nanmean(arr, axis=0)
            s = np.nanstd(arr, axis=0)
            valid = nok >= R // 2          # plot where >= half realizations cross
            if valid.any():
                lbl = f"{p}: $\\sigma_{2 if ring == 1 else 1}$ flips"
                ax.plot(TS[valid], m[valid], ls, marker=mk, ms=3.5,
                        color=colors[p], label=lbl)
                ax.fill_between(TS[valid], (m - s)[valid], (m + s)[valid],
                                color=colors[p], alpha=0.15)
                summary[(p, ring)] = (m, s, nok)
    ax.axvline(T_CHARGE / 2, color="k", lw=0.7, ls=":",
               label=f"equal effort t=T/2={T_CHARGE/2}")
    ax.axvline(T_CHARGE, color="gray", lw=0.7, ls=":",
               label=f"charge exhausted t=T={T_CHARGE}")
    ax.set_xlabel("t"); ax.set_ylabel("shared fraction f")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title(f"Arrow-frustration phase boundary f*(t)  "
                 f"(N={N}, T_charge={T_CHARGE}, ±1σ over {R} realizations)")
    # region annotations (F: sigma_2 flips above the F curve)
    ax.annotate("mirror arrow intact", xy=(0.18, 0.90), fontsize=8, color="C0")
    ax.annotate("mirror arrow reversed\n(above F curve)", xy=(0.95, 0.9),
                fontsize=8, color="C0")
    ax.annotate("forward arrow reversed\n(above M/H curves)", xy=(0.35, 0.30),
                fontsize=8, color="C3")
    ax.legend(fontsize=7, loc="center right")
    fig.tight_layout()
    fig.savefig(FIG + "fig_multi_phase_diagram.png", dpi=150)
    plt.close(fig)

    # boundary values and the equal-effort comparison
    k0, kh, kT = 0, np.argmin(np.abs(TS - T_CHARGE / 2)), \
        np.argmin(np.abs(TS - T_CHARGE))
    print("\n  boundary values (nan = no crossing in >= half realizations):")
    for (p, ring), (m, s, nok) in summary.items():
        name = "sigma_2" if ring == 1 else "sigma_1"
        def cell(k):
            return (f"{m[k]:.3f}±{s[k]:.3f}" if nok[k] >= R // 2 else "none")
        print(f"    {p} ({name} flips): f*(t={TS[k0]:.2f})={cell(k0)}  "
              f"f*(t=T/2={TS[kh]:.2f})={cell(kh)}  "
              f"f*(t~T={TS[kT]:.2f})={cell(kT)}")
    print("\n  EQUAL-EFFORT question (t = T/2): compare F (mirror flips) vs "
          "M (forward flips) above; if both f* remain far apart at t=T/2, the "
          "F-vs-M asymmetry is structural, not a head-start artefact.")


if __name__ == "__main__":
    main()
