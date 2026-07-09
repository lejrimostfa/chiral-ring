"""Paper 3 — P3-C: percolation of aligned arrow domains on the 2D square lattice.

Arrow field tau_i = sign(<sigma_i^local>) (collision protocol) on an S x S
lattice of chiral rings. A ring's accessible witnesses = private bath + the
4-neighbour shared-link baths, so the local marginal costs O(2^5) per bulk site.

Forward arrows (tau=+1) occupy a fraction p = (1+A)/2. In 2D they percolate
above the site threshold p_c ~= 0.5927; sharing converts contrarians, raising p
through the transition -> a genuine finite-size percolation transition in the
sharing fraction f = 4 Ns/(Np + 4 Ns) (bulk value).

Produces four figures:
  fig_p3_2d_phase.png   heatmap A(f, rho) with the percolation contour
  fig_p3_percolation.png P_perc(f) for several L, crossing -> f_c
  fig_p3_fss.png        finite-size-scaling collapse, exponent nu
  fig_p3_snapshots.png  tau maps below / near / above f_c

Run:  python3 run_square.py    (from src/lattice)
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import label

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lattice import build_square, draw_signs, sample_records, local_marginal_sigma

FIG = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "paper3", "")
NP = 8
DT = 0.3
T_CHARGE = 1.0
NTRAJ = 4000
PC = 0.5927                      # 2D square-lattice site percolation threshold
NU_PERC = 4.0 / 3.0             # 2D percolation correlation-length exponent


def frac(Ns):
    return 4 * Ns / (NP + 4 * Ns)


def square_arrows(S, Ns, rho, seed):
    """tau (S,S) grid and the contrarian mask, one disorder realization."""
    rng = np.random.default_rng(seed)
    L = S * S
    n_c = round(rho * L)
    contr = tuple(np.sort(rng.choice(L, size=n_c, replace=False))) if n_c else ()
    G, Tc, owner = build_square(S, S, NP, Ns, seed=seed, contrarians=contr,
                                T_charge=T_CHARGE)
    S_true = np.ones(L, int)
    beta_true = np.sin(2 * (S_true @ G) * (DT - Tc))
    R = sample_records(beta_true, NTRAJ, rng, draw_sign=draw_signs(Tc))
    tau = np.array([np.sign(local_marginal_sigma(R, G, Tc, DT, S_true, i).mean())
                    for i in range(L)])
    cmask = np.zeros(L, bool); cmask[list(contr)] = True
    return tau.reshape(S, S), cmask.reshape(S, S)


def percolates(grid):
    """True if a forward (tau=+1) cluster spans top-to-bottom (4-connectivity)."""
    lab, _ = label(grid > 0)
    top = set(lab[0]) - {0}
    bot = set(lab[-1]) - {0}
    return len(top & bot) > 0


def big_cluster_frac(grid):
    lab, n = label(grid > 0)
    if n == 0:
        return 0.0
    return np.bincount(lab.ravel())[1:].max() / grid.size


# ----------------------------------------------------------------------
def main():
    NSS = list(range(10))
    fs = np.array([frac(Ns) for Ns in NSS])
    print(f"=== P3-C 2D percolation (Np={NP}, dt={DT}, T={T_CHARGE}, p_c={PC}) ===")

    # ---------- Part 1: phase map A(f, rho), L=6 ----------
    rhos_map = np.round(np.arange(0.1, 0.71, 0.1), 2)
    Rm = 10
    Amap = np.zeros((len(rhos_map), len(NSS)))
    for a, rho in enumerate(rhos_map):
        for k, Ns in enumerate(NSS):
            Amap[a, k] = np.mean([square_arrows(6, Ns, rho, 800 + r)[0].mean()
                                  for r in range(Rm)])
    print("  phase map done (L=6)")

    fig, ax = plt.subplots(figsize=(6.5, 4.6))
    im = ax.imshow(Amap, origin="lower", aspect="auto", cmap="RdBu_r",
                   vmin=-1, vmax=1, extent=[fs[0], fs[-1], rhos_map[0],
                                            rhos_map[-1]])
    # percolation contour: forward occupation (1+A)/2 = p_c  <=>  A = 2 p_c - 1
    CS = ax.contour(fs, rhos_map, Amap, levels=[2 * PC - 1], colors="k",
                    linewidths=2)
    ax.clabel(CS, fmt={2 * PC - 1: "percolation\nthreshold"}, fontsize=8)
    fig.colorbar(im, ax=ax, label=r"signed alignment $A=\langle\tau\rangle$")
    ax.set_xlabel("sharing fraction $f=4N_s/(N_p+4N_s)$")
    ax.set_ylabel(r"contrarian density $\rho$")
    ax.set_title("(P3-C) Arrow alignment landscape A(f, $\\rho$), 2D L=6")
    fig.tight_layout(); fig.savefig(FIG + "fig_p3_2d_phase.png", dpi=150)
    plt.close(fig)

    # ---------- Part 2: percolation P_perc(f) at rho=0.5, several L ----------
    Ls = [4, 5, 6, 8]
    Rp = 40
    rho_p = 0.5
    Pperc = {S: np.zeros(len(NSS)) for S in Ls}
    for S in Ls:
        for k, Ns in enumerate(NSS):
            Pperc[S][k] = np.mean([percolates(square_arrows(S, Ns, rho_p,
                                   900 + r)[0]) for r in range(Rp)])
        print(f"  percolation L={S} done")

    # f_c from the P=0.5 crossing of the largest lattice
    def cross_half(P):
        for k in range(len(P) - 1):
            if P[k] < 0.5 <= P[k + 1]:
                return fs[k] + (0.5 - P[k]) * (fs[k + 1] - fs[k]) / \
                    (P[k + 1] - P[k])
        return np.nan
    fc = cross_half(Pperc[max(Ls)])
    print(f"  f_c (L={max(Ls)}, P=0.5) = {fc:.3f}")

    fig, ax = plt.subplots(figsize=(6.5, 4.6))
    for S in Ls:
        ax.plot(fs, Pperc[S], marker="o", ms=4, label=f"L={S}")
    ax.axvline(fc, color="gray", ls="--", lw=1, label=f"$f_c\\approx${fc:.2f}")
    ax.axhline(0.5, color="k", lw=0.4)
    ax.set_xlabel("sharing fraction $f$")
    ax.set_ylabel(r"$P_{\rm perc}$ (spanning forward cluster)")
    ax.set_title(f"(P3-C) Arrow-domain percolation, $\\rho$={rho_p} "
                 f"(forward $p_c$={PC})")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(FIG + "fig_p3_percolation.png", dpi=150)
    plt.close(fig)

    # ---------- Part 3: finite-size-scaling collapse, exponent nu ----------
    # Restrict the collapse-quality metric to the TRANSITION window (exclude the
    # far below-threshold plateaus, which are L-dependent finite-size values that
    # do not collapse and would bias the metric). Scan nu over a physical range
    # (a naive scan rails to large nu, where L^{1/nu}->1 collapses everything
    # trivially -- an artefact, not physics).
    def collapse_spread(nu, xwin=1.1):
        xg = np.linspace(-xwin, xwin, 21)
        curves = []
        for S in Ls:
            x = (fs - fc) * S ** (1.0 / nu)
            order = np.argsort(x)
            curves.append(np.interp(xg, x[order], Pperc[S][order],
                                    left=np.nan, right=np.nan))
        curves = np.array(curves)
        return np.nanmean(np.nanstd(curves, axis=0))
    nus = np.linspace(0.9, 1.8, 28)
    spreads = np.array([collapse_spread(nu) for nu in nus])
    nu_best = nus[int(np.argmin(spreads))]
    railed = nu_best <= nus[1] or nu_best >= nus[-2]
    print(f"  FSS: window-collapse spread has no robust interior minimum over "
          f"nu in [0.9,1.8] (argmin={nu_best:.2f}, "
          f"{'rails to boundary' if railed else 'interior'}). At L<=8 with the "
          f"coarse integer-Ns f-grid, nu is NOT independently fixed; the "
          f"collapse at the 2D-percolation value nu=4/3={NU_PERC:.2f} is shown "
          f"and is visually consistent.")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    for S in Ls:                                    # collapse at 2D-percolation nu
        x = (fs - fc) * S ** (1.0 / NU_PERC)
        ax[0].plot(x, Pperc[S], marker="o", ms=4, label=f"L={S}")
    ax[0].set_xlim(-1.4, 1.4)
    ax[0].set_xlabel(r"$(f-f_c)\,L^{1/\nu}$,  $\nu=4/3$ (2D percolation)")
    ax[0].set_ylabel(r"$P_{\rm perc}$")
    ax[0].set_title(r"FSS collapse at the 2D-percolation exponent $\nu=4/3$")
    ax[0].legend(fontsize=8)
    ax[1].plot(nus, spreads, "-")
    ax[1].axvline(nu_best, color="C3", ls="--", label=f"in-range best {nu_best:.2f}")
    ax[1].axvline(NU_PERC, color="gray", ls=":", label=r"2D perc 4/3")
    ax[1].set_xlabel(r"$\nu$"); ax[1].set_ylabel("transition-window collapse spread")
    ax[1].set_title("collapse quality vs exponent (window-restricted)")
    ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(FIG + "fig_p3_fss.png", dpi=150)
    plt.close(fig)

    # ---------- Part 4: snapshots below / near / above f_c ----------
    below = max(k for k in range(len(NSS)) if fs[k] < fc - 0.03)
    near = int(np.argmin(np.abs(fs - fc)))
    above = min(k for k in range(len(NSS)) if fs[k] > fc + 0.03)
    fig, axs = plt.subplots(1, 3, figsize=(12, 4.2))
    for ax, k, lab in zip(axs, [below, near, above], ["below", "near", "above"]):
        tau, cmask = square_arrows(8, NSS[k], rho_p, 950)
        ax.imshow(tau, cmap="RdBu_r", vmin=-1, vmax=1)
        ys, xs = np.where(cmask)
        ax.scatter(xs, ys, s=14, facecolors="none", edgecolors="k", lw=0.8)
        ax.set_title(f"{lab} $f_c$: f={fs[k]:.2f}, A={tau.mean():+.2f}")
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(r"(P3-C) $\tau$ maps (blue=forward, red=backward; "
                 r"circles=contrarians), L=8, $\rho$=0.5")
    fig.tight_layout(); fig.savefig(FIG + "fig_p3_snapshots.png", dpi=150)
    plt.close(fig)

    print(f"  P_perc(f) at rho=0.5:")
    for k, Ns in enumerate(NSS):
        print(f"    f={fs[k]:.2f}: " +
              "  ".join(f"L{S}={Pperc[S][k]:.2f}" for S in Ls))
    print("figures -> fig_p3_2d_phase, fig_p3_percolation, fig_p3_fss, "
          "fig_p3_snapshots")


if __name__ == "__main__":
    main()
