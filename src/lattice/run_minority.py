"""Paper 3 — P3-D: weighted minority + consolidation.

Two deliverables:
  (1) rho_c(f): the critical contrarian density that reverses the GLOBAL arrow
      (A(f,rho)=0). Tests the "shared witnesses vote double" idea against the
      naive rho_c = 1/2.
  (2) A clean finite-size-scaling of the percolation transition, now driven by
      the CONTINUOUS control parameter rho (fine steps 1/L^2) at fixed sharing,
      with larger lattices L in {8,12,16}. This pins the exponent nu that the
      coarse f-grid of P3-C could not (retroactively strengthening P3-C).

Figures: figures/paper3/fig_p3_rhoc.png, fig_p3_fss_rho.png
Run:  python3 run_minority.py    (from src/lattice)
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import label

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lattice import (build_square, draw_signs, sample_records,
                     local_marginal_sigma, torus_percolation)

FIG = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "paper3", "")
NP = 8
DT = 0.3
T_CHARGE = 1.0
NTRAJ = 1500
PC = 0.5927
NU_PERC = 4.0 / 3.0


def frac(Ns):
    return 4 * Ns / (NP + 4 * Ns)


def square_tau(S, Ns, rho, seed):
    rng = np.random.default_rng(seed)
    L = S * S
    n_c = round(rho * L)
    contr = tuple(np.sort(rng.choice(L, size=n_c, replace=False))) if n_c else ()
    G, Tc, owner = build_square(S, S, NP, Ns, seed=seed, contrarians=contr,
                                T_charge=T_CHARGE)
    S_true = np.ones(L, int)
    beta = np.sin(2 * (S_true @ G) * (DT - Tc))
    R = sample_records(beta, NTRAJ, rng, draw_sign=draw_signs(Tc))
    tau = np.array([np.sign(local_marginal_sigma(R, G, Tc, DT, S_true, i).mean())
                    for i in range(L)])
    return tau.reshape(S, S)


def percolates(grid):
    """Forward cluster wraps the torus (periodic BC)."""
    return torus_percolation(grid > 0)[0]


def zero_cross(x, y, target=0.0):
    for k in range(len(y) - 1):
        if (y[k] - target) * (y[k + 1] - target) <= 0 and y[k] != y[k + 1]:
            return x[k] + (target - y[k]) * (x[k + 1] - x[k]) / (y[k + 1] - y[k])
    return np.nan


def main():
    print(f"=== P3-D minority + FSS (Np={NP}, dt={DT}, T={T_CHARGE}) ===")

    # ---------- (1) rho_c(f): global-arrow reversal ----------
    Nss = [0, 2, 3, 4, 5, 6, 7, 8, 10]
    rhos = np.round(np.arange(0.30, 0.93, 0.03), 3)
    Rr, Sr = 25, 8
    rho_c, diverged = [], []
    for Ns in Nss:
        Acurve = np.array([np.mean([square_tau(Sr, Ns, rho, 300 + r).mean()
                                    for r in range(Rr)]) for rho in rhos])
        rc = zero_cross(rhos, Acurve, 0.0)
        # nan + still-forward at the densest tested rho => rho_c beyond range:
        # the forward arrow wins for EVERY tested contrarian density (divergence)
        div = np.isnan(rc) and Acurve[-1] > 0
        rho_c.append(rhos[-1] if div else rc)
        diverged.append(div)
        tag = f">{rhos[-1]:.2f} (forward wins for all tested rho)" if div \
            else f"{rc:.3f}"
        print(f"  f={frac(Ns):.2f} (Ns={Ns}): rho_c(A=0) = {tag}")
    rho_c = np.array(rho_c); diverged = np.array(diverged)
    fsN = np.array([frac(Ns) for Ns in Nss])

    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    fin = ~diverged
    ax.plot(fsN[fin], rho_c[fin], "o-", color="C0",
            label=r"$\rho_c(f)$ (global arrow reversal)")
    if diverged.any():
        ax.plot(fsN[diverged], rho_c[diverged], "^", color="C3", ms=10,
                label=r"$\rho_c$ diverges (forward wins $\forall\rho$)")
        for x in fsN[diverged]:
            ax.annotate("", xy=(x, rhos[-1] + 0.03), xytext=(x, rhos[-1] - 0.02),
                        arrowprops=dict(arrowstyle="->", color="C3"))
    ax.axhline(0.5, color="k", ls="--", lw=1, label=r"naive $\rho_c=1/2$")
    ax.fill_between(fsN[fin], 0.5, rho_c[fin], where=(rho_c[fin] > 0.5),
                    alpha=0.15, color="C0",
                    label="forward survives a contrarian majority")
    ax.set_xlabel("sharing fraction $f$")
    ax.set_ylabel(r"critical contrarian density $\rho_c$")
    ax.set_ylim(0.47, rhos[-1] + 0.05)
    ax.set_title("(P3-D) Shared witnesses vote double: "
                 r"$\rho_c(f)$ climbs off 1/2 and diverges")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout(); fig.savefig(FIG + "fig_p3_rhoc.png", dpi=150)
    plt.close(fig)

    # ---------- (2) clean FSS: rho-driven percolation at fixed f ----------
    Ns_fss = 2                                   # f = 0.5
    Ls = [8, 12, 16]
    rho_g = np.round(np.linspace(0.30, 0.55, 12), 4)
    Rp = 30
    Pperc = {S: np.zeros(len(rho_g)) for S in Ls}
    for S in Ls:
        for k, rho in enumerate(rho_g):
            Pperc[S][k] = np.mean([percolates(square_tau(S, Ns_fss, rho, 400 + r))
                                   for r in range(Rp)])
        print(f"  FSS L={S} done")
    # rho_c^perc from the P=0.5 crossing of the largest lattice
    rcp = zero_cross(rho_g, Pperc[max(Ls)], 0.5)
    print(f"  rho_c^perc (L={max(Ls)}, f={frac(Ns_fss):.2f}) = {rcp:.3f}")

    def spread(nu, win=1.2):
        xg = np.linspace(-win, win, 21)
        cur = []
        for S in Ls:
            x = (rho_g - rcp) * S ** (1.0 / nu)
            o = np.argsort(x)
            cur.append(np.interp(xg, x[o], Pperc[S][o], left=np.nan, right=np.nan))
        return np.nanmean(np.nanstd(np.array(cur), axis=0))
    nus = np.linspace(0.9, 2.0, 34)
    sp = np.array([spread(nu) for nu in nus])
    nu_best = nus[int(np.argmin(sp))]
    railed = nu_best <= nus[1] or nu_best >= nus[-2]
    print(f"  FSS best-collapse nu = {nu_best:.2f} "
          f"({'RAILS' if railed else 'interior min'}); 2D perc 4/3={NU_PERC:.2f}")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    for S in Ls:
        x = (rho_g - rcp) * S ** (1.0 / nu_best)
        ax[0].plot(x, Pperc[S], "o-", ms=4, label=f"L={S}")
    ax[0].set_xlabel(r"$(\rho-\rho_c)\,L^{1/\nu}$")
    ax[0].set_ylabel(r"$P_{\rm perc}$")
    ax[0].set_title(f"rho-driven FSS collapse, $\\nu$={nu_best:.2f}")
    ax[0].legend(fontsize=8)
    ax[1].plot(nus, sp, "-")
    ax[1].axvline(nu_best, color="C3", ls="--", label=f"best {nu_best:.2f}")
    ax[1].axvline(NU_PERC, color="gray", ls=":", label="2D perc 4/3")
    ax[1].set_xlabel(r"$\nu$"); ax[1].set_ylabel("collapse spread")
    ax[1].set_title("collapse quality vs exponent")
    ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(FIG + "fig_p3_fss_rho.png", dpi=150)
    plt.close(fig)

    print("  P_perc(rho) at f=%.2f:" % frac(Ns_fss))
    for k, rho in enumerate(rho_g):
        print(f"    rho={rho:.3f}: " +
              "  ".join(f"L{S}={Pperc[S][k]:.2f}" for S in Ls))
    print("figures -> fig_p3_rhoc.png, fig_p3_fss_rho.png")
    # stash key numbers for the report
    print("SUMMARY rho_c(f):", dict(zip([round(f, 2) for f in fsN],
                                        [round(x, 3) for x in rho_c])))
    print(f"SUMMARY nu_best={nu_best:.2f} rho_c_perc={rcp:.3f}")


if __name__ == "__main__":
    main()
