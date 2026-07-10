"""Paper 3 — reinforcement: universality-class check (nu f-independence + beta/gamma).

Strengthens the two weak points of the exponent claim:
  (1) nu measured at TWO sharing fractions f in {0.5, 0.71} (not one), L in {8,12,16},
      via the rho-driven percolation collapse -> is nu ~ 4/3 f-independent?
  (2) beta/nu and gamma/nu by finite-size scaling AT criticality, from the SAME
      tau fields: P_inf(rho_c, L) ~ L^{-beta/nu}; S_max(L) ~ L^{gamma/nu}.
2D site percolation: nu=4/3, beta=5/36, gamma=43/18 => beta/nu=5/48=0.104,
gamma/nu=43/24=1.79.

HONEST SCOPE: L<=16 gives order-of-magnitude exponents, not 3-digit precision
(that needs L>=32 / HPC). This tests CONSISTENCY with 2D percolation.

Run:  python3 run_universality.py    (from src/lattice)
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
NP, DT, T_CHARGE, NTRAJ = 8, 0.3, 1.0, 1500
NU_P, BN_P, GN_P = 4.0 / 3.0, 5.0 / 48.0, 43.0 / 24.0   # 2D percolation refs
LS = [8, 12, 16]
R = 40


def frac(Ns):
    return 4 * Ns / (NP + 4 * Ns)


def tau_grid(S, Ns, rho, seed):
    rng = np.random.default_rng(seed)
    L = S * S
    n_c = round(rho * L)
    contr = tuple(np.sort(rng.choice(L, size=n_c, replace=False))) if n_c else ()
    G, Tc, _ = build_square(S, S, NP, Ns, seed=seed, contrarians=contr,
                            T_charge=T_CHARGE)
    St = np.ones(L, int)
    beta = np.sin(2 * (St @ G) * (DT - Tc))
    Rr = sample_records(beta, NTRAJ, rng, draw_sign=draw_signs(Tc))
    tau = np.array([np.sign(local_marginal_sigma(Rr, G, Tc, DT, St, i).mean())
                    for i in range(L)]).reshape(S, S)
    return tau


def stats(grid):
    """(perc, Pinf, S) with torus (periodic-BC) connectivity and the wrapping
    percolation criterion, consistent with the periodic build_square."""
    return torus_percolation(grid > 0)


def zc(x, y, tgt):
    for k in range(len(y) - 1):
        if (y[k] - tgt) * (y[k + 1] - tgt) <= 0 and y[k] != y[k + 1]:
            return x[k] + (tgt - y[k]) * (x[k + 1] - x[k]) / (y[k + 1] - y[k])
    return np.nan


def collapse_nu(fs_rho, Pp, rc):
    def spread(nu):
        xg = np.linspace(-1.1, 1.1, 21)
        cur = []
        for S in LS:
            x = (fs_rho - rc) * S ** (1.0 / nu)
            o = np.argsort(x)
            cur.append(np.interp(xg, x[o], Pp[S][o], left=np.nan, right=np.nan))
        return np.nanmean(np.nanstd(np.array(cur), axis=0))
    nus = np.linspace(0.9, 1.9, 26)
    sp = [spread(nu) for nu in nus]
    return nus[int(np.argmin(sp))]


def main():
    print("=== Universality reinforcement (nu f-independence + beta/gamma) ===")
    cases = {0.5: 2, 0.71: 5}                       # f -> Ns
    rho_grids = {0.5: np.round(np.linspace(0.28, 0.50, 11), 4),
                 0.71: np.round(np.linspace(0.36, 0.58, 11), 4)}
    results = {}
    for f, Ns in cases.items():
        rhos = rho_grids[f]
        Pp = {S: np.zeros(len(rhos)) for S in LS}
        Pi = {S: np.zeros(len(rhos)) for S in LS}
        Sm = {S: np.zeros(len(rhos)) for S in LS}
        for S in LS:
            for k, rho in enumerate(rhos):
                pp, pi, sm = [], [], []
                for r in range(R):
                    per, PInf, Sus = stats(tau_grid(S, Ns, rho, 700 + r))
                    pp.append(per); pi.append(PInf); sm.append(Sus)
                Pp[S][k] = np.mean(pp); Pi[S][k] = np.mean(pi); Sm[S][k] = np.mean(sm)
            print(f"  f={f} L={S} done")
        rc = zc(rhos, Pp[max(LS)], 0.5)
        nu = collapse_nu(rhos, Pp, rc)
        # beta/nu: P_inf at rho_c vs L (interp each L to rho_c)
        Pinf_c = np.array([np.interp(rc, rhos, Pi[S]) for S in LS])
        bn = -np.polyfit(np.log(LS), np.log(np.clip(Pinf_c, 1e-6, None)), 1)[0]
        # gamma/nu: peak susceptibility vs L
        Smax = np.array([Sm[S].max() for S in LS])
        gn = np.polyfit(np.log(LS), np.log(Smax), 1)[0]
        results[f] = dict(rhos=rhos, Pp=Pp, Pi=Pi, Sm=Sm, rc=rc, nu=nu,
                          Pinf_c=Pinf_c, bn=bn, Smax=Smax, gn=gn)
        print(f"  f={f}: rho_c={rc:.3f}, nu={nu:.2f} (perc 4/3={NU_P:.2f}), "
              f"beta/nu={bn:.3f} (perc {BN_P:.3f}), gamma/nu={gn:.2f} (perc {GN_P:.2f})")

    # ---- figure ----
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.3))
    for f, c in zip(cases, ["C0", "C1"]):
        r = results[f]
        for S, mk in zip(LS, ["o", "s", "^"]):
            x = (r["rhos"] - r["rc"]) * S ** (1.0 / r["nu"])
            ax[0].plot(x, r["Pp"][S], mk + "-", color=c, ms=4, alpha=0.8,
                       label=f"f={f} L={S}" if S == LS[0] else None)
    ax[0].set_xlabel(r"$(\rho-\rho_c)L^{1/\nu}$"); ax[0].set_ylabel(r"$P_{\rm perc}$")
    ax[0].set_title(f"(a) nu collapse: f=0.5->{results[0.5]['nu']:.2f}, "
                    f"f=0.71->{results[0.71]['nu']:.2f} (4/3={NU_P:.2f})")
    ax[0].legend(fontsize=7)
    for f, c in zip(cases, ["C0", "C1"]):
        r = results[f]
        ax[1].plot(np.log(LS), np.log(r["Pinf_c"]), "o-", color=c,
                   label=f"f={f}: β/ν={r['bn']:.3f}")
    xl = np.log(LS)
    ax[1].plot(xl, xl * (-BN_P) + (np.log(results[0.5]['Pinf_c'][0]) + BN_P * xl[0]),
               "k:", label=f"perc slope -{BN_P:.3f}")
    ax[1].set_xlabel("ln L"); ax[1].set_ylabel(r"ln $P_\infty(\rho_c)$")
    ax[1].set_title("(b) order parameter: β/ν"); ax[1].legend(fontsize=7)
    for f, c in zip(cases, ["C0", "C1"]):
        r = results[f]
        ax[2].plot(np.log(LS), np.log(r["Smax"]), "o-", color=c,
                   label=f"f={f}: γ/ν={r['gn']:.2f}")
    ax[2].plot(xl, xl * GN_P + (np.log(results[0.5]['Smax'][0]) - GN_P * xl[0]),
               "k:", label=f"perc slope {GN_P:.2f}")
    ax[2].set_xlabel("ln L"); ax[2].set_ylabel(r"ln $S_{\max}$")
    ax[2].set_title("(c) susceptibility: γ/ν"); ax[2].legend(fontsize=7)
    fig.suptitle("Universality reinforcement — 2D percolation (nu=4/3, "
                 "beta/nu=0.104, gamma/nu=1.79); L<=16")
    fig.tight_layout(); fig.savefig(FIG + "fig_p3_universality.png", dpi=150)
    plt.close(fig)

    print("figures -> fig_p3_universality.png")
    for f in cases:
        r = results[f]
        print(f"SUMMARY f={f}: nu={r['nu']:.2f} beta/nu={r['bn']:.3f} "
              f"gamma/nu={r['gn']:.2f}")


if __name__ == "__main__":
    main()
