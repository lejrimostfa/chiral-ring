"""Paper 3 — exponent convergence / stability study.

Addresses the finite-size instability of nu found in run_universality.py
(nu_eff = 1.33..1.74 depending on the rho-grid). Three fixes:

  (1) ROBUST ESTIMATOR. Instead of minimizing a fragile collapse spread, use the
      slope method: at the percolation point, dP_perc/dp|_{p_c} ~ L^{1/nu}. The
      slope of a monotone P_perc(p) is well defined and low-variance.
  (2) DISCRETIZATION FIX. Sweep the INTEGER contrarian count n_c and use the
      exact density p = n_c/L^2 (no round(rho*L) ambiguity that misaligns curves
      across L).
  (3) L_max SWEEP + extrapolation. Effective exponents from consecutive sizes,
      plotted vs 1/L, extrapolated to L -> infinity.

Also beta/nu from P_inf(p_c,L) ~ L^{-beta/nu} and gamma/nu from S_max(L) ~
L^{gamma/nu}. 2D site percolation: nu=4/3, beta/nu=5/48=0.104, gamma/nu=43/24=1.79.

CONFIGURABLE (top constants). Reduced local run: LS=[12,16,24]. HPC: add 32,48
and raise R, NTRAJ. Raw data saved to an .npz so the analysis/figure can be
re-made without resampling.

Run:  python3 run_convergence.py    (from src/lattice)
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
DATA = os.path.join(os.path.dirname(__file__), "convergence_data.npz")

# ---- CONFIG (raise LS/R/NTRAJ for an HPC run) ----
NP, DT, T_CHARGE = 8, 0.3, 1.0
NS = 2                       # f = 4*2/(8+8) = 0.5
LS = [12, 16, 24]            # HPC: [16, 24, 32, 48]
R = 50
NTRAJ = 2000                 # HPC: 4000
P_LO, P_HI = 0.34, 0.46      # density window around p_c ~ 0.40
N_DENS = 9
NU_P, BN_P, GN_P = 4.0 / 3.0, 5.0 / 48.0, 43.0 / 24.0


def tau_grid(S, n_c, seed):
    rng = np.random.default_rng(seed)
    L = S * S
    contr = tuple(np.sort(rng.choice(L, size=n_c, replace=False))) if n_c else ()
    G, Tc, _ = build_square(S, S, NP, NS, seed=seed, contrarians=contr,
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


def measure():
    """P_perc, P_inf, S vs exact density p=n_c/L^2, per L. Returns dict."""
    out = {}
    for S in LS:
        L2 = S * S
        ncs = np.unique(np.round(np.linspace(P_LO, P_HI, N_DENS) * L2).astype(int))
        dens = ncs / L2
        Pp, Pi, Sm = [], [], []
        for nc in ncs:
            pp, pi, sm = [], [], []
            for r in range(R):
                per, PInf, Sus = stats(tau_grid(S, int(nc), 700 + r))
                pp.append(per); pi.append(PInf); sm.append(Sus)
            Pp.append(np.mean(pp)); Pi.append(np.mean(pi)); Sm.append(np.mean(sm))
        out[S] = dict(dens=dens, Pp=np.array(Pp), Pi=np.array(Pi), Sm=np.array(Sm))
        print(f"  L={S} done ({len(ncs)} densities x {R} real)")
    return out


def analyze(out):
    def cross(d, y, tgt):
        for k in range(len(y) - 1):
            if (y[k] - tgt) * (y[k + 1] - tgt) <= 0 and y[k] != y[k + 1]:
                return d[k] + (tgt - y[k]) * (d[k + 1] - d[k]) / (y[k + 1] - y[k])
        return np.nan
    pc, slope, Pinf_c, Smax = {}, {}, {}, {}
    for S in LS:
        d, Pp = out[S]["dens"], out[S]["Pp"]
        pc[S] = cross(d, Pp, 0.5)
        # slope from a linear fit over the steep region (0.15 < P < 0.85)
        m = (Pp > 0.12) & (Pp < 0.88)
        slope[S] = abs(np.polyfit(d[m], Pp[m], 1)[0]) if m.sum() >= 2 else np.nan
        Pinf_c[S] = np.interp(pc[S], d, out[S]["Pi"])
        Smax[S] = out[S]["Sm"].max()
    Ls = np.array(LS, float)
    inv_nu = np.polyfit(np.log(Ls), np.log([slope[S] for S in LS]), 1)[0]
    nu = 1.0 / inv_nu
    bn = -np.polyfit(np.log(Ls), np.log([Pinf_c[S] for S in LS]), 1)[0]
    gn = np.polyfit(np.log(Ls), np.log([Smax[S] for S in LS]), 1)[0]
    # effective nu from consecutive size pairs (convergence)
    nu_eff, inv_L = [], []
    for a in range(len(LS) - 1):
        L1, L2 = LS[a], LS[a + 1]
        inv_nu_e = np.log(slope[L2] / slope[L1]) / np.log(L2 / L1)
        nu_eff.append(1.0 / inv_nu_e); inv_L.append(2.0 / (L1 + L2))
    return dict(pc=pc, slope=slope, Pinf_c=Pinf_c, Smax=Smax, nu=nu, bn=bn, gn=gn,
                nu_eff=np.array(nu_eff), inv_L=np.array(inv_L))


def main():
    print(f"=== exponent convergence study (LS={LS}, R={R}, NTRAJ={NTRAJ}) ===")
    out = measure()
    np.savez(DATA, LS=LS, **{f"L{S}_{k}": out[S][k]
                             for S in LS for k in out[S]})
    a = analyze(out)
    print(f"  rho_c(L): " + ", ".join(f"L{S}={a['pc'][S]:.3f}" for S in LS))
    print(f"  slope(L): " + ", ".join(f"L{S}={a['slope'][S]:.2f}" for S in LS))
    print(f"  ALL-L fit: nu={a['nu']:.2f} (perc {NU_P:.2f}), "
          f"beta/nu={a['bn']:.3f} (perc {BN_P:.3f}), "
          f"gamma/nu={a['gn']:.2f} (perc {GN_P:.2f})")
    print(f"  nu_eff(consecutive pairs) = {np.round(a['nu_eff'],2).tolist()} "
          f"at 1/L = {np.round(a['inv_L'],3).tolist()}")

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.3))
    for S in LS:
        ax[0].plot(out[S]["dens"], out[S]["Pp"], "o-", ms=4, label=f"L={S}")
    ax[0].axhline(0.5, color="k", lw=0.4)
    ax[0].set_xlabel(r"exact density $p=n_c/L^2$"); ax[0].set_ylabel(r"$P_{\rm perc}$")
    ax[0].set_title("(a) P_perc vs exact density (crossing $\\to p_c$)")
    ax[0].legend(fontsize=8)
    xl = np.log(np.array(LS, float))
    ax[1].plot(xl, np.log([a["slope"][S] for S in LS]), "o-",
               label=f"slope $\\to 1/\\nu$, $\\nu$={a['nu']:.2f}")
    ax[1].plot(xl, xl / NU_P + (np.log(a["slope"][LS[0]]) - xl[0] / NU_P), "k:",
               label=r"perc $1/\nu$=0.75")
    ax[1].set_xlabel("ln L"); ax[1].set_ylabel("ln |slope|")
    ax[1].set_title("(b) slope method for $\\nu$"); ax[1].legend(fontsize=8)
    if len(a["nu_eff"]):
        ax[2].plot(a["inv_L"], a["nu_eff"], "o-", label=r"$\nu_{\rm eff}$ (pairs)")
    ax[2].axhline(NU_P, color="k", ls=":", label="perc 4/3")
    ax[2].set_xlabel("1/L"); ax[2].set_ylabel(r"$\nu_{\rm eff}$")
    ax[2].set_xlim(left=0)
    ax[2].set_title("(c) convergence: extrapolate 1/L→0"); ax[2].legend(fontsize=8)
    fig.suptitle(f"Exponent convergence (slope method, exact density; LS={LS}, "
                 f"R={R}) — beta/nu={a['bn']:.2f}, gamma/nu={a['gn']:.2f}")
    fig.tight_layout(); fig.savefig(FIG + "fig_p3_convergence.png", dpi=150)
    plt.close(fig)
    print("figures -> fig_p3_convergence.png ; data -> convergence_data.npz")
    print(f"SUMMARY nu={a['nu']:.2f} beta/nu={a['bn']:.3f} gamma/nu={a['gn']:.2f} "
          f"nu_eff={np.round(a['nu_eff'],2).tolist()}")


if __name__ == "__main__":
    main()
