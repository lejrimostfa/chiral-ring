"""Paper 3 — P3-B: arrow contagion on the 1D chain (collision protocol).

Sweep L in {8,16,32}, sharing fraction f = Ns/(Np+2Ns) via Ns in {0..6},
contrarian density rho in {0.1..0.5}, R=10 disorder realizations, n_traj=2e4.
Arrow of a ring: tau_i = sign(<sigma_i^local>) in the collision protocol (the
signed accumulation rate; notes/paper3/arrow_definition.md). Continuous test:
A(f=0) = 1 - 2 rho exactly.

Outputs (printed + two figures):
  A_signed(f, rho), fraction of contrarians converted, mean tau-domain length,
  correlation length xi(f).
Figures: figures/paper3/fig_p3_chain_alignment.png, fig_p3_domains_chain.png

Run:  python3 run_chain.py    (from src/lattice)
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lattice import build_chain, draw_signs, sample_records, local_marginal_sigma

FIG = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "paper3", "")
NP = 3
DT = 0.3
T_CHARGE = 1.0
NTRAJ = 20_000
R = 10
LS = [8, 16, 32]
NSS = list(range(7))
RHOS = [0.1, 0.2, 0.3, 0.4, 0.5]


def frac(Ns):
    return Ns / (NP + 2 * Ns)


def chain_arrows(L, Ns, rho, seed):
    """tau (L,) and the contrarian index list, one disorder realization."""
    rng = np.random.default_rng(seed)
    n_c = round(rho * L)
    contr = tuple(np.sort(rng.choice(L, size=n_c, replace=False))) if n_c else ()
    G, Tc, owner = build_chain(L, NP, Ns, seed=seed, contrarians=contr,
                               T_charge=T_CHARGE)
    S_true = np.ones(L, int)
    beta_true = np.sin(2 * (S_true @ G) * (DT - Tc))
    Rrec = sample_records(beta_true, NTRAJ, rng, draw_sign=draw_signs(Tc))
    tau = np.array([np.sign(local_marginal_sigma(Rrec, G, Tc, DT, S_true, i).mean())
                    for i in range(L)])
    return tau, np.array(contr, int)


def domain_lengths(tau):
    """Lengths of maximal same-sign runs along the chain."""
    edges = np.where(np.diff(tau) != 0)[0] + 1
    return np.diff(np.concatenate(([0], edges, [len(tau)])))


def main():
    print("=== P3-B chain contagion (collision, Np=%d, dt=%.1f, T=%.1f) ===" %
          (NP, DT, T_CHARGE))
    # data[(L,rho)][Ns] -> dict of arrays over realizations
    A = {(L, rho): np.zeros((len(NSS), R)) for L in LS for rho in RHOS}
    conv = {(L, rho): np.zeros((len(NSS), R)) for L in LS for rho in RHOS}
    dom = {(L, rho): np.zeros((len(NSS), R)) for L in LS for rho in RHOS}
    for L in LS:
        for rho in RHOS:
            for k, Ns in enumerate(NSS):
                for r in range(R):
                    tau, contr = chain_arrows(L, Ns, rho, 6000 + r)
                    A[(L, rho)][k, r] = tau.mean()
                    conv[(L, rho)][k, r] = (tau[contr] > 0).mean() if len(contr) else 1.0
                    dom[(L, rho)][k, r] = domain_lengths(tau).mean()
        print(f"  L={L} done")

    # continuous baseline test: A(f=0) = 1 - 2 * (contrarian count)/L exactly.
    # (The reference is the DISCRETIZED density round(rho*L)/L, not rho itself.)
    print("  baseline test A(Ns=0) vs 1 - 2*round(rho*L)/L (L=32):")
    worst = 0.0
    for rho in RHOS:
        ref = 1 - 2 * round(rho * 32) / 32
        a0 = A[(32, rho)][0].mean()
        worst = max(worst, abs(a0 - ref))
        print(f"    rho={rho}: A={a0:+.3f}  ref={ref:+.3f}")
    print(f"  -> max deviation {worst:.2e} ({'PASS' if worst < 1e-9 else 'CHECK'})")

    fs = np.array([frac(Ns) for Ns in NSS])
    # per-ring convention (Eq. 1): in 1D (bulk degree 2) f_ring = 2 * f_bond.
    # All axes are in the per-ring convention to match the 2D figures.
    fs_ring = 2 * fs

    # ---------- figure 1: alignment + conversion (L=32) ----------
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    for rho in RHOS:
        m, s = A[(32, rho)].mean(1), A[(32, rho)].std(1)
        ax[0].errorbar(fs_ring, m, yerr=s, marker="o", ms=4, capsize=2,
                       label=f"$\\rho$={rho}")
        ax[0].axhline(1 - 2 * rho, color="gray", ls=":", lw=0.6)
        ax[1].errorbar(fs_ring, conv[(32, rho)].mean(1), yerr=conv[(32, rho)].std(1),
                       marker="s", ms=4, capsize=2, label=f"$\\rho$={rho}")
    ax[0].axhline(1, color="green", ls=":", lw=0.8)
    ax[0].set_xlabel(r"sharing fraction $f=2N_s/(N_p+2N_s)$ (per-ring)")
    ax[0].set_ylabel(r"signed alignment $A=\langle\tau_i\rangle$")
    ax[0].set_title("Signed alignment on the chain (L=32)")
    ax[0].legend(fontsize=8)
    ax[1].set_xlabel(r"sharing fraction $f$ (per-ring)")
    ax[1].set_ylabel("fraction of contrarians converted")
    ax[1].set_title("Contrarian conversion versus sharing (L=32)")
    ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(FIG + "fig_p3_chain_alignment.png", dpi=150)
    plt.close(fig)

    # ---------- figure 2: domains + correlation length ----------
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    for rho in RHOS:
        ax[0].errorbar(fs_ring, dom[(32, rho)].mean(1), yerr=dom[(32, rho)].std(1),
                       marker="o", ms=4, capsize=2, label=f"$\\rho$={rho}")
    ax[0].set_xlabel(r"sharing fraction $f$ (per-ring)")
    ax[0].set_ylabel(r"mean $\tau$-domain length")
    ax[0].set_title("Mean domain length versus sharing (L=32)")
    ax[0].legend(fontsize=8)
    for L in LS:
        ax[1].errorbar(fs_ring, dom[(L, 0.3)].mean(1), yerr=dom[(L, 0.3)].std(1),
                       marker="^", ms=5, capsize=2, label=f"L={L}")
    ax[1].set_xlabel(r"sharing fraction $f$ (per-ring)")
    ax[1].set_ylabel(r"mean $\tau$-domain length  [$\rho$=0.3]")
    ax[1].set_title(r"Domain length versus size ($L=8/16/32$, $\rho$=0.3)")
    ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(FIG + "fig_p3_domains_chain.png", dpi=150)
    plt.close(fig)

    # ---------- printed summary ----------
    print("  A_signed(f) at rho=0.3, L=32 (domain length = 1D correlation scale):")
    for k, Ns in enumerate(NSS):
        print(f"    f={fs[k]:.2f}: A={A[(32,0.3)][k].mean():+.3f}  "
              f"conv={conv[(32,0.3)][k].mean():.2f}  "
              f"domain_len={dom[(32,0.3)][k].mean():.1f}")
    print("figures -> fig_p3_chain_alignment.png, fig_p3_domains_chain.png")


if __name__ == "__main__":
    main()
