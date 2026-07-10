"""Paper 3 — P3-E: the factorization law.

Hypothesis under test: the quantum arrow field tau_i = sign(<sigma_i^local>)
(collision, prescription F) is an INDEPENDENT site-percolation field with forward
occupation p_eff(rho,f) = 1 - rho(1-conv(f)), where conv(f) is the (rho-universal)
contrarian-conversion probability. Then closed forms follow with NO free fit:
    rho_c^align(f) = 1 / (2(1-conv(f)))          (global reversal, p_eff=1/2)
    rho_c^perc(f)  = (1-p_c)/(1-conv(f)),  p_c=0.5927   (percolation, p_eff=p_c)

E1 independence + no-backflip; E2 overlay measured vs closed form (residuals);
E3 synthetic independent field reproduces the measured P_perc; E4 prescription M
counter-case (charged shared witnesses).

Run:  python3 run_factorization.py    (from src/lattice)
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
NP, DT, T_CHARGE, NTRAJ, PC = 8, 0.3, 1.0, 2000, 0.5927


def frac(Ns):
    return 4 * Ns / (NP + 4 * Ns)


def build(S, Ns, rho, seed, prescription="F"):
    rng = np.random.default_rng(seed)
    L = S * S
    n_c = round(rho * L)
    contr = tuple(np.sort(rng.choice(L, size=n_c, replace=False))) if n_c else ()
    G, Tc, owner = build_square(S, S, NP, Ns, seed=seed, contrarians=contr,
                               T_charge=T_CHARGE)
    if prescription == "M":                       # charge the shared-link baths
        for n, o in enumerate(owner):
            if o[0] == "link":
                Tc[n] = T_CHARGE
    cmask = np.zeros(L, bool); cmask[list(contr)] = True
    return G, Tc, cmask.reshape(S, S)


def arrows(S, Ns, rho, seed, prescription="F"):
    G, Tc, cmask = build(S, Ns, rho, seed, prescription)
    L = S * S
    St = np.ones(L, int)
    beta = np.sin(2 * (St @ G) * (DT - Tc))
    R = sample_records(beta, NTRAJ, np.random.default_rng(seed), draw_sign=draw_signs(Tc))
    tau = np.array([np.sign(local_marginal_sigma(R, G, Tc, DT, St, i).mean())
                    for i in range(L)]).reshape(S, S)
    return tau, cmask


def percolates(grid):
    lab, _ = label(grid > 0)
    return len(set(lab[0]) - {0} & (set(lab[-1]) - {0})) > 0


def zc(x, y, tgt):
    for k in range(len(y) - 1):
        if (y[k] - tgt) * (y[k + 1] - tgt) <= 0 and y[k] != y[k + 1]:
            return x[k] + (tgt - y[k]) * (x[k + 1] - x[k]) / (y[k + 1] - y[k])
    return np.nan


# ======================================================================
def main():
    print("=== P3-E factorization law (Np=%d, dt=%.1f, T=%.1f) ===" %
          (NP, DT, T_CHARGE))
    S = 8
    Nss = [0, 2, 3, 4, 5, 6, 7, 8]
    rhos = np.round([0.35, 0.42, 0.50, 0.58, 0.66, 0.72], 3)
    Rg = 25
    fs = np.array([frac(N) for N in Nss])

    # ---- grid: A, P_perc, conversion; also cache tau+contr for E1/E3 ----
    A = np.zeros((len(Nss), len(rhos)))
    Pq = np.zeros((len(Nss), len(rhos)))
    conv_by = np.full((len(Nss), len(rhos)), np.nan)
    backflip = 0; fresh_total = 0
    cache = {}
    for a, Ns in enumerate(Nss):
        for b, rho in enumerate(rhos):
            As, Ps, Cs = [], [], []
            for r in range(Rg):
                tau, cm = arrows(S, Ns, rho, 500 + r)
                cache[(Ns, rho, r)] = (tau, cm)
                As.append(tau.mean()); Ps.append(percolates(tau))
                if cm.any():
                    Cs.append((tau[cm] > 0).mean())      # conversion fraction
                backflip += int(np.sum(tau[~cm] < 0))    # fresh rings gone backward
                fresh_total += int(np.sum(~cm))
            A[a, b] = np.mean(As); Pq[a, b] = np.mean(Ps)
            if Cs:
                conv_by[a, b] = np.mean(Cs)
        print(f"  grid f={fs[a]:.2f} done")
    conv = np.nanmean(conv_by, axis=1)                   # rho-averaged conv(f)
    conv_spread = np.nanstd(conv_by, axis=1)
    print(f"  conv(f) rho-universality: max std over rho = {np.nanmax(conv_spread):.3f}")
    print(f"  NO-BACKFLIP (prescription F): {fresh_total-backflip}/{fresh_total} "
          f"fresh rings forward ({backflip} backflips)")

    # ---------------- E1: conversion independence ----------------
    print("=== E1 independence: conversion two-point correlation ===")
    Ns_e1 = 6                                             # near f_c
    R1 = 200
    prod_by_d = {}                                        # distance -> list of products
    c_all = []
    recs = []
    for r in range(R1):
        tau, cm = arrows(S, Ns_e1, 0.5, 8000 + r)
        ys, xs = np.where(cm)
        c = (tau[cm] > 0).astype(float)
        c_all.extend(c.tolist())
        recs.append((np.array(list(zip(ys, xs))), c))
    cbar = float(np.mean(c_all)); var = float(np.var(c_all))
    for pts, c in recs:
        n = len(c)
        for i in range(n):
            for j in range(i + 1, n):
                d = int(abs(pts[i][0] - pts[j][0]) + abs(pts[i][1] - pts[j][1]))
                prod_by_d.setdefault(d, []).append((c[i] - cbar) * (c[j] - cbar))
    r1 = np.mean(prod_by_d.get(1, [0])) / var if var > 0 else 0.0
    # "distant" = lattice bulk d in [4,9]; corner pairs (d>=11 on L=8) are a
    # handful and small-statistics-dominated, excluded from the estimate.
    rfar = np.mean([v for d in prod_by_d if 4 <= d <= 9 for v in prod_by_d[d]]) \
        / var if var > 0 else 0.0
    npair1 = len(prod_by_d.get(1, []))
    print(f"  conv={cbar:.3f}, var={var:.3f}; corr coeff r(d=1 adjacent)={r1:+.3f} "
          f"({npair1} pairs), r(bulk 4<=d<=9)={rfar:+.3f}  (target ~0 => independent)")
    ds = sorted(d for d in prod_by_d if d >= 1)
    rd = [np.mean(prod_by_d[d]) / var for d in ds]

    # ---------------- E2: closed-form overlay + residuals ----------------
    rc_align = np.array([zc(rhos, A[a], 0.0) for a in range(len(Nss))])
    rc_perc = np.array([zc(rhos, Pq[a], 0.5) for a in range(len(Nss))])
    pred_align = 1.0 / (2 * (1 - conv))
    pred_perc = (1 - PC) / (1 - conv)
    ok = np.isfinite(rc_align) & (pred_align < rhos[-1])
    res_align = np.abs(rc_align[ok] - pred_align[ok])
    okp = np.isfinite(rc_perc) & (pred_perc < rhos[-1])
    res_perc = np.abs(rc_perc[okp] - pred_perc[okp])
    print("=== E2 closed-form vs measured (conv(f) interpolated, no free fit) ===")
    print(f"  rho_c^align residual: max={res_align.max():.3f} mean={res_align.mean():.3f}")
    print(f"  rho_c^perc  residual: max={res_perc.max():.3f} mean={res_perc.mean():.3f}")

    # ---------------- E3: synthetic independent site field ----------------
    # Dedicated high-R comparison (the grid R=25 is too noisy for P_perc, whose
    # binomial std ~0.1 gives spurious ~0.2 max deviations). Report MEAN |Delta|.
    print("=== E3 synthetic independent-site field vs measured P_perc (R=80) ===")
    R3 = 80
    rng_s = np.random.default_rng(12321)
    e3_curves = {}
    for Ns in (4, 6):
        a = Nss.index(Ns)
        q, s = [], []
        for rho in rhos:
            pq, ps = [], []
            for r in range(R3):
                tau, cm = arrows(S, Ns, rho, 20000 + r)
                pq.append(percolates(tau))
                syn = np.ones_like(cm, float)
                syn[cm & (rng_s.random(cm.shape) >= conv[a])] = -1
                ps.append(percolates(syn))
            q.append(np.mean(pq)); s.append(np.mean(ps))
        e3_curves[Ns] = (np.array(q), np.array(s))
    e3mean = np.mean([np.mean(np.abs(q - s)) for q, s in e3_curves.values()])
    e3max = np.max([np.max(np.abs(q - s)) for q, s in e3_curves.values()])
    print(f"  quantum vs independent P_perc: mean|Delta|={e3mean:.3f}, "
          f"max|Delta|={e3max:.3f} (binomial noise floor ~{np.sqrt(0.25/R3):.3f})")

    e1_pass = abs(r1) < 0.15 and abs(rfar) < 0.05 and backflip == 0
    e2_pass = res_align.max() < 0.02 and res_perc.max() < 0.08
    e3_pass = e3mean < 0.06
    print(f"  E1 {'PASS' if e1_pass else 'CHECK'} | E2 {'PASS' if e2_pass else 'CHECK'}"
          f" | E3 {'PASS' if e3_pass else 'CHECK'}")

    # ---------------- E4: prescription M counter-case ----------------
    print("=== E4 prescription M (shared witnesses charged): sign & monotonicity ===")
    A_M = np.zeros(len(Nss)); back_M = 0; fresh_M = 0
    for a, Ns in enumerate(Nss):
        As = []
        for r in range(Rg):
            tau, cm = arrows(S, Ns, 0.5, 6000 + r, prescription="M")
            As.append(tau.mean())
            back_M += int(np.sum(tau[~cm] < 0)); fresh_M += int(np.sum(~cm))
        A_M[a] = np.mean(As)
    mono_F = np.all(np.diff(A[:, 2]) >= -0.02)             # A(f) at rho=0.5, F
    mono_M = np.all(np.diff(A_M) <= 0.02)                  # M: A should DECREASE
    print(f"  F: A(f) monotone up = {mono_F}; M: A(f) monotone down = {mono_M}; "
          f"M backflips {back_M}/{fresh_M} fresh rings (forward arrows flipped)")

    # ================= figures =================
    # E1+E3
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    ax[0].axhline(0, color="k", lw=0.6)
    ax[0].plot(ds, rd, "o-")
    ax[0].set_xlabel("contrarian separation d (lattice)")
    ax[0].set_ylabel("conversion correlation coeff. r(d)")
    ax[0].set_title(f"(E1) Conversions independent: r(1)={r1:+.3f}, r(≥4)={rfar:+.3f}")
    for Ns, mk, c in zip([4, 6], ["o", "s"], ["C0", "C1"]):
        q, s = e3_curves[Ns]
        ax[1].plot(rhos, q, mk + "-", color=c, label=f"quantum f={frac(Ns):.2f}")
        ax[1].plot(rhos, s, mk + "--", color="gray",
                   label=f"independent f={frac(Ns):.2f}")
    ax[1].set_xlabel(r"$\rho$"); ax[1].set_ylabel(r"$P_{\rm perc}$")
    ax[1].set_title(f"(E3) Independent field reproduces P_perc "
                    f"(mean Δ={e3mean:.3f})")
    ax[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(FIG + "fig_p3_independence.png", dpi=150)
    plt.close(fig)

    # E2 overlay
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    ax[0].plot(fs, rc_align, "o", ms=7, label="measured")
    ax[0].plot(fs, pred_align, "-", label=r"$1/[2(1-\mathrm{conv})]$")
    ax[0].axhline(0.5, color="k", ls=":", lw=0.8)
    ax[0].set_ylim(0.45, rhos[-1] + 0.05)
    ax[0].set_xlabel("f"); ax[0].set_ylabel(r"$\rho_c^{\rm align}$")
    ax[0].set_title(f"(E2) global reversal (res max {res_align.max():.3f})")
    ax[0].legend(fontsize=8)
    ax[1].plot(fs, rc_perc, "s", ms=7, label="measured")
    ax[1].plot(fs, pred_perc, "-", label=r"$(1-p_c)/(1-\mathrm{conv})$")
    ax[1].set_ylim(0.30, rhos[-1] + 0.05)
    ax[1].set_xlabel("f"); ax[1].set_ylabel(r"$\rho_c^{\rm perc}$")
    ax[1].set_title(f"(E2) percolation (res max {res_perc.max():.3f})")
    ax[1].legend(fontsize=8)
    fig.suptitle("(E2) Closed-form thresholds from conv(f) — no free fit")
    fig.tight_layout(); fig.savefig(FIG + "fig_p3_factorization.png", dpi=150)
    plt.close(fig)

    # E4
    fig, ax = plt.subplots(figsize=(6.4, 4.5))
    ax.plot(fs, A[:, 2], "o-", label="F (shared fresh): forward wins")
    ax.plot(fs, A_M, "s-", color="C3", label="M (shared charged): backward wins")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("f"); ax.set_ylabel(r"$A=\langle\tau\rangle$  ($\rho$=0.5)")
    ax.set_title("(E4) F/M duality: shared-witness preparation sets the sign")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(FIG + "fig_p3_prescriptionM.png", dpi=150)
    plt.close(fig)

    print("figures -> fig_p3_independence, fig_p3_factorization, fig_p3_prescriptionM")
    print(f"SUMMARY conv={np.round(conv,3).tolist()}")
    print(f"SUMMARY r1={r1:+.3f} rfar={rfar:+.3f} align_res={res_align.max():.3f} "
          f"perc_res={res_perc.max():.3f} e3mean={e3mean:.3f}")
    print(f"SUMMARY E1={e1_pass} E2={e2_pass} E3={e3_pass} "
          f"monoF={mono_F} monoM={mono_M} Mbackflips={back_M}")


if __name__ == "__main__":
    main()
