"""Paper 3 — P3-F/F2: sharpened independence test  r(d>=2) == 0  (torus).

The factorization law predicts the conversion-correlation coefficient r(d) = 0
EXACTLY for lattice separation d >= 2 (not merely "small"): two contrarians at
distance >= 2 share no witness. The only permitted term is the nearest-neighbour
r(1) (a shared bond). This script tests that prediction directly, with:

  (1) PERIODIC separation  d_eff = min(|dy|,L-|dy|) + min(|dx|,L-|dx|)  (torus
      minimum image), consistent with the periodic build_square. Pairs that were
      spuriously "far" under an open-lattice Manhattan distance fold back.
  (2) BOOTSTRAP error bars on r(d): resample the disorder realizations with
      replacement and recompute r(d); the spread is the 1-sigma bar.
  (3) An explicit PER-DISTANCE VERDICT: r(1) > 0 significant (shared bond);
      r(2), r(3) compatible with 0 at 2 sigma  ->  PASS.

Also regenerates the companion E3 panel (independent-site field reproduces the
measured P_perc) so the manuscript figure fig_p3_independence_v2 carries both the
sharpened E1 and E3. Data are cached to independence_v2.npz.

Run:  python3 run_independence_v2.py    (from src/lattice, periodic engine)
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lattice import torus_percolation
from run_factorization import arrows, frac        # periodic engine, prescription F

FIG = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "paper3", "")
DATA = os.path.join(os.path.dirname(__file__), "independence_v2.npz")

# ---- CONFIG ----
S = 8
NS_E1 = 6                 # near the transition: enough converted contrarians
RHO_E1 = 0.5
R1 = 200                  # realizations for E1
NBOOT = 400               # bootstrap resamples
DMAX = 4                  # report d_eff up to L/2 (periodic-image convention)
RHOS_E3 = np.round([0.35, 0.42, 0.50, 0.58, 0.66, 0.72], 3)
R3 = 80
SEED0_E1, SEED0_CONV, SEED0_E3, SEED_SYN = 8000, 500, 20000, 12321


def d_eff(p, q, L):
    """Torus minimum-image Manhattan separation between grid points p=(y,x)."""
    dy = abs(int(p[0]) - int(q[0])); dx = abs(int(p[1]) - int(q[1]))
    return min(dy, L - dy) + min(dx, L - dx)


def collect_E1():
    """Per realization: contrarian positions and converted (tau=+1) flags."""
    reals = []
    for r in range(R1):
        tau, cm = arrows(S, NS_E1, RHO_E1, SEED0_E1 + r)
        ys, xs = np.where(cm)
        pts = np.stack([ys, xs], 1)
        c = (tau[cm] > 0).astype(float)
        reals.append((pts, c))
    return reals


def rd_from(reals, dmax):
    """r(d) for d=1..dmax from a list of realizations, torus d_eff. Pooled cbar,
    var; correlation coefficient per distance."""
    c_all = np.concatenate([c for _, c in reals]) if reals else np.array([0.0])
    cbar = float(c_all.mean()); var = float(c_all.var())
    num = {d: [] for d in range(1, dmax + 1)}
    for pts, c in reals:
        n = len(c)
        for i in range(n):
            for j in range(i + 1, n):
                d = d_eff(pts[i], pts[j], S)
                if 1 <= d <= dmax:
                    num[d].append((c[i] - cbar) * (c[j] - cbar))
    r = np.array([np.mean(num[d]) / var if (var > 0 and num[d]) else np.nan
                  for d in range(1, dmax + 1)])
    npair = np.array([len(num[d]) for d in range(1, dmax + 1)])
    return r, npair, cbar, var


def bootstrap_rd(reals, dmax, nboot):
    rng = np.random.default_rng(2024)
    R = len(reals)
    boots = np.full((nboot, dmax), np.nan)
    for b in range(nboot):
        idx = rng.integers(0, R, R)
        rb, _, _, _ = rd_from([reals[k] for k in idx], dmax)
        boots[b] = rb
    return np.nanstd(boots, axis=0)


def measure_conv(Ns):
    """rho-averaged single-ring contrarian conversion probability at sharing Ns."""
    vals = []
    for rho in RHOS_E3:
        for r in range(25):
            tau, cm = arrows(S, Ns, rho, SEED0_CONV + r)
            if cm.any():
                vals.append((tau[cm] > 0).mean())
    return float(np.mean(vals))


def E3_curves():
    """quantum vs synthetic-independent P_perc(rho) for Ns in {4,6}."""
    rng_s = np.random.default_rng(SEED_SYN)
    out = {}
    for Ns in (4, 6):
        conv = measure_conv(Ns)
        q, s = [], []
        for rho in RHOS_E3:
            pq, ps = [], []
            for r in range(R3):
                tau, cm = arrows(S, Ns, rho, SEED0_E3 + r)
                pq.append(torus_percolation(tau > 0)[0])
                syn = np.ones_like(cm, float)
                syn[cm & (rng_s.random(cm.shape) >= conv)] = -1
                ps.append(torus_percolation(syn > 0)[0])
            q.append(np.mean(pq)); s.append(np.mean(ps))
        out[Ns] = (conv, np.array(q), np.array(s))
        print(f"  E3 Ns={Ns} f={frac(Ns):.2f} conv={conv:.3f} done")
    return out


def main():
    print(f"=== P3-F/F2 sharpened independence (torus, S={S}, f={frac(NS_E1):.2f}) ===")
    reals = collect_E1()
    r, npair, cbar, var = rd_from(reals, DMAX)
    err = bootstrap_rd(reals, DMAX, NBOOT)
    print(f"  conv={cbar:.3f}, var={var:.3f}  ({R1} realizations)")
    ds = np.arange(1, DMAX + 1)
    for k, d in enumerate(ds):
        z = r[k] / err[k] if err[k] > 0 else np.nan
        print(f"  r(d={d}) = {r[k]:+.3f} +/- {err[k]:.3f}  "
              f"({npair[k]} pairs, {z:+.1f} sigma)")
    # verdicts
    r1_sig = (r[0] > 0) and (r[0] > 2 * err[0])
    far_zero = all(abs(r[k]) <= 2 * err[k] for k in range(1, DMAX))  # d>=2
    verdict = r1_sig and far_zero
    print(f"  VERDICT: r(1)>0 significant = {r1_sig}; "
          f"r(d>=2) within 2 sigma of 0 = {far_zero}  ->  "
          f"{'PASS' if verdict else 'CHECK'}")

    e3 = E3_curves()
    e3mean = np.mean([np.mean(np.abs(q - s)) for _, q, s in e3.values()])
    floor = np.sqrt(0.25 / R3)
    print(f"  E3 mean|Delta| = {e3mean:.3f}  (binomial floor ~{floor:.3f})")

    np.savez(DATA, ds=ds, r=r, err=err, npair=npair, cbar=cbar, var=var,
             **{f"E3_{Ns}_q": e3[Ns][1] for Ns in e3},
             **{f"E3_{Ns}_s": e3[Ns][2] for Ns in e3},
             conv4=e3[4][0], conv6=e3[6][0], rhos_e3=RHOS_E3)

    # ---- figure: neutral titles ----
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    ax[0].axhline(0, color="k", lw=0.6)
    ax[0].errorbar(ds, r, yerr=err, fmt="o-", capsize=3, color="C0")
    ax[0].axvspan(1.5, DMAX + 0.5, alpha=0.08, color="C2")
    ax[0].set_xlabel(r"contrarian separation $d_{\rm eff}$ (torus min-image)")
    ax[0].set_ylabel(r"conversion correlation $r(d)$")
    ax[0].set_xticks(ds)
    ax[0].set_title("Conversion correlations: only the shared bond survives")
    ax[0].annotate(f"$r(1)={r[0]:+.2f}$\n(shared bond)", (1, r[0]),
                   textcoords="offset points", xytext=(12, 6), fontsize=8)
    ax[0].annotate(r"$r(d\geq 2)\approx 0$ (independent)",
                   (2.4, 0), textcoords="offset points", xytext=(0, 14),
                   fontsize=8, color="C2")
    for Ns, mk, c in zip((4, 6), ("o", "s"), ("C0", "C1")):
        conv, q, s = e3[Ns]
        ax[1].plot(RHOS_E3, q, mk + "-", color=c, label=f"quantum f={frac(Ns):.2f}")
        ax[1].plot(RHOS_E3, s, mk + "--", color="gray",
                   label=f"independent f={frac(Ns):.2f}")
    ax[1].set_xlabel(r"$\rho$"); ax[1].set_ylabel(r"$P_{\rm perc}$")
    ax[1].set_title(f"Independent-site field reproduces $P_{{\\rm perc}}$ "
                    f"(mean $\\Delta$={e3mean:.3f})")
    ax[1].legend(fontsize=7)
    fig.suptitle("Conversion correlations and the independent-field test "
                 "(torus, prescription F)")
    fig.tight_layout(); fig.savefig(FIG + "fig_p3_independence_v2.png", dpi=150)
    plt.close(fig)
    print("figure -> fig_p3_independence_v2.png ; data -> independence_v2.npz")
    print(f"SUMMARY r1={r[0]:+.3f}+/-{err[0]:.3f} r2={r[1]:+.3f}+/-{err[1]:.3f} "
          f"r3={r[2]:+.3f}+/-{err[2]:.3f} verdict={'PASS' if verdict else 'CHECK'} "
          f"e3mean={e3mean:.3f}")


if __name__ == "__main__":
    main()
