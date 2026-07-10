"""Paper 3 — P3-F/F3: replot fig_p3_convergence.png from cached data with a
NEUTRAL title (no headline exponent in the title or panel legends). Reads
convergence_data.npz produced by run_convergence.py; no resampling.

Run:  python3 replot_convergence.py    (from src/lattice)
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_convergence import analyze, NU_P, DATA, FIG

d = np.load(DATA)
LS = [int(s) for s in d["LS"]]
out = {S: {k: d[f"L{S}_{k}"] for k in ("dens", "Pp", "Pi", "Sm")} for S in LS}
a = analyze(out)

fig, ax = plt.subplots(1, 3, figsize=(14, 4.3))
for S in LS:
    ax[0].plot(out[S]["dens"], out[S]["Pp"], "o-", ms=4, label=f"L={S}")
ax[0].axhline(0.5, color="k", lw=0.4)
ax[0].set_xlabel(r"exact density $p=n_c/L^2$"); ax[0].set_ylabel(r"$P_{\rm perc}$")
ax[0].set_title(r"(a) spanning probability (crossing $\to \rho_c$)")
ax[0].legend(fontsize=8)
xl = np.log(np.array(LS, float))
ax[1].plot(xl, np.log([a["slope"][S] for S in LS]), "o-", label="slope estimator")
ax[1].plot(xl, xl / NU_P + (np.log(a["slope"][LS[0]]) - xl[0] / NU_P), "k:",
           label=r"percolation slope $1/\nu=3/4$")
ax[1].set_xlabel("ln L"); ax[1].set_ylabel("ln |slope|")
ax[1].set_title("(b) slope scaling at the crossing"); ax[1].legend(fontsize=8)
if len(a["nu_eff"]):
    ax[2].plot(a["inv_L"], a["nu_eff"], "o-", label=r"$\nu_{\rm eff}$ (pairwise)")
ax[2].axhline(NU_P, color="k", ls=":", label="percolation 4/3")
ax[2].set_xlabel("1/L"); ax[2].set_ylabel(r"$\nu_{\rm eff}$")
ax[2].set_xlim(left=0)
ax[2].set_title(r"(c) pairwise $\nu_{\rm eff}$ — trend, not extrapolation")
ax[2].legend(fontsize=8)
fig.suptitle("Exponent consistency checks (slope method, exact density; "
             "periodic boundaries)")
fig.tight_layout(); fig.savefig(FIG + "fig_p3_convergence.png", dpi=150)
plt.close(fig)
print(f"replotted fig_p3_convergence.png (neutral title); "
      f"nu={a['nu']:.2f} beta/nu={a['bn']:.3f} gamma/nu={a['gn']:.2f} "
      f"(shown on axes, not in title)")
