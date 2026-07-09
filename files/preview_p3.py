"""Preview v2: arrow = direction of record ACCUMULATION.
tau_i = sign( d<sigma_i^local>/dt )  -- forward if ring i's accessible
information is growing, backward if being consumed. Contagion test."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lattice_p3 import *

L, Np, Ns, t0, dt, T_charge = 10, 3, 2, 0.8, 0.06, 1.6
contr = (2, 5, 8)
lams = np.linspace(0, 1.6, 9)
n_traj, seeds = 4000, (0, 1, 2)
S_all = all_branches(L); true_idx = 0

def mean_local_sigma(G, Tc, owner, t, seed):
    beta = beta_matrix(G, Tc, t, S_all)
    rng = np.random.default_rng(seed)
    R = sample_records(beta[true_idx], n_traj, rng)
    out = []
    for i in range(L):
        nbrs = [j for j in (i-1, i) if 0 <= j < L-1]
        sl = local_marginal_sigma(R, G, Tc, t, owner, S_all, true_idx, i, nbrs)
        out.append(sl.mean())
    return np.array(out)

curves, flips = [], []
for seed in seeds:
    A_signed, frac_flip = [], []
    for lam in lams:
        G, Tc, owner = build_chain(L, Np, Ns, seed=seed, contrarians=contr,
                                   T_charge=T_charge, shared_strength=lam)
        m_plus  = mean_local_sigma(G, Tc, owner, t0+dt, 100+seed)
        m_minus = mean_local_sigma(G, Tc, owner, t0-dt, 200+seed)
        tau = np.sign(m_plus - m_minus)          # arrow = accumulation direction
        A_signed.append(tau.mean())
        frac_flip.append((tau[list(contr)] > 0).mean())
    curves.append(A_signed); flips.append(frac_flip)

curves, flips = np.array(curves), np.array(flips)
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].errorbar(lams, curves.mean(0), yerr=curves.std(0), marker='o', capsize=3)
ax[0].axhline(1-2*len(contr)/L, ls=':', c='gray', label=r'no contagion: $1-2\rho$')
ax[0].axhline(1, ls=':', c='green', label='full alignment')
ax[0].set_xlabel(r'shared coupling strength $\lambda_s$')
ax[0].set_ylabel(r'signed alignment $\langle\tau_i\rangle$')
ax[0].set_title(f'Arrow contagion (chain L={L}, 30% contrarians)')
ax[0].legend(fontsize=8)
ax[1].errorbar(lams, flips.mean(0), yerr=flips.std(0), marker='s', color='crimson', capsize=3)
ax[1].set_xlabel(r'$\lambda_s$'); ax[1].set_ylabel('contrarians converted')
ax[1].set_title('Backward arrows flipped forward by shared records')
plt.tight_layout(); plt.savefig('fig_p3_preview.png', dpi=140)
print("lam:      ", np.round(lams,2).tolist())
print("A_signed: ", np.round(curves.mean(0),3).tolist())
print("converted:", np.round(flips.mean(0),3).tolist())
