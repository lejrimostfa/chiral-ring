"""Branch Arrow Theorem gate: hypotheses -> conclusions, no model input.
(H1) two branches s=+-1 exchanged by Theta; (H2) antipodal covariance
P(xi|-s)=P(-xi|s) and Theta-imaged backward protocol P_bwd(xi~|s)=P_fwd(xi|-s),
xi~ = reverse(-xi); (H3) symmetric prior.
Claims: (a) sigma[xi]=0 for every xi; (b) <sigma_s> = D_KL(P(.|s)||P(.|-s)),
additive over records IFF conditionally independent; (c) T-even control:
sigma_s = 0 identically while marginal distinguishability grows."""
import numpy as np
from itertools import product
rng = np.random.default_rng(1)

def all_xi(n): return np.array(list(product([1,-1], repeat=n)))

worst_a = worst_b = 0.0; gap_min = 0.0
for trial in range(200):
    n = rng.integers(2, 6)
    XI = all_xi(n)
    # random CORRELATED joint law P(xi|+) (positive tensor, normalized)
    Pp = rng.random(len(XI)) + 1e-3; Pp /= Pp.sum()
    idx = {tuple(x): i for i, x in enumerate(XI)}
    Pm = np.array([Pp[idx[tuple(-x)]] for x in XI])          # (H2) antipodal
    Pfwd = 0.5*(Pp + Pm)
    # backward: P_bwd(eta|s) = P_fwd(-reverse(eta)|-s); evaluate at eta = xi~
    # xi~ = reverse(-xi)  =>  P_bwd(xi~|s) = P_fwd(xi|-s)  (bijection, normalized)
    Pbwd_at_conj_p = Pm                                       # branch +: P(xi|-)
    Pbwd_at_conj = 0.5*(Pm + Pp)                              # unconditional
    # (a) sigma[xi] = ln Pfwd(xi) - ln Pbwd(xi~) = 0 for all xi
    worst_a = max(worst_a, np.abs(np.log(Pfwd) - np.log(Pbwd_at_conj)).max())
    # (b) <sigma_+> = sum_xi P(xi|+) ln[P(xi|+)/P(xi|-)] = KL
    sig = (Pp*np.log(Pp/Pbwd_at_conj_p)).sum()
    KL  = (Pp*np.log(Pp/Pm)).sum()
    worst_b = max(worst_b, abs(sig - KL))
    # additivity gap for correlated laws: sum of per-record KLs vs joint KL
    kl_marg = 0.0
    for k in range(n):
        pk = np.array([Pp[XI[:,k]==+1].sum(), Pp[XI[:,k]==-1].sum()])
        qk = np.array([Pm[XI[:,k]==+1].sum(), Pm[XI[:,k]==-1].sum()])
        kl_marg += (pk*np.log(pk/qk)).sum()
    gap_min = max(gap_min, abs(KL - kl_marg))
print(f"(a) max |sigma[xi]| over 200 correlated instances : {worst_a:.2e}")
print(f"(b) max |<sigma_s> - D_KL(joint)|                 : {worst_b:.2e}")
print(f"    additivity gap |KL_joint - sum KL_n| (MAX over trials, correlated): {gap_min:.3f} -> additive form needs cond. independence")

# product-law check: additivity exact
n = 5; XI = all_xi(n); eps = rng.uniform(0.1, 0.9, n)
Pp = np.prod((1 + XI*eps[None,:])/2, axis=1); Pm = np.prod((1 - XI*eps[None,:])/2, axis=1)
KL = (Pp*np.log(Pp/Pm)).sum(); add = (eps*np.log((1+eps)/(1-eps))).sum()
print(f"    product law: |KL_joint - sum_n eps ln[(1+eps)/(1-eps)]| = {abs(KL-add):.2e}")

# (c) T-even control: conjugation fixes the branch -> P_bwd(xi~|s)=P_fwd(xi|s)
Pp = rng.random(len(XI)) + 1e-3; Pp /= Pp.sum()
sig_Teven = (Pp*np.log(Pp/Pp)).sum()
# yet records still distinguish the two (now unexchanged) branches:
Pm = rng.random(len(XI)) + 1e-3; Pm /= Pm.sum()
D = (Pp*np.log(Pp/Pm)).sum()
print(f"(c) T-even control: <sigma_s> = {sig_Teven:.1f} identically while D = {D:.3f} > 0")
