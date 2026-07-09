"""Paper 3 — verification gates for the lattice engine (sprint P3-A).

verify_note.py / verify_multi.py pattern. G1-G4 are the prototype gates ported
onto the production module; G5-G6 are the new collision-protocol gates. NO
figure is produced anywhere in Paper 3 until all six print PASS.

  G1  Delta = I + C_bwd on the lattice-generated law (Paper-2 theorem)   1e-12
  G2  f=0 -> Paper-1 Theorem 2 exact (<sigma_i> = D_i)                   1e-12
  G3  L=3 brute-force state vector vs branch construction                1e-10
  G4  local marginal == global marginal (sign + correlation)             sampling
  G5  COLLISION: exact additivity (deficit = 0) + <sigma_i> = D_i,
      disjoint baths                                                     1e-12
  G6  COLLISION L=2: reproduces the Paper-2 C5 f=0 limit,
      <sigma_forward> = +D, <sigma_mirror> = -D                          1e-12

Run:  python3 verify_lattice.py    (from src/lattice)
"""
import os, sys
from itertools import product
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lattice import (build_chain, all_branches, beta_matrix, collision_beta,
                     draw_signs, sample_records, productions,
                     local_marginal_sigma, local_D, _flip_index)


# ---------------------------------------------------------------- exhaustive
def _exhaustive_prod(beta, S_all, draw_sign):
    """Exact <sigma_S> and Option-A marginals <sigma_i> over ALL 2^N records.
    Records generated from the draw law (fresh +1 / charged -1), scored with the
    forward likelihood. Returns (mean sigma_S, mean sigma_i per ring)."""
    B, N = beta.shape
    xs = np.array(list(product([1, -1], repeat=N)))               # (X, N)
    Pf = np.clip(np.prod((1 + xs[:, None, :] * beta[None, :, :]) / 2, axis=2),
                 1e-300, None)                                    # (X, B) forward
    Pdraw = np.prod((1 + draw_sign[None, None, :] * xs[:, None, :]
                     * beta[None, :, :]) / 2, axis=2)             # (X, B) generating
    flip = _flip_index(S_all)
    L = S_all.shape[1]
    Pm = {(i, s): Pf[:, S_all[:, i] == s].mean(axis=1)
          for i in range(L) for s in (1, -1)}
    sigS = np.empty(B); sigi = np.empty((B, L))
    for b in range(B):
        w = Pdraw[:, b]
        sigS[b] = float((w * (np.log(Pf[:, b]) - np.log(Pf[:, flip[b]]))).sum())
        for i in range(L):
            si = S_all[b, i]
            sigi[b, i] = float((w * (np.log(Pm[(i, si)])
                                     - np.log(Pm[(i, -si)]))).sum())
    return sigS.mean(), sigi.mean(axis=0)


# ---------------------------------------------------------------- G1
def gate_G1_decomposition(tol=1e-12):
    G, Tc, owner = build_chain(L=2, Np=1, Ns=2, seed=2)
    t = 0.7; S_all = all_branches(2); beta = beta_matrix(G, Tc, t, S_all)
    N = beta.shape[1]
    xs = np.array(list(product([1, -1], repeat=N)))
    P = np.prod((1 + xs[:, None, :] * beta[None, :, :]) / 2, axis=2)
    Pxi = P.mean(axis=1)
    P1 = {s: P[:, S_all[:, 0] == s].mean(axis=1) for s in (1, -1)}
    P2 = {s: P[:, S_all[:, 1] == s].mean(axis=1) for s in (1, -1)}
    flip = _flip_index(S_all)
    Delta = I = Cb = 0.0
    for b in range(4):
        s1, s2 = S_all[b]
        w = P[:, b] / 4
        sS = np.log(P[:, b] / P[:, flip[b]])
        si = np.log(P1[s1] / P1[-s1]) + np.log(P2[s2] / P2[-s2])
        J = np.log(P[:, b] * Pxi / (P1[s1] * P2[s2]))
        Delta += (w * (sS - si)).sum()
        I += (w * J).sum()
        Cb += -(P[:, flip[b]] / 4 * J).sum()
    err = abs(Delta - (I + Cb))
    ok = err < tol and Delta >= -tol
    print(f"[{'PASS' if ok else 'FAIL'}] G1 Delta=I+C_bwd on lattice law: "
          f"err={err:.2e}, Delta={Delta:.4f}, I={I:.4f}, C_bwd={Cb:.4f}")
    return ok


# ---------------------------------------------------------------- G2
def gate_G2_paper1(tol=1e-12):
    G, Tc, owner = build_chain(L=2, Np=4, Ns=0, seed=1)
    t = 0.8; S_all = all_branches(2)
    worst = 0.0
    for i in (0, 1):
        Wi = [n for n, o in enumerate(owner) if o == ('priv', i)]
        eps = np.abs(np.sin(2 * G[i, Wi] * t))
        D = np.sum(eps * np.log((1 + eps) / (1 - eps)))
        E = 0.0
        for r in product([1, -1], repeat=len(Wi)):
            r = np.array(r)
            p = np.prod((1 + r * eps) / 2)
            E += p * np.sum(np.log((1 + r * eps) / (1 - r * eps)))
        worst = max(worst, abs(E - D))
    ok = worst < tol
    print(f"[{'PASS' if ok else 'FAIL'}] G2 f=0 -> Paper-1 Th2 exact: {worst:.2e}")
    return ok


# ---------------------------------------------------------------- G3
def gate_G3_bruteforce(tol=1e-10):
    from scipy.linalg import expm
    G, Tc, owner = build_chain(L=3, Np=1, Ns=1, seed=3)
    L, N = G.shape; t = 0.6
    sz = np.diag([1., -1.]).astype(complex); I2 = np.eye(2, dtype=complex)

    def kron(ops):
        out = np.array([[1.]], dtype=complex)
        for o in ops:
            out = np.kron(out, o)
        return out
    Zi = [kron([sz if k == i else I2 for k in range(L)]) for i in range(L)]
    H = np.zeros((2 ** L * 2 ** N,) * 2, dtype=complex)
    for n in range(N):
        zn = kron([sz if k == n else I2 for k in range(N)])
        A = sum(G[i, n] * Zi[i] for i in range(L))
        H += np.kron(A, zn)
    plus = np.ones(2) / np.sqrt(2)
    psi0 = kron([plus.reshape(2, 1)] * L).ravel()
    psi0 = np.kron(psi0, kron([plus.reshape(2, 1)] * N).ravel())
    psi_t = expm(-1j * H * t) @ psi0
    S_all = all_branches(L)
    psi_ana = np.zeros_like(psi_t)
    for b, S in enumerate(S_all):
        a = S @ G
        bath = kron([(expm(-1j * a[n] * sz * t) @ plus).reshape(2, 1)
                     for n in range(N)]).ravel()
        ring = np.zeros(2 ** L, dtype=complex)
        idx = int(''.join('0' if s == 1 else '1' for s in S), 2)
        ring[idx] = 2 ** (-L / 2)
        psi_ana += np.kron(ring, bath)
    err = np.abs(psi_t - psi_ana).max()
    ok = err < tol
    print(f"[{'PASS' if ok else 'FAIL'}] G3 L=3 brute force vs branches: {err:.2e}")
    return ok


# ---------------------------------------------------------------- G4
def gate_G4_local_vs_global(L=4, n_traj=4000, seed=5):
    G, Tc, owner = build_chain(L=L, Np=3, Ns=2, seed=seed)
    t = 0.8; S_all = all_branches(L); beta = beta_matrix(G, Tc, t, S_all)
    rng = np.random.default_rng(seed); true_idx = 0
    R = sample_records(beta[true_idx], n_traj, rng)
    _, sig_g = productions(R, beta, S_all, true_idx)
    S_true = S_all[true_idx]
    agree, corr = [], []
    for i in range(L):
        sl = local_marginal_sigma(R, G, Tc, t, S_true, i)
        agree.append(np.sign(sig_g[:, i].mean()) == np.sign(sl.mean()))
        corr.append(np.corrcoef(sig_g[:, i], sl)[0, 1])
    ok = all(agree)
    print(f"[{'PASS' if ok else 'FAIL'}] G4 local vs global marginal: "
          f"sign-agree={ok}, mean corr={np.mean(corr):.3f} "
          f"{[f'{c:.3f}' for c in corr]}")
    return ok


# ---------------------------------------------------------------- G5
def gate_G5_collision_additivity(dt=0.3, tol=1e-12):
    """Collision protocol, disjoint baths: exact additivity (deficit = 0) and
    <sigma_i> = D_i (closed form). All fresh witnesses, fixed window dt."""
    G, Tc, owner = build_chain(L=2, Np=4, Ns=0, seed=1)          # disjoint
    S_all = all_branches(2)
    beta = collision_beta(G, Tc, dt, S_all)
    sigS, sigi = _exhaustive_prod(beta, S_all, draw_signs(Tc))
    deficit = sigS - sigi.sum()
    D = np.array([local_D(G, Tc, dt, i) for i in range(2)])
    err_add = abs(deficit)
    err_D = np.max(np.abs(sigi - D))
    ok = err_add < tol and err_D < tol
    print(f"[{'PASS' if ok else 'FAIL'}] G5 collision additivity: "
          f"deficit={err_add:.2e}, |<sigma_i>-D_i|={err_D:.2e}")
    return ok


# ---------------------------------------------------------------- G6
def gate_G6_collision_paper2(dt=0.3, T=1.0, tol=1e-12):
    """Collision L=2, ring 0 fresh + ring 1 charged (contrarian), disjoint:
    reproduces the Paper-2 C5 f=0 limit, <sigma_0>=+D_0, <sigma_1>=-D_1."""
    G, Tc, owner = build_chain(L=2, Np=4, Ns=0, seed=1,
                               contrarians=(1,), T_charge=T)
    S_all = all_branches(2)
    beta = collision_beta(G, Tc, dt, S_all)
    sigS, sigi = _exhaustive_prod(beta, S_all, draw_signs(Tc))
    D = np.array([local_D(G, Tc, dt, i) for i in range(2)])
    expected = np.array([+D[0], -D[1]])                          # fwd +D, mirror -D
    err = np.max(np.abs(sigi - expected))
    ok = err < tol
    print(f"[{'PASS' if ok else 'FAIL'}] G6 collision L=2 == Paper-2 C5 f=0: "
          f"<sigma>={np.round(sigi,4)} vs [+D0,-D1]={np.round(expected,4)} "
          f"err={err:.2e}")
    return ok


if __name__ == "__main__":
    print("=== PAPER 3 GATES (lattice engine) ===")
    results = [gate_G1_decomposition(), gate_G2_paper1(), gate_G3_bruteforce(),
               gate_G4_local_vs_global(), gate_G5_collision_additivity(),
               gate_G6_collision_paper2()]
    allok = all(results)
    print(f"\n{sum(results)}/6 gates PASS -> "
          f"{'CLEARED for figures' if allok else 'BLOCKED (no figures)'}")
    sys.exit(0 if allok else 1)
