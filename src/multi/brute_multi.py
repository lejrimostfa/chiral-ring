"""Brute-force validation of model_multi.py (Paper 2, work-order step 3).

Full state vector, M=2 rings (dim 3) + N bath qubits: dim 9*2^N (N=8 -> 2304).
Pattern follows ../evolve.py; the single-ring mirror check reuses evolve.py
directly (Paper-1 code, not reimplemented).

Checks (machine precision unless marked statistical):
  V1  product init: reduced ring 4x4 == rho_joint (Gamma)                exact
  V2  per-branch bath sigma_y bias == sin(2 a_n(S) t)                     exact
  V3  brute-force record sampling P_i == closed-form sample_joint         statistical
  V4  entangled init (|+-|+|-+>)/sqrt2: concurrence == |Gamma_(+-,-+)|    exact
  V5  mirror prep (single ring): reduced coherence == prod cos(2g(T-t))   exact

Run:  python3 brute_multi.py    (from src/multi)
"""
import os, sys
import numpy as np
from scipy.sparse import kron, identity, diags, csr_matrix

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from model import h_system, current_op, chiral_states, I0
import evolve as evolve1   # Paper-1 single-ring brute force (reused for mirror)
import model_multi as mm


# ----------------------------------------------------------------------
# M=2 ring pair + bath
# ----------------------------------------------------------------------
def build_hamiltonian(g_mat):
    """H = HS1 + HS2 + (I1/I0)⊗B1 + (I2/I0)⊗B2, dim 9*2^N."""
    M, N = g_mat.shape
    assert M == 2
    dimE = 2 ** N
    HS = csr_matrix(h_system())
    Iop = csr_matrix(current_op() / I0)
    I3 = identity(3)
    IE = identity(dimE)
    bits = ((np.arange(dimE)[:, None] >> np.arange(N)[::-1][None, :]) & 1)
    diagB1 = ((1 - 2 * bits) * g_mat[0][None, :]).sum(axis=1)
    diagB2 = ((1 - 2 * bits) * g_mat[1][None, :]).sum(axis=1)
    B1, B2 = diags(diagB1), diags(diagB2)
    return (kron(HS, kron(I3, IE)) + kron(I3, kron(HS, IE))
            + kron(Iop, kron(I3, B1)) + kron(I3, kron(Iop, B2)))


def initial_state(N, kind="product"):
    kp, km = chiral_states()
    sup = (kp + km) / np.sqrt(2)
    if kind == "product":
        psiS = np.kron(sup, sup)
    elif kind == "entangled":                      # (|+-> + |-+>)/sqrt2
        psiS = (np.kron(kp, km) + np.kron(km, kp)) / np.sqrt(2)
    else:
        raise ValueError(kind)
    psiE = np.ones(2 ** N) / np.sqrt(2 ** N)
    return np.kron(psiS, psiE)


def evolve(H, psi0, ts):
    from scipy.sparse.linalg import expm_multiply
    out, psi, t_prev = [], psi0.copy(), 0.0
    for t in ts:
        if t > t_prev:
            psi = expm_multiply(-1j * (t - t_prev) * H, psi)
            t_prev = t
        out.append(psi.copy())
    return out


def chiral_basis4():
    """The four |s1 s2> vectors (9-dim) in mm.all_branches(2) order:
    (+,+),(+,-),(-,+),(-,-)."""
    kp, km = chiral_states()
    ring = {1: kp, -1: km}
    S_all = mm.all_branches(2)
    return np.array([np.kron(ring[s1], ring[s2]) for s1, s2 in S_all])  # (4,9)


def reduced_chiral(psi, N):
    """4x4 reduced ring density matrix in the chiral basis (== rho_joint)."""
    M9 = psi.reshape(9, 2 ** N)
    rho9 = M9 @ M9.conj().T
    Bmat = chiral_basis4()               # (4,9)
    return Bmat.conj() @ rho9 @ Bmat.T   # <a|rho|b>


def conditional_bath(psi, N, S):
    """Bath state conditional on the rings being in branch S (project, normalise)."""
    kp, km = chiral_states()
    ring = {1: kp, -1: km}
    braS = np.kron(ring[S[0]], ring[S[1]]).conj()   # <S|
    M9 = psi.reshape(9, 2 ** N)
    v = braS @ M9
    nv = np.linalg.norm(v)
    return v / nv if nv > 1e-14 else v


def sy_expectations(vec, N):
    """<sigma_y> on each of N qubits for a (possibly unnormalised) 2^N vector."""
    T = vec.reshape((2,) * N)
    SY = np.array([[0, -1j], [1j, 0]], complex)
    out = np.empty(N)
    for n in range(N):
        Tn = np.moveaxis(T, n, 0).reshape(2, -1)
        rho1 = Tn @ Tn.conj().T
        tr = np.real(np.trace(rho1))
        out[n] = np.real(np.trace(SY @ rho1)) / tr if tr > 1e-14 else 0.0
    return out


def concurrence_wootters(rho4):
    """Wootters concurrence of a 2-qubit (4x4) density matrix."""
    sy = np.array([[0, -1j], [1j, 0]], complex)
    Y = np.kron(sy, sy)
    R = rho4 @ Y @ rho4.conj() @ Y
    ev = np.sqrt(np.clip(np.real(np.linalg.eigvals(R)), 0, None))
    ev = np.sort(ev)[::-1]
    return max(0.0, ev[0] - ev[1] - ev[2] - ev[3])


# ----------------------------------------------------------------------
# Checks
# ----------------------------------------------------------------------
def v1_coherence(g_mat, ts):
    H = build_hamiltonian(g_mat)
    states = evolve(H, initial_state(len(g_mat[0]), "product"), ts)
    N = g_mat.shape[1]
    err = 0.0
    for t, psi in zip(ts, states):
        rho_bf = reduced_chiral(psi, N)
        rho_cf = mm.rho_joint(g_mat, t)
        err = max(err, np.max(np.abs(rho_bf - rho_cf)))
    ok = err < 1e-10
    print(f"  V1  reduced ring 4x4 == rho_joint(Gamma)   max|err|={err:.2e}  "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def v2_bias(g_mat, ts):
    H = build_hamiltonian(g_mat)
    states = evolve(H, initial_state(len(g_mat[0]), "product"), ts)
    N = g_mat.shape[1]
    err = 0.0
    for t, psi in zip(ts, states):
        for S in mm.all_branches(2):
            bath = conditional_bath(psi, N, S)
            sy = sy_expectations(bath, N)
            sy_cf = mm.record_bias(g_mat, S, t)
            err = max(err, np.max(np.abs(sy - sy_cf)))
    ok = err < 1e-9
    print(f"  V2  per-branch bath <sigma_y> == sin(2 a_n(S)t) max|err|={err:.2e}  "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def v3_records(g_mat, t, n_traj=60_000):
    """Sample sigma_y records from the EXACT state; posterior P_i vs closed form."""
    H = build_hamiltonian(g_mat)
    N = g_mat.shape[1]
    psi = evolve(H, initial_state(N, "product"), [t])[0]
    # rotate each qubit sigma_y eigenbasis -> computational (as ../evolve.py)
    V = np.array([[1, 1], [1j, -1j]], complex) / np.sqrt(2)
    T = psi.reshape((3, 3) + (2,) * N)
    for n in range(N):
        T = np.tensordot(T, V.conj().T, axes=([2 + n], [1]))
        T = np.moveaxis(T, -1, 2 + n)
    A = T.reshape(9, 2 ** N)
    probs = (np.abs(A) ** 2).sum(axis=0)
    probs /= probs.sum()
    rng = np.random.default_rng(12345)
    picks = rng.choice(2 ** N, size=n_traj, p=probs)
    bit = ((picks[:, None] >> np.arange(N)[::-1][None, :]) & 1)
    r = (1 - 2 * bit).astype(float)                 # sigma_y outcome +1/-1
    # closed-form posterior over 4 branches
    S_all = mm.all_branches(2)
    mu = np.clip(np.sin(2.0 * (S_all @ g_mat) * t), -1 + 1e-12, 1 - 1e-12)
    logL = np.log1p(r[:, None, :] * mu[None, :, :]).sum(axis=2)
    logL -= logL.max(axis=1, keepdims=True)
    post = np.exp(logL); post /= post.sum(axis=1, keepdims=True)
    m_i = post @ S_all.astype(float)
    P_i_bf = np.mean(np.abs(m_i), axis=0)
    sem = np.std(np.abs(m_i), axis=0) / np.sqrt(n_traj)
    cf = mm.sample_joint(g_mat, t, n_traj=n_traj, rng=np.random.default_rng(999))
    P_i_cf = cf["P_i"]
    ok = np.all(np.abs(P_i_bf - P_i_cf) < 5 * (sem + 1 / np.sqrt(n_traj)))
    print(f"  V3  brute-force P_i={np.round(P_i_bf,4)} vs closed {np.round(P_i_cf,4)}"
          f"  {'PASS' if ok else 'FAIL'}")
    return ok


def v4_entangled(g_mat, ts):
    H = build_hamiltonian(g_mat)
    states = evolve(H, initial_state(len(g_mat[0]), "entangled"), ts)
    N = g_mat.shape[1]
    err = 0.0
    for t, psi in zip(ts, states):
        rho4 = reduced_chiral(psi, N)
        C_bf = concurrence_wootters(rho4)
        C_cf = mm.concurrence_pm(g_mat, t)
        err = max(err, abs(C_bf - C_cf))
    ok = err < 1e-9
    print(f"  V4  entangled concurrence == |Gamma_(+-,-+)|  max|err|={err:.2e}  "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def v5_mirror(g, T, ts):
    """Single-ring mirror: reduced coherence == prod cos(2 g (T-t)); records
    consumed (distinguishability decreasing). Reuses evolve.py (Paper 1)."""
    N = len(g)
    H = evolve1.build_hamiltonian(g)
    psi0 = evolve1.initial_state(N)
    psiT = evolve1.evolve(H, psi0, [T])[0]
    # Theta = I_ring (x) (i sigma_y)^N  o  K   (baseline pure-dephasing conjugation)
    psi_conj = np.conj(psiT).reshape((3,) + (2,) * N)
    iSY = np.array([[0, 1], [-1, 0]], complex)      # i*sigma_y
    Tt = psi_conj
    for n in range(N):
        Tt = np.tensordot(Tt, iSY.T, axes=([1 + n], [1]))
        Tt = np.moveaxis(Tt, -1, 1 + n)
    psi_mir0 = Tt.reshape(3 * 2 ** N)
    states = evolve1.evolve(H, psi_mir0, ts)
    kp, km = chiral_states()
    err, Ds = 0.0, []
    for t, psi in zip(ts, states):
        M3 = psi.reshape(3, 2 ** N)
        rho3 = M3 @ M3.conj().T
        coh = 2.0 * np.real(kp.conj() @ rho3 @ km)  # 2*rho_{+-}
        coh_cf = float(np.prod(np.cos(2 * g * (T - t))))
        err = max(err, abs(coh - coh_cf))
        eps = np.clip(np.abs(np.sin(2 * g * (T - t))), 0, 1 - 1e-12)
        Ds.append(float(np.sum(eps * np.log((1 + eps) / (1 - eps)))))
    decreasing = all(Ds[i + 1] <= Ds[i] + 1e-12 for i in range(len(Ds) - 1))
    ok = err < 1e-9 and decreasing
    print(f"  V5  mirror coherence == prod cos(2g(T-t))  max|err|={err:.2e}  "
          f"D decreasing={decreasing}  {'PASS' if ok else 'FAIL'}")
    return ok


def run_all():
    print("=== brute_multi.py validation (M=2, N=8) ===")
    rng = np.random.default_rng(2024)
    N = 8
    g_mat = rng.uniform(0.1, 0.5, size=(2, N))       # generic shared bath
    ts = [0.2, 0.5, 0.9, 1.4]
    results = [
        v1_coherence(g_mat, ts),
        v2_bias(g_mat, ts),
        v3_records(g_mat, 0.6),
        v4_entangled(g_mat, ts),
        # mirror in the monotonic-consumption regime (2*g_max*T < pi/2), so the
        # "records consumed" claim is a genuine monotone decrease of D.
        v5_mirror(rng.uniform(0.05, 0.35, size=N), T=1.5,
                  ts=np.linspace(0, 1.5, 12)),
    ]
    allok = all(results)
    print(f"ALL: {'PASS' if allok else 'FAIL'}")
    return allok


if __name__ == "__main__":
    sys.exit(0 if run_all() else 1)
