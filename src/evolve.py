"""Brute-force validation: full state vector on dim 3*2^N (spec §8)."""
import numpy as np
from scipy.sparse import kron, identity, diags, csr_matrix
from scipy.sparse.linalg import expm_multiply
from model import h_system, current_op, chiral_states, I0

def build_hamiltonian(g):
    N = len(g)
    HS = csr_matrix(h_system())
    Iop = csr_matrix(current_op() / I0)
    dimE = 2 ** N
    # B = sum_n g_n sz^(n), diagonal in computational basis
    bits = ((np.arange(dimE)[:, None] >> np.arange(N)[::-1][None, :]) & 1)
    diagB = ((1 - 2 * bits) * g[None, :]).sum(axis=1)
    B = diags(diagB)
    return kron(HS, identity(dimE)) + kron(Iop, B)

def initial_state(N):
    kp, km = chiral_states()
    psiS = (kp + km) / np.sqrt(2)
    psiE = np.ones(2 ** N) / np.sqrt(2 ** N)  # |+x>^N
    return np.kron(psiS, psiE)

def evolve(H, psi0, ts):
    """Return states at the requested times (expm_multiply over intervals)."""
    out, psi, t_prev = [], psi0.copy(), 0.0
    for t in ts:
        if t > t_prev:
            psi = expm_multiply(-1j * (t - t_prev) * H, psi)
            t_prev = t
        out.append(psi.copy())
    return out

def rdm_SF(psi, N, frag):
    """Reduced density matrix of S + fragment (list of qubit indices)."""
    T = psi.reshape((3,) + (2,) * N)
    frag = list(frag)
    comp = [n for n in range(N) if n not in frag]
    perm = [0] + [1 + n for n in frag] + [1 + n for n in comp]
    M = np.transpose(T, perm).reshape(3 * 2 ** len(frag), 2 ** len(comp))
    return M @ M.conj().T

def rdm_S(psi, N):
    M = psi.reshape(3, 2 ** N)
    return M @ M.conj().T

def vn_entropy(rho):
    w = np.linalg.eigvalsh(rho)
    w = w[w > 1e-14]
    return float(-(w * np.log2(w)).sum())

def mutual_information_num(psi, N, frag):
    rho_SF = rdm_SF(psi, N, frag)
    rho_S = rdm_S(psi, N)
    # rho_F: trace S out of rho_SF
    m = len(frag)
    R = rho_SF.reshape(3, 2 ** m, 3, 2 ** m)
    rho_F = np.einsum('imin->mn', R)
    return vn_entropy(rho_S) + vn_entropy(rho_F) - vn_entropy(rho_SF)

def sample_polarization_num(psi, N, n_traj, rng):
    """Measure all bath qubits in sigma_y (exact joint distribution), then
    conditional <I>_xi. Rotate each qubit so sigma_y eigenbasis -> computational."""
    V = np.array([[1, 1], [1j, -1j]], complex) / np.sqrt(2)  # columns |+y>,|-y>
    T = psi.reshape((3,) + (2,) * N)
    for n in range(N):
        T = np.tensordot(T, V.conj().T, axes=([1 + n], [1]))
        T = np.moveaxis(T, -1, 1 + n)
    A = T.reshape(3, 2 ** N)                      # amplitudes site x bitstring
    probs = (np.abs(A) ** 2).sum(axis=0)
    probs = probs / probs.sum()
    picks = rng.choice(2 ** N, size=n_traj, p=probs)
    Iop = current_op()
    pol = np.empty(n_traj)
    for i, b in enumerate(picks):
        v = A[:, b]; v = v / np.linalg.norm(v)
        pol[i] = np.real(np.vdot(v, Iop @ v)) / I0
    return float(np.mean(np.abs(pol))), pol
