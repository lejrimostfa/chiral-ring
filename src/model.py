"""Chiral ring model: system operators and unit tests (spec §1)."""
import numpy as np

J = 1.0
I0 = np.sqrt(3.0) * J

def h_system():
    """3-site ring hopping Hamiltonian."""
    H = np.zeros((3, 3), complex)
    for j in range(3):
        H[j, (j + 1) % 3] = -J
        H[(j + 1) % 3, j] = -J
    return H

def current_op():
    """Current (chirality) operator, T-odd."""
    I = np.zeros((3, 3), complex)
    for j in range(3):
        I[(j + 1) % 3, j] = 1j * J
        I[j, (j + 1) % 3] = -1j * J
    return I

def chiral_states():
    """|+> and |-> : quasimomentum k = +-2pi/3 states."""
    ks = [2 * np.pi / 3, -2 * np.pi / 3]
    vecs = []
    for k in ks:
        v = np.exp(1j * k * np.arange(3)) / np.sqrt(3)
        vecs.append(v)
    return vecs[0], vecs[1]

def run_unit_tests():
    HS, I = h_system(), current_op()
    kp, km = chiral_states()
    # [H_S, I] = 0
    assert np.allclose(HS @ I - I @ HS, 0), "[H_S, I] != 0"
    # eigenvalues of I on chiral states: +-I0 ; k=0 state: 0
    assert np.allclose(I @ kp, I0 * kp), "I|+> != +I0|+>"
    assert np.allclose(I @ km, -I0 * km), "I|-> != -I0|->"
    k0 = np.ones(3) / np.sqrt(3)
    assert np.allclose(I @ k0, 0), "I|k=0> != 0"
    # degeneracy of the doublet
    assert np.allclose(HS @ kp, J * kp) and np.allclose(HS @ km, J * km)
    # T = complex conjugation in site basis: T|+> = |->, T I T^-1 = -I
    assert np.allclose(np.conj(kp), km), "T|+> != |->"
    assert np.allclose(np.conj(I), -I), "I not T-odd under K_ring"
    # initial state: T-invariant superposition, <I> = 0
    psi0 = (kp + km) / np.sqrt(2)
    assert abs(np.vdot(psi0, I @ psi0)) < 1e-14, "<I>_0 != 0"
    ph = np.vdot(np.conj(psi0), psi0)
    assert abs(abs(ph) - 1) < 1e-12, "psi0 not T-invariant up to phase"
    print("model.py unit tests: ALL PASSED")

if __name__ == "__main__":
    run_unit_tests()
