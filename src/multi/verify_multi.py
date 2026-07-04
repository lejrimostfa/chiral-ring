"""Independent verification suite — deficit decomposition theorem (task B.2).

Pattern follows ../verify_note.py. RECODED FROM SCRATCH: no import of
model_multi / analyze_c2 / run_scenarios; every quantity below is built
directly from the definitions in notes/multi_definitions.md §6ter, by explicit
enumeration (M=2).

Parameterization: for M=2 any flip-symmetric product record law has per-record
bias EXACTLY linear in the branch variables (§6ter lemma; the odd functions on
{±1}² are spanned by s1 and s2):
    x_n(S) = s1*u_n + s2*v_n,   validity |u_n|+|v_n| <= 1.
Physical collision biases: u_n = sin(2 g1n t) cos(2 g2n t),
                            v_n = cos(2 g1n t) sin(2 g2n t).

CHECK 1: Delta = TC + C_bwd (Eq. 22-24) to machine precision, on three
         configurations (disjoint, shared f=0.5, DFS identical) x three times;
         sub-checks: disjoint gives Delta = TC = C_bwd = 0 identically.
CHECK 2: small-t limit. (a) SYMBOLIC (sympy), N=1 and N=2, generic couplings:
         leading order Delta = 4 (sum_n u_n v_n)^2, TC = (1/2)(sum_n u_n v_n)^2,
         hence TC/Delta -> 1/8 exactly at leading order (lambda^4).
         (b) NUMERIC confirmation of the closed forms at N=10, t=0.02.
CHECK 3: sign survey of C_bwd: exhaustive N=1 grid + random product laws
         N in {2,3,4}. Reports the minimum found and any counterexample.

Run:  python3 verify_multi.py    (from src/multi)
"""
import itertools
import sys
import numpy as np

LN = np.log


# ----------------------------------------------------------------------
# From-scratch exact machinery (M=2)
# ----------------------------------------------------------------------
BRANCHES = [(1, 1), (1, -1), (-1, 1), (-1, -1)]


def record_probs(U, V):
    """P(xi|S) table, shape (4, 2^N), by explicit product. Biases x_n = s1 u + s2 v."""
    N = len(U)
    xis = list(itertools.product((1, -1), repeat=N))
    P = np.empty((4, 2 ** N))
    for b, (s1, s2) in enumerate(BRANCHES):
        for j, xi in enumerate(xis):
            p = 1.0
            for r, u, v in zip(xi, U, V):
                p *= (1.0 + r * (s1 * u + s2 * v)) / 2.0
            P[b, j] = p
    return P


def deficit_terms(U, V):
    """(Delta, TC, C_bwd) exactly, straight from the §6ter definitions."""
    P = np.clip(record_probs(U, V), 1e-300, None)     # (4, 2^N)
    anti = [3, 2, 1, 0]                               # -S for BRANCHES order
    Pm = P.mean(axis=0)                               # P(xi)
    P_s1 = {+1: P[[0, 1]].mean(axis=0), -1: P[[2, 3]].mean(axis=0)}
    P_s2 = {+1: P[[0, 2]].mean(axis=0), -1: P[[1, 3]].mean(axis=0)}

    # Delta = <sigma_S> - <sigma_1> - <sigma_2>
    D_joint = np.mean([np.sum(P[b] * (LN(P[b]) - LN(P[anti[b]])))
                       for b in range(4)])
    D1 = 0.5 * sum(np.sum(P_s1[s] * (LN(P_s1[s]) - LN(P_s1[-s])))
                   for s in (+1, -1))
    D2 = 0.5 * sum(np.sum(P_s2[s] * (LN(P_s2[s]) - LN(P_s2[-s])))
                   for s in (+1, -1))
    Delta = D_joint - D1 - D2

    # TC via posterior entropies (independent route, no J shortcut)
    post = (P / 4.0) / Pm[None, :]
    post = np.clip(post, 1e-300, 1.0)
    H_S = -np.sum(post * LN(post), axis=0)
    p1 = np.clip(post[0] + post[1], 1e-300, 1 - 1e-15)   # P(s1=+1|xi)
    p2 = np.clip(post[0] + post[2], 1e-300, 1 - 1e-15)
    H1 = -(p1 * LN(p1) + (1 - p1) * LN(1 - p1))
    H2 = -(p2 * LN(p2) + (1 - p2) * LN(1 - p2))
    TC = float(np.sum(Pm * (H1 + H2 - H_S)))

    # C_bwd from the exact expression (24)
    J = np.empty_like(P)
    for b, (s1, s2) in enumerate(BRANCHES):
        J[b] = LN(P[b]) + LN(Pm) - LN(P_s1[s1]) - LN(P_s2[s2])
    C_bwd = -np.mean([np.sum(P[anti[b]] * J[b]) for b in range(4)])
    return float(Delta), TC, float(C_bwd)


def uv_from_couplings(g1, g2, t):
    return np.sin(2 * g1 * t) * np.cos(2 * g2 * t), \
        np.cos(2 * g1 * t) * np.sin(2 * g2 * t)


# ----------------------------------------------------------------------
# CHECK 1: the identity on three configurations
# ----------------------------------------------------------------------
def check1():
    print("=" * 60)
    print("CHECK 1: Delta = TC + C_bwd (three configs, machine precision)")
    rng = np.random.default_rng(21)
    N = 8
    ok = True
    configs = {}
    # disjoint: first half ring-1 only, second half ring-2 only
    g1 = np.zeros(N); g2 = np.zeros(N)
    g1[:N // 2] = rng.uniform(0.1, 0.4, N // 2)
    g2[N // 2:] = rng.uniform(0.1, 0.4, N // 2)
    configs["disjoint"] = (g1, g2)
    configs["shared f=0.5"] = (
        np.where(np.arange(N) < 3 * N // 4, rng.uniform(0.1, 0.4, N), 0.0),
        np.where((np.arange(N) < N // 2) | (np.arange(N) >= 3 * N // 4),
                 rng.uniform(0.1, 0.4, N), 0.0))
    gid = rng.uniform(0.1, 0.4, N)
    configs["DFS identical"] = (gid, gid.copy())
    for name, (a, b) in configs.items():
        for t in (0.3, 0.7, 1.2):
            U, V = uv_from_couplings(a, b, t)
            Delta, TC, C = deficit_terms(U, V)
            err = abs(Delta - (TC + C))
            ok &= err < 1e-10
            if name == "disjoint":
                ok &= abs(Delta) < 1e-10 and abs(TC) < 1e-10 and abs(C) < 1e-10
        print(f"  {name:14s}: Delta={Delta:9.5f} TC={TC:8.5f} "
              f"C_bwd={C:9.5f}  |identity err|={err:.2e}")
    print(f"  disjoint sub-check: Delta = TC = C_bwd = 0 identically")
    print(f"  ... PASS: {ok}")
    return ok


# ----------------------------------------------------------------------
# CHECK 2: small-t limit — symbolic (sympy) + numeric confirmation
# ----------------------------------------------------------------------
def check2():
    print("=" * 60)
    print("CHECK 2: small-t leading order and the 1/8 constant")
    import sympy as sp

    def symbolic_leading(Nsym):
        lam = sp.symbols("lam", positive=True)
        us = sp.symbols(f"u0:{Nsym}", real=True)
        vs = sp.symbols(f"v0:{Nsym}", real=True)
        xis = list(itertools.product((1, -1), repeat=Nsym))
        P = {}
        for S in BRANCHES:
            for xi in xis:
                p = sp.Rational(1)
                for r, u, v in zip(xi, us, vs):
                    p *= (1 + r * lam * (S[0] * u + S[1] * v)) / 2
                P[(S, xi)] = sp.expand(p)
        Pm = {xi: sum(P[(S, xi)] for S in BRANCHES) / 4 for xi in xis}
        Ps1 = {(s, xi): sum(P[(S, xi)] for S in BRANCHES if S[0] == s) / 2
               for s in (1, -1) for xi in xis}
        Ps2 = {(s, xi): sum(P[(S, xi)] for S in BRANCHES if S[1] == s) / 2
               for s in (1, -1) for xi in xis}

        order = 6

        def ser(e):
            return sp.series(e, lam, 0, order).removeO()

        # Delta = E_fwd[ ln P(xi|S) - ln P(xi|-S)
        #               - ln Ps1(s1)/Ps1(-s1) - ln Ps2(s2)/Ps2(-s2) ]
        Delta = 0
        for S in BRANCHES:
            mS = (-S[0], -S[1])
            for xi in xis:
                integrand = (sp.log(P[(S, xi)]) - sp.log(P[(mS, xi)])
                             - sp.log(Ps1[(S[0], xi)]) + sp.log(Ps1[(-S[0], xi)])
                             - sp.log(Ps2[(S[1], xi)]) + sp.log(Ps2[(-S[1], xi)]))
                Delta += sp.Rational(1, 4) * P[(S, xi)] * ser(integrand)
        Delta4 = sp.simplify(sp.expand(Delta).coeff(lam, 4))

        # TC = E_fwd[J],  J = ln P(xi|S) + ln P(xi) - ln Ps1 - ln Ps2
        TC = 0
        for S in BRANCHES:
            for xi in xis:
                J = (sp.log(P[(S, xi)]) + sp.log(Pm[xi])
                     - sp.log(Ps1[(S[0], xi)]) - sp.log(Ps2[(S[1], xi)]))
                TC += sp.Rational(1, 4) * P[(S, xi)] * ser(J)
        TC4 = sp.simplify(sp.expand(TC).coeff(lam, 4))

        w = sum(u * v for u, v in zip(us, vs))     # the collective invariant
        okD = sp.simplify(Delta4 - 4 * w ** 2) == 0
        okT = sp.simplify(TC4 - sp.Rational(1, 2) * w ** 2) == 0
        # lambda^2 coefficients must vanish (deficit starts at 4th order)
        okD2 = sp.simplify(sp.expand(Delta).coeff(lam, 2)) == 0
        okT2 = sp.simplify(sp.expand(TC).coeff(lam, 2)) == 0
        return okD and okD2, okT and okT2

    okD1, okT1 = symbolic_leading(1)
    print(f"  N=1 symbolic: Delta_lead = 4(u v)^2 : {okD1};  "
          f"TC_lead = (1/2)(u v)^2 : {okT1}")
    okD2, okT2 = symbolic_leading(2)
    print(f"  N=2 symbolic: Delta_lead = 4(u1 v1 + u2 v2)^2 : {okD2};  "
          f"TC_lead = (1/2)(u1 v1 + u2 v2)^2 : {okT2}")
    print("  => TC/Delta -> 1/8 EXACTLY at leading order (lambda^4), generic")
    print("     couplings (no equal-coupling assumption), collective invariant")
    print("     w = sum_n u_n v_n; both sides quadratic in w (N^2, not N).")

    # numeric confirmation at N=10, t=0.02
    rng = np.random.default_rng(11)
    g1 = rng.uniform(0.1, 0.4, 10); g2 = rng.uniform(0.1, 0.4, 10)
    t = 0.02
    U, V = uv_from_couplings(g1, g2, t)
    Delta, TC, C = deficit_terms(U, V)
    w = float(np.sum(U * V))
    relD = abs(Delta - 4 * w ** 2) / (4 * w ** 2)
    relT = abs(TC - w ** 2 / 2) / (w ** 2 / 2)
    okn = relD < 5e-3 and relT < 5e-3 and abs(TC / Delta - 0.125) < 1e-3
    print(f"  N=10 numeric t={t}: Delta vs 4w^2 rel.err={relD:.2e}; "
          f"TC vs w^2/2 rel.err={relT:.2e}; TC/Delta={TC/Delta:.5f}")
    ok = okD1 and okT1 and okD2 and okT2 and okn
    print(f"  ... PASS: {ok}")
    return ok


# ----------------------------------------------------------------------
# CHECK 3: sign survey of C_bwd
# ----------------------------------------------------------------------
def check3():
    print("=" * 60)
    print("CHECK 3: sign survey of C_bwd (violable?)")
    # exhaustive N=1: biases mu(++) = u+v, mu(+-) = u-v range over (-1,1)^2
    grid = np.linspace(-0.98, 0.98, 99)
    worst = np.inf; arg = None
    for mpp in grid:
        for mpm in grid:
            u, v = (mpp + mpm) / 2, (mpp - mpm) / 2
            _, _, C = deficit_terms(np.array([u]), np.array([v]))
            if C < worst:
                worst, arg = C, (mpp, mpm)
    print(f"  N=1 exhaustive grid (99x99): min C_bwd = {worst:.3e} "
          f"at (mu_pp, mu_pm)={np.round(arg, 3)}")
    ok1 = worst > -1e-12

    # random product laws, N in {2,3,4}
    rng = np.random.default_rng(5)
    worstN = np.inf; argN = None
    for N in (2, 3, 4):
        for _ in range(2000):
            mpp = rng.uniform(-0.98, 0.98, N)
            mpm = rng.uniform(-0.98, 0.98, N)
            U, V = (mpp + mpm) / 2, (mpp - mpm) / 2
            _, _, C = deficit_terms(U, V)
            if C < worstN:
                worstN, argN = C, (N, mpp.copy(), mpm.copy())
    print(f"  random product laws N=2..4 (6000 draws): min C_bwd = {worstN:.3e}")
    ok2 = worstN > -1e-12
    if ok1 and ok2:
        print("  no violation found: C_bwd >= 0 on the entire M=2 product-law")
        print("  family tested (exhaustive at N=1). Status: numerically")
        print("  supported conjecture, no proof claimed.")
    else:
        print(f"  COUNTEREXAMPLE FOUND: {argN if not ok2 else arg}")
    print(f"  ... survey completed; violation found: {not (ok1 and ok2)}")
    return True    # the survey itself succeeding is the check


if __name__ == "__main__":
    r1 = check1()
    r2 = check2()
    r3 = check3()
    print("=" * 60)
    allok = r1 and r2 and r3
    print(f"ALL CHECKS: {'PASS' if allok else 'FAIL'}")
    sys.exit(0 if allok else 1)
