# Exact reduction — analytic propositions (Paper 3, Section 4 upgrade)

Companion to sprint P3-F. Goal: upgrade the factorization law from "established
numerically (E1–E3)" to "Propositions A–B proven, Proposition C verified
symbolically + exhaustively". Written before the sprint; the sprint verifies,
it does not rediscover.

## Setting (fixed by the model, restated for the proofs)

- True branch: S_true = (+1, ..., +1). Contrarian status is a BATH preparation
  (private witnesses charged, Tc = T), never a chirality flip.
- Prescription F: every shared-bond witness is fresh (Tc = 0), regardless of
  the types of the two rings on its bond.
- Collision protocol: witness n interacts in one window, is measured once, its
  record r_n is permanent. Given the full branch S, records at distinct
  witnesses are INDEPENDENT (product law) — the conditional-independence
  structure of Papers 1–2, inherited by the lattice engine.
- Accessible set of ring i: xi_i = { records of i's private witnesses } ∪
  { records of witnesses on bonds incident to i }.
- Arrow: tau_i = sign( d<sigma_i^loc>/dt ), with sigma_i^loc the marginal
  log-likelihood ratio built from xi_i only, neighbour chiralities
  marginalized with uniform prior.

KEY STRUCTURAL FACT used by every proof below: tau_i is a deterministic
functional of (xi_i ; couplings G_{i,·} and incident-bond couplings). Nothing
else enters.

## Proposition A (exact 1-dependence ⇒ fact 3, sharpened)

Statement. Under S_true and prescription F, the field {tau_i} is a
1-dependent random field: for any two rings i, j at graph distance d(i,j) ≥ 2,
tau_i and tau_j are INDEPENDENT (exactly, not approximately). The only
possible correlation is between nearest neighbours, mediated by their shared
bond witnesses.

Proof. For d(i,j) ≥ 2 the incident-bond sets of i and j are disjoint, hence
xi_i ∩ xi_j = ∅. Records at distinct witnesses are independent given S (which
is fixed = S_true), and the witness preparations are deterministic functions
of the (fixed) contrarian placement; so xi_i ⟂ xi_j, and tau_i, tau_j are
functions of independent inputs. ∎

Falsifiable sharpening (NEW, test in P3-F): the theory predicts r(d) = 0
IDENTICALLY for every d ≥ 2 — not just "small in the bulk". The measured
r(2) = 0.046 and r(3) = 0.013 must be pure sampling noise. Test: error bars
on r(2), r(3) (pair-count based, bootstrap over realizations); PASS iff
compatible with zero at 2 sigma. The nearest-neighbour r(1) = +0.112 is the
only permitted nonzero value, and its mechanism is explicit: two adjacent
contrarians share Ns fresh witnesses, so their conversion functionals share
Ns input records.

## Proposition B (conv is a function of f alone ⇒ fact 2)

Statement. Under S_true and prescription F, P(tau_i = -1 | i contrarian) is
independent of rho and of the placement of the other contrarians; it depends
only on the ring's own witness portfolio — Np charged private witnesses and
deg·Ns fresh shared witnesses, i.e. on f (and the fixed coupling
distribution). Hence conv(f) is rho-universal by construction, and likewise
the fresh-ring flip probability depends on f alone.

Proof. The law of xi_i is determined by (i) the couplings of i's incident
witnesses and (ii) their preparations. Under F, the preparation of every
incident SHARED witness is "fresh" irrespective of the neighbour's type; the
preparation of the private witnesses is set by i's own type. The chiralities
entering a_n(S_true) on incident bonds are all +1 regardless of contrarian
placement. So the law of xi_i — hence of tau_i — carries no dependence on
rho or on the other rings' types. ∎

Consequence recorded for the manuscript: the measured std_rho = 0.007 of
conv across densities is a CHECK of the engine, not evidence for the law —
the law is exact, and the residual 0.007 is the sampling floor.

## Proposition C (fresh never flips ⇒ fact 1) — the one needing care

Statement (conjectured identity, to verify). Under S_true, prescription F,
and dt in the production window (all biases positive: a_n(S_true) > 0 and
sin(2 a_n dt) > 0 for every incident witness), the expected increment of
sigma_i^loc per new record is strictly positive for a fresh ring; hence
<sigma_i^loc> is strictly increasing and tau_i = +1 with probability 1.

Why it is NOT automatic (the honest subtlety). sigma_i^loc uses the marginal
model P(xi|s_i) = 2^{-deg} sum over neighbour chiralities, while the sampling
law is the single branch Q = P(·|S_true). The per-record increment is
E_Q[ ln P(r|xi, s_i=+) / P(r|xi, s_i=-) ], an expectation of a log-ratio of
posterior predictives under a law that is one COMPONENT of the mixture — not
a KL divergence, so positivity is not generic. It should hold here because
Q is the maximal-bias component (all neighbours +1) and all biases are
positive, making P(r=+1|xi, s_i=+) ≥ P(r=+1|xi, s_i=-) pointwise; a full
proof plausibly goes through a likelihood-ratio-ordering (stochastic
dominance) argument on the posterior predictives.

Verification plan (P3-F, symbolic — the Paper-2 1/8-law method):
  C1. sympy, deg=1 (one neighbour, generic couplings g_p, g_s, generic Np=1,
      Ns=1, generic dt): expand the per-record increment exactly; check
      positivity identically in the parameters (or exhibit the condition).
  C2. sympy, deg=2, generic couplings: same. If the increment is a sum of
      manifestly nonnegative terms, promote C to "proven at deg ≤ 2,
      verified exhaustively beyond".
  C3. Exhaustive numeric: random couplings, deg=3,4, 10^4 draws — minimum of
      the increment (target: > 0, machine-zero margin like Paper-2 C_bwd).
  C4. The exhaustive production count (35600/35600 forward) is then the
      final corroboration, not the claim.

## What the manuscript gains (integration map)

If A, B pass trivially (they are proofs, only the r(2)/r(3) noise test is
empirical) and C1–C3 land:

- Theorem 1 (factorization) is restated with facts (2)–(3) PROVEN and fact
  (1) proven-at-small-deg + verified: the wording "exact reduction" in the
  title is then fully earned.
- Section 4.2 becomes "Two proofs and one verified identity" — the three
  structural facts with their one-paragraph proofs inline.
- E1 changes role: from evidence to CHECK, with the sharpened prediction
  r(d≥2) ≡ 0 tested (new error bars on r(2), r(3)).
- The epistemic table in the Discussion moves the factorization from
  "Numerical, strong" to "Exact (A, B) / verified identity (C)".
- The r(1) correction gets its exact mechanism sentence: Ns shared input
  records between adjacent conversion functionals — and (optional, nice)
  a prediction: r(1) should DECREASE as Np grows at fixed f (shared inputs
  diluted); one cheap run can test it.

## Status

A, B: proven above (review the two proofs, then they go verbatim into the
manuscript). C: conjectured identity + verification plan. Nothing here is
established until verify_lattice.py's gates pass in the same session and
C1–C3 run green.

## VERDICTS (session Claude, engine lattice_p3.py, gates 4/4 PASS in-session)

- LEMMA (new, exact): the local marginal FACTORIZES OVER BONDS,
  sigma_i^loc = sum_priv sigma_p + sum_bonds sigma_b, because neighbours are
  marginalized independently and witnesses partition by bond. Verified
  against local_marginal_sigma at 6.7e-16 (L=3, deg=2). Consequence:
  Proposition C for ANY coordination reduces to a single bond.
- C1 (K=0->1, symbolic): increment = b * ln[(2+b)/(2-b)] >= 0, closed form,
  manifestly positive (sympy).
- C1 (K=1->2, symbolic grid): min = +4.4e-09 on 400x400 over (0,1)^2,
  attained only at the boundary b_new -> 0 (zero-information record), the
  correct equality locus.
- C3 (randomized): min increment +4.8e-07 over 20000 configs, K <= 6;
  full private+two-bond sequences: min +6.7e-05 over 2000 sequences.
- Proposition C: VERIFIED IDENTITY (closed form at one record; positive
  everywhere tested beyond). Open analytic item: closed-form proof for
  arbitrary K per bond.
- Propositions A, B: proven above (no change).
- Script: verify_C1.py (ship into src/lattice/, rerun as a gate in-repo).

## F0 verdict (partial, decisive where it matters)

- 2D CONFIRMED: the measured f grid {0.5, 0.6, 0.67, 0.71, 0.75, 0.78, 0.80}
  reproduces deg*Ns/(Np+deg*Ns) EXACTLY for Ns=2..8 (Np=8, deg=4).
  Manuscript Eq. (1) correct for every quantitative claim (all 2D).
- 1D OPEN AND AT STAKE: the reported range f~0.36-0.40 cannot exist on the
  per-ring grid (Np=3 gives 0, 0.4, 0.571, ...). Likely run_chain.py labeled
  f = Ns/(Np+2Ns). AUDIT run_chain.py; if confirmed, relabel the 1D axis
  per-ring and REVISIT the manuscript sentence comparing the 1D threshold
  to the pair f* ~ 0.40 (per-ring it may sit near 0.73).

## F2 verdict (post-sprint P3-F, periodic torus engine)

Proposition A's sharpened prediction r(d>=2) = 0 IDENTICALLY is CONFIRMED
with error bars: r(2) = -0.004 +/- 0.008, r(3) ~ 0 (torus engine, commit
9acc8ca). The old bulk value -0.025 (open-boundary engine, no bars) is
superseded. Note for the record: the boundary fix also moved beta/nu from
0.27 to 0.107 (~exact vs 0.104) and E3 from 0.040 to 0.022 -- the open
boundaries were the dominant systematic, exactly as a 1-dependent field
with edge rings of reduced degree would predict (edge rings had deg<4,
hence a locally different f and conv).
