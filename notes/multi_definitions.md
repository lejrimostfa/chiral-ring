# Multi-Ring Model — Definitions and Conjectures (Paper 2)

Companion to `chiral_ring_spec.md` and `backward_process.md` (Paper 1). Paper 1
is the foundation and is **not modified**. All Paper-2 code lives in `src/multi/`.

Status: definitions fixed before coding (work-order step 1). Conjectures C1–C5
are stated here as targets; verdicts live in `results_multi.md`.

**Terminology (interpretation-agnostic).** "Branch" = a component of the
einselected decomposition of the joint ring state; the chirality tuple
`S = (s_1,…,s_M)` labels it. Any *branch-relative* quantity is a conditional
expectation in the ordinary Bayesian sense (condition on the classical branch
variable, which is conserved). No Everettian commitment is made or needed; an
observer-internal reading (a record-holder conditions on the facts they hold) is
sufficient and is the one used in Paper 1, §5.

---

## 1. Setup

- `M` chiral rings (principal case `M = 2`), each a Paper-1 3-site ring. Ring `i`
  carries a conserved chirality `s_i = ±1` (eigenvalue of its current operator
  `Î_i`, values `±I₀`, `I₀ = √3 J`). Because `[H_S, Î_i] = 0` for every ring,
  the joint tuple `S = (s_1,…,s_M) ∈ {±1}^M` is a fixed classical branch
  variable, drawn once at `t = 0` with uniform prior `P(S) = 2^{-M}`.
- Bath: `N` qubits, coupling matrix `g_{i,n} ≥ 0` (`i` ring, `n` qubit). Default
  bath state `|+x⟩^{⊗N}`.
- Conditional on `S`, qubit `n` sees a pure-dephasing rotation
  `U_n(S,t) = e^{-i a_n(S) t σ_z}`, with the **branch phase rate**

      a_n(S) = Σ_i s_i g_{i,n}.                                        (1)

  Its conditional state is `|φ_n(S,t)⟩ = e^{-i a_n(S) t σ_z}|+x⟩`.

Everything below is a *closed form* derived from `|φ_n(S,t)⟩`. Only posterior
*sampling* (drawing records ξ and forming Bayes posteriors) is Monte-Carlo. If a
quantity seems to require simulating dynamics, the branch structure has been
missed (work-order guard-rail).

## 2. Coherences — the object Γ

Single-qubit overlap between two branches:

    ⟨φ_n(S',t)|φ_n(S,t)⟩ = ⟨+x| e^{-i (a_n(S)-a_n(S')) t σ_z} |+x⟩
                          = cos( (a_n(S) − a_n(S')) t )  ≡ c_n(S,S';t).  (2)

Real, so all branch coherences are real. The joint reduced ring state in the
chiral basis is

    ρ_joint(t) = (1/2^M) Σ_{S,S'} Γ_{S,S'}(t) |S⟩⟨S'|,
    Γ_{S,S'}(t) = Π_n c_n(S,S';t),    Γ_{S,S} = 1.                      (3)

`ρ_joint` is a `2^M × 2^M` matrix; off-diagonal `(S,S')` block equals
`Γ_{S,S'}/2^M`. (This is a *reduced* ring–ring object; the global ring+bath state
stays pure under the unmonitored dynamics.)

**Decoherence-free subspace (DFS).** `Γ_{S,S'}(t) = 1 ∀t` iff
`c_n(S,S';t)=1 ∀n,t` iff `a_n(S) = a_n(S') ∀n`, i.e.

    Σ_i (s_i − s'_i) g_{i,n} = 0   for every qubit n.                   (4)

When all rings share **identical couplings** `g_{i,n} = g_n` (∀i), then
`a_n(S) = g_n Σ_i s_i` depends only on the total chirality `Σ_i s_i`. Any two
branches with the same total chirality are mutually decoherence-free. For `M = 2`
the pair `(+,−)` and `(−,+)` both have total chirality `0`: the coherence
`Γ_{(+,−),(−,+)}(t) = 1 ∀t` (protected), while `Γ_{(+,+),(−,−)}` decays maximally.
This is the analytic content behind conjecture **C1**.

## 3. Records — collision readout and distinguishability

Collision protocol (as in `backward_process.md` §2): qubit `n` interacts for a
window and is measured in `σ_y`, outcome `r_n = ±1`. The `σ_y` expectation in
`|φ_n(S,t)⟩` is

    μ_n(S,t) = ⟨φ_n(S,t)| σ_y |φ_n(S,t)⟩ = sin( 2 a_n(S) t ),          (5)

so the record likelihood is

    P(r_n | S) = (1 + r_n μ_n(S,t)) / 2.                                (6)

(For a dedicated collision window Δt the argument is `2 a_n(S) Δt`; we treat the
readout time `t` as the interaction time, the snapshot convention of Paper-1
`exact.py`. Both are the same one-parameter family.)

Global sign flip `−S = (−s_1,…,−s_M)` gives `a_n(−S) = −a_n(S)`, hence
`μ_n(−S) = −μ_n(S)`: the antipodal branch is the perfectly time-reversed record
source, exactly as in Paper 1 (there `M=1`, `−S = −s`).

**Branch-relative stochastic entropy production** (Paper 1, Eq. 6, lifted to the
joint antipode):

    σ_S[ξ] = Σ_n ln [ P(r_n|S) / P(r_n|−S) ]
           = Σ_n ln [ (1 + r_n μ_n(S)) / (1 − r_n μ_n(S)) ].           (7)

Mean over forward trajectories in branch `S`:

    ⟨σ_S⟩ = Σ_n μ_n(S) ln[(1+μ_n(S))/(1−μ_n(S))]
          = Σ_n D_KL( P(·|S) ‖ P(·|−S) ) ≡ D(S,t) ≥ 0.                 (8)

Define per-qubit quality `ε_n(S) = |μ_n(S)| = |sin(2 a_n(S) t)|`, so
`D(S,t) = Σ_n ε_n ln((1+ε_n)/(1−ε_n))`. For `M=1` this reduces to Paper-1
`record_distinguishability` verbatim.

## 4. Posterior, polarization, marginal arrows

Given a record `ξ = (r_1,…,r_N)`, the Bayes posterior over the `2^M` branches is

    P(S | ξ) ∝ P(S) Π_n (1 + r_n μ_n(S)) = 2^{-M} Π_n (1 + r_n μ_n(S)).  (9)

Per-ring branch estimate and polarization (Paper-1 `P = E|2q−1|` lifted):

    m_i(ξ) = E[s_i | ξ] = Σ_S s_i P(S|ξ),   P_i(t) = E_ξ | m_i(ξ) |.    (10)

Ring–ring record correlation (used in scenario c):

    χ_12(ξ) = E[s_1 s_2 | ξ] = Σ_S s_1 s_2 P(S|ξ),  E_ξ[χ_12].          (11)

### 4.1 Marginal entropy production σ_i  (the C2 subtlety)

`σ_S` in (7) is a *joint* arrow: it uses the full branch `S` and its global
antipode `−S`. A per-ring arrow `σ_i` must compare `s_i = +1` against `s_i = −1`
while the other chiralities are *unknown to ring i's record channel*. There is a
genuine modelling choice here; we state the options and fix one.

- **Option A — marginal-likelihood ratio (CHOSEN).** Marginalise the other
  chiralities out of the likelihood with the uniform prior:

      L_i^{±}(ξ) = Σ_{s_{-i}} 2^{-(M-1)} Π_n (1 + r_n μ_n(s_i=±1, s_{-i})),
      σ_i[ξ] = ln [ L_i^{+}(ξ) / L_i^{-}(ξ) ] · sign convention below.  (12)

  To make `σ_i` branch-relative (an arrow *for the branch actually realised*,
  matching (7)), we evaluate it in branch `S` and orient it by `s_i`:

      σ_i^{(S)}[ξ] = s_i · ln [ L_i^{+}(ξ) / L_i^{-}(ξ) ].               (13)

  For `M=1`, `L_i^{±}` has no sum and (13) collapses to (7): consistent.

- **Option B — plug-in / conditional.** Condition on the *true* other
  chiralities and flip only `s_i`:
  `σ_i^{plug}[ξ|S] = Σ_n ln[P(r_n|S)/P(r_n|S^{(i)})]`, `S^{(i)} =` S with `s_i`
  flipped. This is clean but **not** additive-decomposable against a shared bath
  in the way C2 needs: `Σ_i σ_i^{plug}` double-counts qubits coupled to several
  rings through the joint antipode, and its deficit against `σ_S` has no clean
  information-theoretic name.

**Choice and justification.** We adopt Option A/(13). Rationale: (i) it reduces
to Paper 1 at `M=1`; (ii) each `σ_i` is then the log-likelihood-ratio of an
*observable* marginal (the record's information about `s_i` alone, integrating
out what it cannot resolve), so `⟨σ_i⟩ = ` the marginal KL distinguishability
`D_i`, the natural per-ring analogue of (8); (iii) the deficit below acquires a
clean name. Option B is recorded here and in a code comment as the alternative;
we do not pursue it (C2 stop-loss: one session).

### 4.2 The additivity deficit (C2)

    Δ(t) ≡ ⟨σ_S⟩ − Σ_i ⟨σ_i⟩ = D(S,t) − Σ_i D_i(t).                   (14)

With disjoint baths (no qubit couples to two rings) the joint likelihood
factorises over rings, the marginals are exact, and `Δ = 0` (scenario a — an
*exactness test*, not a result). With a shared bath, qubits carry joint
information about several chiralities; the marginal channels lose the
cross-ring part. **Conjecture C2:** `Δ(t) ≥ 0`, and equals the *multi-information*
(total correlation) that the record ensemble carries among the rings,

    Δ(t) =? I_multi[records] = Σ_i H(σ_i-channel) − H(joint σ-channel),  (15)

estimated on the sampled trajectories (symmetric binning). C2's job is to test
(15) numerically and report the residual; the sign `Δ ≥ 0` is the firm part.

## 5. Mirror preparation (backward boundary condition)

Following `backward_process.md`: a *mirror* bath for ring `i` is prepared in the
Θ-conjugate of the **final** state of a forward evolution of duration `T`. In
closed form the per-qubit record quality **decreases**:

    ε_n^{mir}(t) = | sin( 2 a_n(S) (T − t) ) |,   t ∈ [0,T].            (16)

At `t=0` the mirror bath is "fully charged" (`ε_n = |sin(2 a_n(S) T)|`); as `t`
grows the records are **consumed**, `ε_n → 0` at `t=T`. Operationally: draw
records from the *backward* law `P_bwd(r_n|S) = (1 − r_n ε_n^{mir})/2` (Paper 1,
Eq. 3) while still scoring σ with the forward law (7). Then

    ⟨σ_i^{mir}⟩(t) = − Σ_n ε_n^{mir}(t) ln[(1+ε_n^{mir})/(1−ε_n^{mir})]
                   = − D_i^{mir}(t) ≤ 0.                                (17)

So a mirror ring runs its arrow **backwards**: entropy is destroyed, records
consumed. This is a legitimate boundary condition, not a dynamical law change;
the resulting global initial state is ring–bath **entangled** and is documented
as such (it is the Θ-image of a forward-evolved, hence entangled, state). Test
(c) of the model module checks (17) at `f=0` to sampling precision.

## 6. Inter-ring entangled initial state (for C3)

Optional ring-pair initial state (instead of the product of on-ring
superpositions):

    |Ψ_rings(0)⟩ = ( |+−⟩ + |−+⟩ ) / √2   (chiral basis, M=2).          (18)

This lives entirely in the zero-total-chirality sector `{(+,−),(−,+)}`. Its
`4×4` joint coherence matrix is read off from Γ (3): the only dynamical
off-diagonal is `Γ_{(+,−),(−,+)}(t)`. The **concurrence** of the reduced
two-ring state `ρ_12(t)` is, for this X-form state, the closed form

    C(t) = |Γ_{(+,−),(−,+)}(t)|    (normalised as below),               (19)

i.e. concurrence is *literally the protected coherence*. Define the **relative
distinguishability** of the two occupied branches,

    D_rel(t) = Σ_n |sin( (a_n(+,−) − a_n(−,+)) t )| · ln[(1+·)/(1−·)]
             = record distinguishability of (+,−) vs (−,+),            (20)

which controls how fast a bath could tell the two branches apart. **Conjecture
C3:** `C(t)` is a monotone decreasing function of `D_rel(t)` — entanglement
decays exactly as records distinguishing the two branches form — and in the DFS
(identical couplings, where `a_n(+,−) = a_n(−,+)`, so `D_rel ≡ 0`) the
concurrence is protected: `C(t) = 1 ∀t`. Slogan: *entanglement survives exactly
where no records form.*

## 7. The five conjectures (targets)

- **C1 — DFS objectivity.** Identical couplings ⟹ the `(+,−)↔(−,+)` coherences
  survive (Γ = 1) while `(+,+)↔(−,−)` decays; the **total** current
  objectifies (redundant records of `Σ_i s_i`), the **individual** chiralities
  never do. Test: `Γ` curves + per-ring `P_i` staying low while total-current
  polarization saturates.

- **C2 — additivity deficit.** `Δ(t) = ⟨σ_S⟩ − Σ_i⟨σ_i⟩ ≥ 0`, conjectured equal
  to the record multi-information (15). Marginal `σ_i` via Option A/(13).
  STOP-LOSS: one session on the formalisation; then measure and report the
  residual as-is.

- **C3 — entanglement lives off the records.** For initial state (18),
  `C(t) = |Γ_{(+,−),(−,+)}(t)|` is a decreasing function of `D_rel(t)`, and
  survives indefinitely in the DFS. "Entanglement survives exactly where no
  records form."

- **C4 — detailed fluctuation theorem.** At trajectory level
  `p(σ)/p(−σ) = e^{σ}`; `σ<0` trajectories exist with exponential suppression.
  Test: symmetric-bin histogram of `σ_S[ξ]` at three times; slope-1 check of
  `ln[p(σ)/p(−σ)]`; fraction of `σ<0` vs `D`. (Integral form `⟨e^{−σ}⟩=1` is the
  Paper-1 sanity check.)

- **C5 — arrow frustration (open question).** Ring 1 on a fresh forward bath,
  ring 2 on a mirror bath (records to consume), sharing a fraction `f` of qubits.
  Shared qubits cannot simultaneously accumulate and consume distinguishability.
  Competing hypotheses, none privileged: (a) the majority arrow wins;
  (b) both `σ_i` attenuate symmetrically; (c) an `f`-dependent mixed regime. **No
  firm prediction** — the observed regime, described without forcing an
  interpretation, is the result.

## 8. Reproduction and consistency anchors (must hold by construction)

- `M=1`: (2)→`∏cos(2gt)`, (5)→`sin(2gt)`, (6)/(7)/(8) → `exact.py`
  quantities verbatim (machine precision). [model test a]
- DFS (identical couplings): `Γ_{(+,−),(−,+)} ≡ 1`, `D_rel ≡ 0`, `C ≡ 1`.
  [model test b]
- Mirror at `f=0`: `⟨σ_i^{mir}⟩(t) = −D_i^{mir}(t)` (17). [model test c]
- Disjoint baths: `Δ ≡ 0` (14). [scenario a, exactness]
