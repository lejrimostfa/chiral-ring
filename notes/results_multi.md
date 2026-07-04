# Multi-Ring Model — Results (Paper 2)

Verdicts on conjectures C1–C5 from `multi_definitions.md`, with the number that
decides each. Configuration unless noted: M=2 rings, N=64 bath qubits, collision
readout, R=10 disorder realizations, n_traj=20k. Foundation (Paper 1, `src/`,
`notes/backward_process.md`) unchanged. All figures regenerable with a single
command (`python3 src/multi/run_scenarios.py <key>` or `collapse_multi.py`).

**Validation gate (must pass before reading results).**
- `model_multi.py` unit tests PASS: (a) M=1 reproduces `exact.py` to machine
  precision; (b) DFS exact; (c) mirror `f=0` gives `⟨σ⟩ = −D^mir` within 5 SEM.
- `brute_multi.py` (full state vector, N=8) PASS: reduced coherence == Γ
  (err 3.7e-16), per-branch bias == sin(2 a_n(S) t) (5.6e-16), concurrence ==
  |Γ_(+-,-+)| (8.4e-15), mirror coherence == ∏cos(2g(T−t)) (1.1e-15).
- Control scenario (a), disjoint baths: additivity deficit, correlation and
  multi-information all zero to machine precision (max 1.1e-14). This is the
  exactness test, not a result.

---

## C1 — DFS objectivity of the total current — **CONFIRMED**

Figures: `fig_multi_b_dfs_coherence.png`, `fig_multi_b_objectivity.png`.

With identical couplings (`g_{1,n}=g_{2,n}`), the zero-total-chirality coherence
is protected and the max-total coherence decoheres:

    |Γ_{(+,-),(-,+)}|  min over t = 1.0000   (protected, DFS)
    |Γ_{(+,+),(-,-)}|  min over t = 0.0000   (decoheres)

Objectivity at t=3.0:

| coupling | P_1 (individual) | P_total (recordable sector) |
|----------|------------------|------------------------------|
| identical (DFS) | **0.501** (chance) | **1.000** (objectifies) |
| random          | **1.000** (objectifies) | — |

**The number that decides it:** individual polarization stays at chance (0.501)
under identical couplings but saturates (1.000) under random couplings, while the
total current objectifies (1.000) in both. The bath records only the total
chirality `Σ_i s_i`; the individual chiralities never objectify in the DFS.

---

## C2 — additivity deficit = multi-information — **INFIRMED (literal); firm part holds**

Figures: `fig_multi_c_deficit.png`, `fig_multi_d_sweep_f.png`.

Firm, confirmed part: `Δ(t) = ⟨σ_S⟩ − Σ_i⟨σ_i⟩ ≥ 0` everywhere (exactly 0 for
disjoint baths; positive and monotone increasing in the shared fraction f:
`0 → 0.33 → 2.29 → … → 43.7` nats as `f: 0 → 1`).

Refuted part: the identity `Δ = record multi-information` **fails**, and it must,
for a structural reason:

    Δ(final) = 44.35 ± 5.48 nats     (extensive; grows with shared-qubit count)
    any multi-information of M=2 chirality bits ≤ M ln2 = 1.386 nats   (bounded)
    Δ / (posterior total correlation)  ≈ 2.5×10^3   (median over t)

Δ is an **extensive excess distinguishability** (joint minus marginal KL over the
records), unbounded in N. Any genuine multi-information of the M chirality bits
is bounded by `M ln2`. They cannot be equal. (Note: at large t the posterior
total correlation itself → 0 as posteriors sharpen, while Δ keeps growing — the
mismatch widens.) The marginal `σ_i` uses Option A (marginal-likelihood ratio,
fixed in `multi_definitions.md §4.1`); under this choice the deficit measures the
record information that is joint-only, not a bounded correlation. **Verdict:**
`Δ ≥ 0` is real and turns on with bath sharing; the "= multi-information"
identity is a category error (extensive vs bounded). Per the C2 stop-loss we do
not re-open the σ_i definition; the negative result stands as reported.

---

## C3 — entanglement lives off the records — **CONFIRMED**

Figure: `fig_multi_c3_concurrence.png` (the C3 figure).

Initial ring-pair state `(|+-⟩+|-+⟩)/√2`; concurrence `C(t) = |Γ_{(+-),(-+)}|`.

    identical (DFS): D_rel(t) ≡ 0  (max 0.0e+00),  C(t) ≡ 1.0000  (protected ∀t)
    random couplings: C(t) monotone decreasing in D_rel(t), C(final) = 0.0000
                      (monotone-in-D_rel check = True)

**The number that decides it:** in the DFS `D_rel ≡ 0` and `C ≡ 1` (entanglement
survives indefinitely); with generic couplings `C` is a monotone decreasing
function of the relative distinguishability `D_rel`, reaching 0 exactly as the
records that separate the two branches saturate. *Entanglement survives exactly
where no records form.*

---

## C4 — detailed fluctuation theorem — **CONFIRMED**

Figure: `fig_multi_f_fluctuation.png` (n_traj = 4×10^5).

Symmetric-bin test of `ln[p(σ)/p(−σ)]` vs σ (well-sampled bins only):

| t | ⟨σ_S⟩ | D | fraction σ<0 | FT slope |
|---|-------|---|---------------|----------|
| 0.10 | 0.718 | 0.719 | 0.313 | **1.006** |
| 0.20 | 2.875 | 2.876 | 0.199 | **1.003** |
| 0.35 | 8.807 | 8.814 | 0.125 | **1.019** |

**The number that decides it:** the detailed-FT slope is 1.00 ± 0.02 across
times (theory: exactly 1). `⟨σ_S⟩ = D` at every time (Paper-1 identity, lifted).
`σ<0` trajectories exist and are exponentially suppressed as D grows (fraction
0.31 → 0.12, and → 0.012 already by t=1.0): the arrow-reversing trajectories are
present with `e^{σ}` weight, exactly as predicted.

---

## C5 — arrow frustration — **OPEN; observed regime = f-dependent with a threshold**

Figures: `fig_multi_e_frustration.png`, `fig_multi_e_crossing.png`.

**Modelling choice (stated, not derived).** Ring 1 forward (fresh bath), ring 2
mirror (charged for T, records to consume), sharing fraction f. A shared qubit
cannot simultaneously be fresh (ring 1 is writing it) and mirror-charged
(ring 2 wants to consume it). *Prescription F:* shared qubits are fresh and read
forward by both rings; only ring-2-private qubits carry the mirror charge. Other
prescriptions (shared qubits mirror-charged, or a 50/50 split) would give
different curves; C5 has no privileged answer, so we report what this one
produces. Limits are validated: f=0 gives `σ_1=+D_1` (forward) and `σ_2=−D_2^mir`
(mirror), matching the model tests.

Observed (σ in nats, final time; T=1.5):

| f | σ_1 (forward) | σ_2 (mirror) | σ_2 zero-crossing in t |
|---|---------------|--------------|------------------------|
| 0.00 | +40.0 | consumes, → 0 at t≈1.45 | t≈1.45 |
| 0.25 | +38.5 | +17.3 | t≈0.80 |
| 0.50 | +35.7 | +25.4 | t≈0.65 |
| 1.00 | +28.1 | +29.8 | none (accumulates throughout) |

Two effects, neither matching a pure hypothesis:
- The forward ring's arrow **attenuates** with sharing (`σ_1: 40 → 28`): shared
  qubits carry joint information, diluting each ring's marginal arrow.
- The mirror ring's arrow **reverses**: fresh shared qubits inject forward
  entropy production that overwhelms the mirror consumption. At the intermediate
  time t=0.5 the mirror arrow annuls at a sharing threshold

      f* = 0.761   (σ_2 < 0 consuming below f*, σ_2 > 0 accumulating above).

**Verdict:** not hypothesis (a) "majority wins" and not (b) "symmetric
attenuation" — the observed regime is (c) **f-dependent with a threshold f\***:
below f* the mirror ring still runs its arrow backward; above f* the shared fresh
qubits flip it forward. The forward ring meanwhile weakens monotonically. We
report this regime without asserting it is the unique physics — under a different
shared-qubit prescription the threshold would move; the *existence* of a
sharing-controlled reversal is the robust observation.

---

## Data collapse (step 5) — collapse SURVIVES sharing, BREAKS in the DFS

Figure: `fig_multi_collapse.png`. Pooled marginal `(D_i=⟨σ_i⟩, P_i)` over
coupling type, time, ring index, realization (1008 generic points).

    GENERIC collapse (shared f=0..1 + random): RMS residual = 0.0200,
        max = 0.0795 (bootstrap 95% CI on max [0.0735, 0.0795])
    per-config RMS: shared f=0.00..1.00 and random all in [0.010, 0.026]
    identical (DFS): RMS = 0.393   <- breaks the collapse

**The number that decides it:** the Paper-1 co-transition collapse survives bath
sharing at RMS 0.020 (all shared/random configs on one master curve), but the
DFS leaves it at RMS 0.39. The break is not noise: in the DFS the marginal
distinguishability `D_i` accrues while `P_i` stays pinned at chance (individual
chirality is unrecordable), so `(D_i, P_i)` runs along `P_i ≈ 0.5` off the curve.
Collapse-breaking is itself the decoherence-free-subspace signature, tying step 5
back to C1 and C3.

---

## Summary table

| Conjecture | Verdict | Deciding number |
|-----------|---------|-----------------|
| C1 DFS objectivity | **Confirmed** | P_1=0.50 (indiv, chance) vs P_total=1.0 (identical); P_1=1.0 (random) |
| C2 deficit = multi-info | **Infirmed (literal)**; Δ≥0 firm | Δ=44.4 nats ≫ M ln2 = 1.39 nats (extensive vs bounded) |
| C3 entanglement off records | **Confirmed** | DFS C≡1 (D_rel≡0); random C↓0 monotone in D_rel |
| C4 detailed FT | **Confirmed** | FT slope = 1.00±0.02; ⟨σ⟩=D; σ<0 fraction 0.31→0.01 |
| C5 arrow frustration | **Open → regime (c)** | forward σ_1: 40→28; mirror σ_2 flips at f*=0.76 |
| Collapse survival | Survives; DFS excepted | generic RMS=0.020 vs DFS RMS=0.39 |

Negative results are results: C2's identity is false for a clean structural
reason (extensive deficit vs bounded multi-information), and the collapse break
in the DFS is a positive diagnostic of where records fail to form.
