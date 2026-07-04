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

## C2 — Collective (synergistic) irreversibility — **REINTERPRETED & CONFIRMED (extensive synergy); literal multi-information identity INFIRMED**

Figures: `fig_multi_c_deficit.png`, `fig_multi_d_sweep_f.png`,
`fig_multi_c2_exact.png`, `fig_multi_deficit_scaling.png`.
Scripts: `run_scenarios.py c`, `analyze_c2.py`, `run_deficit_scaling.py`.

**Firm part (confirmed):** `Δ(t) = ⟨σ_S⟩ − Σ_i⟨σ_i⟩ ≥ 0` everywhere; exactly 0
for disjoint baths (machine precision); monotone increasing in the shared
fraction (`0 → 0.33 → … → 43.7` nats as `f: 0 → 1` at N=64, t=2).

**Exact decomposition (enumeration, verified to 1e-14 across seeds, times, f):**

    Δ = I(s₁;s₂|ξ) + C_bwd

- `I(s₁;s₂|ξ) ≤ ln 2` is the genuine record-induced multi-information — the
  object the original conjecture pointed at. It IS inside Δ, but bounded.
- `C_bwd` is the antipodal cross co-information (the pointwise record synergy
  `i(ξ;s₁;s₂)` scored on the true branch but averaged over records drawn from
  the antipodal branch). It is a cross-entropy-type term with no bounded
  Shannon reading — this is the extensive part.
- Universal constant: the bounded fraction `I(s₁;s₂|ξ)/Δ → 0.1249 ≈ 1/8` as
  `t → 0`, independent of disorder seed and of f (tested f ∈ {0.3, 0.6, 1.0},
  4 seeds); it decays to ~0.01 at large t. The conjectured multi-information
  never dominates the deficit.

**Scaling (the synergistic-irreversibility reading, confirmed):**

    f=1, t=1.0, N ∈ {16,32,64,128}:  Δ = 0.6876·N − 0.385,  R² = 0.99995

Δ is **linear-extensive in the total number of shared records** (the intercept
is negligible against the ~88-nat scale at N=128). At fixed N=64, sweeping the
shared count k mixes shared and private records and shows onset curvature
(`R² ≈ 0.95`, β < 0): the deficit is *not* additively decomposable qubit-by-qubit
in a mixed bath —

**Per-record decomposition: NOT exact.** With
`δ_n = KL_n^joint − Σ_i KL_n^(i)` (closed form, exactly 0 on private qubits),
enumeration at N=10 gives relative residuals `|Δ − Σ_n δ_n|/Δ` of 0.20–0.85
with a sign change: `Σ_n δ_n` is a closed-form extensive *proxy*, not an
identity. The residual is the intra-marginal record-correlation term
`Σ_i [Σ_n KL_n^(i) − D_i]`. No new closed form is claimed in the notes.

**Verdict:** the deficit is a **collective (synergistic) irreversibility** —
an entropy production carried only by the joint record channel, extensive in
the number of shared records (R²=0.99995), containing the record
multi-information as its bounded (≤ ln 2, small-t fraction 1/8) component. The
literal identity `Δ = multi-information` stays infirmed. Option A for σ_i
unchanged (stop-loss respected).

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

## C5 — arrow frustration — **RESOLVED: sharing-controlled arrow reversal; the preparation of the shared witnesses decides which arrow wins**

Figures: `fig_multi_e_frustration.png`, `fig_multi_e_crossing.png`,
`fig_multi_c5_prescriptions.png`, `fig_multi_fstar_stability.png`.
Scripts: `run_scenarios.py e`, `run_c5_robustness.py`, `run_fstar_stability.py`.
Prescriptions derived in `multi_definitions.md §6bis` before coding.

**Setup.** Ring 1 forward (fresh bath), ring 2 mirror (charged for T=1.5),
sharing fraction f. A shared qubit cannot simultaneously accumulate and consume;
its preparation is a boundary-condition choice, tested three ways (exclusive
qubits fixed: ring-1-private fresh, ring-2-private mirror-charged):
**F** shared fresh (T_c=0) · **M** shared mirror-charged (T_c=T) · **H** shared
half-charged (T_c=T/2). f=0 limits reproduce `σ_1=+D_1`, `σ_2=−D_2^mir` exactly
(sampler validated bit-for-bit against the step-4 code).

**Robustness result (N=64, t=0.5, R=10):**

| prescription | shared witnesses | f\*(σ₂) | f\*(σ₁) | which arrow flips |
|---|---|---|---|---|
| F | fresh (accumulate)      | **0.761** | none | the **mirror** ring's |
| M | charged (consume)       | none | **0.103** | the **forward** ring's |
| H | half-charged (consume, weaker) | none | **0.683** | the **forward** ring's, late |

The phenomenon — a sharing-controlled reversal of exactly one ring's arrow —
exists under **all three** prescriptions. Which arrow cedes, and at what f, is
set by the shared witnesses' preparation: fresh witnesses drag both arrows
forward (the mirror ring flips); charged witnesses drag both backward (the
forward ring flips); a weaker charge flips the forward ring only at much larger
sharing. This is outcome (ii) of the work order: *the preparation of the shared
witnesses decides which arrow wins.*

**Statistical stability (prescription F, per-realization f\*, R=20):**

| | N=32 | N=64 | N=128 |
|---|---|---|---|
| t=0.5 | 0.763 ± 0.039 | 0.765 ± 0.034 | 0.748 ± 0.040 |
| t=1.0 | 0.128 ± 0.036 | 0.114 ± 0.018 | 0.103 ± 0.014 |

20/20 realizations show a crossing in every cell; all std ≤ 0.04 < 0.1, so the
word **threshold** (not "crossover") is licensed by the pre-registered
criterion. f\* is N-independent at fixed t and self-averaging (std shrinks with
N); it moves with t (0.76 → 0.11) as the mirror charge is consumed.

**Verdict:** not hypothesis (a) "majority wins", not (b) "symmetric
attenuation" — regime (c), now sharpened: a **prescription-dependent, disorder-
stable threshold** at which the arrow of the ring whose witnesses are
outnumbered reverses. The robust, prescription-independent facts: (i) exactly
one arrow flips, never both; (ii) the flip is controlled by f; (iii) the loser
is determined by the shared-witness preparation, not by which ring is "forward"
or "mirror".

---

## Data collapse (step 5) — collapse SURVIVES sharing, BREAKS in the DFS — now with disorder error bars

Figures: `fig_multi_collapse.png`, `fig_multi_collapse_errors.png`.
Scripts: `collapse_multi.py`, `run_collapse_errors.py`.

Pooled fit (step 5, 1008 generic points):

    GENERIC collapse (shared f=0..1 + random): RMS residual = 0.0200,
        max = 0.0795 (bootstrap 95% CI on max [0.0735, 0.0795])
    identical (DFS): RMS = 0.393   <- breaks the collapse

Per-realization master curves (T4, R=10 — each realization builds its own
generic master curve and is scored against it):

    generic RMS = 0.0179 ± 0.0037   (max over realizations 0.0245)
    DFS     RMS = 0.4008 ± 0.0097   (min over realizations 0.3796)
    worst-case separation  min(DFS) / max(generic) = 15.5×

**The number that decides it:** the dichotomy is disorder-stable — the worst
DFS realization sits 15.5× above the worst generic one; the two distributions
do not approach. The break is not noise: in the DFS the marginal
distinguishability `D_i` accrues while `P_i` stays pinned at chance (individual
chirality is unrecordable), so `(D_i, P_i)` runs along `P_i ≈ 0.5` off the
curve. Collapse-breaking is itself the decoherence-free-subspace signature,
tying step 5 back to C1 and C3.

---

## Summary table

| Conjecture | Verdict | Deciding number |
|-----------|---------|-----------------|
| C1 DFS objectivity | **Confirmed** | P_1=0.50 (indiv, chance) vs P_total=1.0 (identical); P_1=1.0 (random) |
| C2 synergistic irreversibility | **Reinterpreted & confirmed**; literal identity infirmed | Δ=0.688·N, R²=0.99995; exact Δ=I(s₁;s₂\|ξ)+C_bwd; I/Δ→1/8 |
| C3 entanglement off records | **Confirmed** | DFS C≡1 (D_rel≡0); random C↓0 monotone in D_rel |
| C4 detailed FT | **Confirmed** | FT slope = 1.00±0.02; ⟨σ⟩=D; σ<0 fraction 0.31→0.01 |
| C5 arrow frustration | **Resolved: reversal robust, threshold prescription-dependent** | F/M/H: f\*=0.761/0.103/0.683; std(f\*)≤0.04 → threshold |
| Collapse survival | Survives; DFS excepted (disorder-stable) | 0.0179±0.0037 vs 0.4008±0.0097 (15.5× worst case) |

Negative results are results: the literal C2 identity fails (and the exact
decomposition says precisely how), the per-record decomposition of Δ fails
(residual documented), and the collapse break in the DFS is a positive
diagnostic of where records fail to form.

---

## What changed vs the first pass

1. **C2 upgraded from "category error" to an exact statement.** The first pass
   compared Δ against the posterior total correlation (bounded) and stopped at
   "extensive vs bounded, cannot match". The re-analysis found the exact
   decomposition `Δ = I(s₁;s₂|ξ) + C_bwd` (machine precision, enumeration), a
   universal small-t bounded fraction `1/8`, and the clean extensivity law
   `Δ = 0.688·N` (R²=0.99995). The conjectured multi-information is *inside* Δ;
   it just never dominates. Renamed: collective (synergistic) irreversibility.
2. **C5's f\* = 0.76 was indeed prescription-dependent — the worry was
   founded.** Under M (charged shared witnesses) it is the *forward* ring that
   flips, at f\* = 0.103; under H at 0.683. What survives all prescriptions:
   exactly one arrow flips, controlled by f, and the loser is chosen by the
   shared-witness preparation. The single-prescription "mirror flips at 0.76"
   statement of the first pass is now the F-row of a table.
3. **"Threshold" vocabulary is now licensed by measurement**, not assumption:
   std(f\*) ≤ 0.04 < 0.1 across 120 realization-cells (20 per cell, all with
   crossings), N-independent and self-averaging.
4. **The collapse dichotomy carries error bars**: 0.0179±0.0037 vs
   0.4008±0.0097, worst-case separation 15.5× — previously a single pooled
   number each.
5. **One negative result added**: the per-record decomposition `Δ = Σ_n δ_n` is
   not exact (residual 0.20–0.85, sign-changing); `Σ_n δ_n` remains a useful
   closed-form extensive proxy (exactly zero on private qubits).
