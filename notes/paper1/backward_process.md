# Backward Process for a Recording Bath — Derivation Note

Companion to `chiral_ring_spec.md`, §5 Tier 2. Status: derived, pending independent hand-check of Eqs. (3), (6), (8).

## 1. Time reversal operator — the choice, made explicit

Take Θ = K_ring ⊗ (iσ_y K)^⊗N, i.e. complex conjugation in the site basis on the ring, physical spin time-reversal on each bath qubit. Then:

- Θ H_S Θ⁻¹ = H_S (real hopping)
- Θ Î Θ⁻¹ = −Î (current is T-odd) ✓
- Θ σ_z Θ⁻¹ = −σ_z, hence Θ H_int Θ⁻¹ = (−Î)⊗(−Σg_nσ_z) = +H_int ✓

**The total dynamics is T-invariant: microreversibility holds.** All time asymmetry lives in states/branches, not in the law. (The naive choice Θ = K on everything makes H_int T-odd and wrecks the construction — this is decision point #1 and the reason the note exists.)

Under Θ: |±⟩ ↔ |∓⟩ (branch flips), |+_x⟩ → −|−_x⟩, σ_y outcomes flip sign (a recorded current direction is itself T-odd data).

## 2. Trajectory protocol: collision-model readout

For clean Markovian trajectory statistics, adopt the collision version: bath qubit n, prepared in |+_x⟩, interacts with the ring during [t_{n−1}, t_n] (Δt = t_n − t_{n−1}) and is projectively measured in the σ_y basis at t_n, outcome r_n = ±1. The unconditional system dynamics keeps the product decoherence-factor structure; measuring the bath never back-acts on ρ_S (Î is conserved).

Since [H_S, Î] = 0, the chirality s = ±1 is a **fixed classical branch variable** drawn once at t=0 with prob ½ each. Conditional on s, qubit n ends at angle 2s g_n Δt from +x in the equatorial plane, and

    P(r_n | s) = (1 + r_n s ε_n)/2,   ε_n ≡ sin(2 g_n Δt)          (1)

The entire forward process is: hidden fair coin s, then independent biased readings. Record quality = ε_n.

    P_fwd(ξ) = ½ Σ_s Π_n (1 + r_n s ε_n)/2,   ξ = (r_1,…,r_M)      (2)

## 3. Backward process

Θ-conjugated protocol: bath qubits prepared in Θ|+_x⟩ ∝ |−_x⟩, evolution ΘUΘ⁻¹ (motion reversal; in backward-branch label s' this reads e^{+is'g_nτσ_z}), measurement with the Θ-imaged detectors **keeping their forward labels** (the Θ-image of the forward r-detector is the physical |−r⟩_y state, still labeled r; the T-oddness of the recorded datum is carried by the trajectory conjugation ξ̃: r̃_n = −r_n, reversed order, used in Eqs. 4 and 6 — NOT by relabeling the apparatus). [Convention fixed 2026-07-03 after numerical verification caught the original "flipped labels" phrasing as inconsistent.] Direct computation (verified to machine precision):

    P_bwd(r̃_n | s') = (1 − r̃_n s' ε_n)/2                            (3)

with r̃_n the backward-read outcome and s' the backward branch. The T-conjugate of forward record ξ is ξ̃ = reversed order, r̃_n = −r_n.

## 4. Result 1 — the unconditional entropy production vanishes identically

    σ[ξ] = ln P_fwd(ξ) − ln P_bwd(ξ̃)
         = ln ½Σ_s Π(1+r_n s ε_n)/2 − ln ½Σ_{s'} Π(1+r_n s' ε_n)/2 = 0   ∀ξ    (4)

(The double sign flip r→−r, s→−s in (3) restores (2) after summing over the symmetric branch prior.) **This is not a bug.** Pure dephasing dissipates nothing; the dynamics is T-invariant; the branch prior is symmetric. The ensemble has no arrow — precisely the spontaneous-symmetry-breaking structure anticipated in spec §9. Any claimed nonzero unconditional σ in the baseline is a coding error; Eq. (4) is a free exactness test for the numerics.

## 5. Result 2 — the arrow is branch-relative, and it equals record information

Decision point #2: what is the physically meaningful backward comparison? An observer *inside* branch s holds s as a fact — records, memories, the Mlodinow–Brun sense of "memory" are all conditioned on s. The operationally accessible reversal is the Θ-conjugated protocol **within the same branch** (the observer cannot conjugate the fact s; Θ would map their world to the other branch). Define the branch-relative stochastic entropy production:

    σ_s[ξ] ≡ ln P(ξ | s) − ln P_bwd(ξ̃ | s)                          (5)

Using (1) and (3):

    σ_s[ξ] = Σ_n ln [ (1 + r_n s ε_n) / (1 − r_n s ε_n) ]            (6)

Its average over forward trajectories in branch s:

    ⟨σ_s⟩ = Σ_n ε_n ln[(1+ε_n)/(1−ε_n)] = Σ_n D_KL( p(r|s) ‖ p(r|−s) ) ≥ 0   (7)

**Interpretation — the central analytic claim of the paper:** the branch-relative entropy production is exactly the accumulated Kullback–Leibler distinguishability between the branch and its time-reverse, i.e. the information content of the records about the chirality. In this model, irreversibility and record-formation are not merely correlated — under the branch-relative definition they are the *same quantity*. Fluctuation theorem: ⟨e^{−σ_s}⟩ = Σ_ξ P_bwd(ξ̃|s) = 1 holds exactly (numerical sanity check).

## 6. Result 3 — all observables collapse onto one variable

Define the accumulated distinguishability D(t) = Σ_{n: t_n ≤ t} D_KL,n. Then:

- **Entropy production:** ⟨σ⟩(t) = D(t) exactly, by (7).
- **Polarization:** conditional on record ξ, the system state has chiral weights (q, 1−q), q = posterior P(s=+1|ξ), so ⟨Î⟩_ξ = I₀(2q−1) and P(t) = E|2q−1| = posterior concentration. By large-deviation theory the concentration rate is governed by the same per-record divergences (Chernoff/KL family). Rises on the same schedule as D(t).
- **Redundancy:** the classical mutual information between s and a fragment's outcomes is Σ_{n∈F}[1 − h((1+ε_n)/2)], monotone in the same ε_n; the quantum I(S:F) of spec §3 upper-bounds it and shares its timescale.
- **Frank's O (spec §5bis):** a symmetrized KL between forward and echo record distributions — same family, directly comparable to (7).

**Consequence for the paper's framing:** in the baseline, the co-transition is a theorem, not a numerical discovery. The numerics' job becomes (i) confirming the identity and the finite-size sharpening, and (ii) **variant B** (H_E ≠ 0, degrading records), where the branch variable is no longer conserved, Eq. (4) fails, unconditional σ ≠ 0, the identity (7) breaks, and the real question — which of redundancy, polarization, irreversibility leads or lags when records decay — becomes open and is the publishable result.

## 7. Hand-check list (do before coding)

1. Verify Θ H_int Θ⁻¹ = +H_int with Θ = K_ring ⊗ (iσ_yK)^⊗N, and that Θ=K-only fails (§1).
2. Re-derive (1) and (3) from the single-qubit rotations; signs of ε_n.
3. Confirm (4) by writing out a 2-record example (M=2) explicitly.
4. Confirm ⟨e^{−σ_s}⟩ = 1 from (6) analytically for M=1.
5. Check the claimed posterior identity ⟨Î⟩_ξ = I₀(2q−1) (off-diagonals of ρ_S|ξ don't contribute since Î is diagonal in the chiral basis).

## 8. Caveats to carry into the paper

- **Variant B microreversibility (added after run 3):** with H_E = Σω_nσ_x, the §1 operator fails (σ_x is odd under iσ_yK). The correct antiunitary is Θ = K_ring ⊗ (σ_x K)^⊗N: it flips σ_z (H_int even ✓), preserves σ_x (H_E even ✓), and fixes |+_x⟩. Microreversibility holds in variant B too; the derivation of §§4–5 carries over with Helstrom-optimal ε_n(t) = √(1−|c_n(t)|²).

- Eq. (4)'s vanishing relies on pure dephasing + symmetric prior. A biased initial superposition (unequal chiral weights) gives unconditional σ ≠ 0 already in baseline — worth one extra figure (bias sweep).
- The branch-relative definition (5) is a choice; defend it via the observer-internal argument (Mlodinow–Brun) and by noting that the Θ-conjugated-branch alternative gives identically zero (Eq. 4 conditional analogue), i.e. carries no information — the arrow *is* the branch-relativity.
- Collision protocol vs. persistent records: redundancy (spec §3) is computed on the unmonitored model, σ on the monitored one; both share the same unconditional S-dynamics. State this explicitly to preempt referee confusion.
