# Chiral Ring Decoherence Model — Technical Specification

**Goal.** Test numerically whether three quantities transition *jointly* (same timescale, same scaling with coupling) in a minimal model:
1. **Redundancy** of chirality information in the environment (Zurek's quantum Darwinism measure)
2. **Chirality polarization** — fixation of a T-odd observable at the branch/trajectory level
3. **Irreversibility** — stochastic entropy production of the decoherence process

Claim under test: *the local fixation of a time-reversal-odd observable, the objectivity of that observable (redundancy), and the irreversibility of the process co-emerge at the same threshold.*

---

## 1. The system: 3-site chiral ring

Single spinless particle on a 3-site ring (minimal system with well-defined circulation). Hilbert space dim 3, site basis {|0⟩, |1⟩, |2⟩}.

**System Hamiltonian** (ħ = 1, J = 1 sets units):

    H_S = -J Σ_{j=0}^{2} ( |j⟩⟨j+1| + |j+1⟩⟨j| ),   indices mod 3

Eigenstates are quasimomentum states |k⟩ = (1/√3) Σ_j e^{ikj} |j⟩ with k ∈ {0, +2π/3, −2π/3} and energies E_k = −2J cos k. The two states k = ±2π/3 are **degenerate** (E = +J) and carry opposite circulation. Call them |+⟩ and |−⟩ (the chiral doublet).

**Current (chirality) operator:**

    Î = iJ Σ_j ( |j+1⟩⟨j| − |j⟩⟨j+1| )

Properties to verify in code at startup (unit tests):
- [H_S, Î] = 0 (current is conserved → exact pointer states, pure decoherence regime)
- Î|±⟩ = ±I₀|±⟩ with I₀ = √3 J; Î|k=0⟩ = 0
- Under time reversal T (antiunitary; complex conjugation in the **site** basis): T|+⟩ = |−⟩, hence T Î T⁻¹ = −Î. **Î is T-odd.**

**Initial system state:** |ψ_S(0)⟩ = (|+⟩ + |−⟩)/√2.
Verify: T-invariant up to global phase (it maps to (|−⟩+|+⟩)/√2), and ⟨Î⟩ = 0.
The k=0 state is excluded by the initial condition and never populated (current conservation) → effective 2-level chiral doublet, but keep the 3-dim space for cleanliness.

## 2. The environment: N-qubit recording bath

N qubits (target N = 16–24; development at N = 8).

**Interaction ("the field that deforms the loop"):**

    H_int = (Î / I₀) ⊗ Σ_{n=1}^{N} g_n σ_z^{(n)}

**Bath self-Hamiltonian:** H_E = 0 in the baseline (pure decoherence). Variant B adds Σ_n ω_n σ_x^{(n)} to study imperfect records.

**Initial bath state:** |+_x⟩^{⊗N} with |+_x⟩ = (|0⟩+|1⟩)/√2, so each qubit precesses conditionally on the chirality and acquires which-path information.

**Couplings:** two configurations, both required:
- (a) Uniform g_n = g — permutation symmetry, fastest checks against analytics
- (b) Random g_n ~ Uniform[g/2, 3g/2] — genericity; this is the configuration for the paper's main figures

**Total:** H = H_S ⊗ 1 + H_int (+ H_E in variant B).

### Exact solvability (use it)

Since [H_S, Î] = 0 and H_int is diagonal in the chiral doublet, the model is **analytically solvable** in the baseline: conditioned on |±⟩, bath qubit n rotates about z with angular frequency ±g_n. The decoherence factor is

    r(t) = ⟨E_−(t)|E_+(t)⟩ = Π_n cos(2 g_n t)

All reduced density matrices of (S ∪ any fragment F) have closed forms (Zurek-style). **Implement both the analytic expressions and the brute-force state-vector evolution; they must agree to machine precision at N = 8.** Analytics then serve to push redundancy curves to large N cheaply; brute force is required only for variant B and for trajectories.

## 3. Observable 1 — Redundancy R_δ(t)

For fragment F of m bath qubits, compute quantum mutual information

    I(S:F) = S(ρ_S) + S(ρ_F) − S(ρ_SF)

- Partial information plots: Ī(m) = average of I(S:F) over n_samples = 50 random fragments per size m (exact average over all fragments when feasible, m ≤ 6 or uniform couplings)
- Redundancy: R_δ = N / m_δ where m_δ is the smallest m with Ī(m) ≥ (1−δ) S(ρ_S). Use δ = 0.1 (standard) and report sensitivity for δ ∈ {0.05, 0.2}
- Output: R_δ(t) on a time grid; identify t_R = first time R_δ ≥ 2 (onset of genuine redundancy) and the shape of the plateau

Validation target: reproduce the classic antisymmetric partial-information plot with plateau at S(ρ_S) (Riedel–Zurek–Zwolak, NJP 14, 083010, 2012).

## 4. Observable 2 — Branch-level chirality polarization

**Critical subtlety (this is the physics, encode it as written):** the unconditional state never polarizes. By symmetry, ⟨Î⟩_ρS(t) = 0 and the decohered ρ_S = ½(|+⟩⟨+| + |−⟩⟨−|) is itself T-invariant. The T-symmetry breaking is **spontaneous and branch-level**, exactly like magnetization in a ferromagnet ensemble. Therefore polarization must be measured on trajectories/branches:

**Method (baseline, exact):** projectively measure each bath qubit in the σ_y basis at time t (the recording basis for this coupling — verify: conditional bath states differ by rotation about z from |+_x⟩, so information lives in the x–y plane; compute the optimal readout basis explicitly rather than assuming). For each of n_traj = 2000 sampled measurement outcomes ξ, obtain the conditional system state ρ_S|ξ and its chirality ⟨Î⟩_ξ.

**Order parameters:**

    P(t)   = E_ξ [ |⟨Î⟩_ξ| ] / I₀            (mean branch polarization, ∈ [0,1])
    A_T(t) = E_ξ [ D_tr( ρ_S|ξ , T ρ_S|ξ T⁻¹ ) ]   (branch-level T-asymmetry, trace distance)

Also record the full histogram of ⟨Î⟩_ξ/I₀ vs time — the expected signature is a unimodal→bimodal transition (this makes the best figure of the paper).

Identify t_P = time at which P(t) reaches ½ its asymptotic value.

## 5. Observable 3 — Entropy production / irreversibility

Two tiers, both required:

**Tier 1 (cheap, ensemble-level):** for pure decoherence there is no heat flow, and the entropy production of the map equals ΔS_S(t) = S(ρ_S(t)) − S(ρ_S(0)). Report it; note in the paper that for this model it is also the decoherence measure, so Tier 1 alone cannot discriminate the claim.

**Tier 2 (the actual result):** trajectory-level stochastic entropy production. Unravel the dynamics as sequential monitoring of bath qubits (qubit n measured at time t_n on a grid; Markov-chain of outcomes). For each trajectory ξ = (outcomes), compute the log-ratio of forward and time-reversed path probabilities:

    σ[ξ] = ln [ P_fwd(ξ) / P_bwd(ξ̃) ]

with the backward process defined by the time-reversed protocol and T-conjugated states (Crooks / Maes–Netočný construction; cite Maes & Netočný, cond-mat/0211252, for σ = log-ratio as the measure of time-reversal breaking). Deliver ⟨σ⟩(t) and the distribution p(σ); check the integral fluctuation theorem ⟨e^{−σ}⟩ = 1 as a numerical sanity test.

Identify t_σ = time at which ⟨σ⟩ reaches ½ its asymptotic value.

**Design note:** the precise definition of the backward process for a recording bath must be written down in a short accompanying note (one page of algebra) *before* coding Tier 2. This is the one genuinely novel formal step of the project; if it stalls, fall back to the two-time distinguishability proxy D_tr(ρ_S(t+τ), ρ̃_S,bwd(t+τ)) and flag the limitation.

## 5bis. Observable 4 (optional, cheap) — Echo-based record asymmetry O

Bridge to the Tlalpan interpretation (Frank, arXiv 2508.19301), which defines a record-asymmetry order parameter O as the symmetrized KL divergence between forward-protocol and echo-protocol record distributions. In this model:

- Forward protocol: evolve to time t, measure bath qubits → distribution p_fwd(r)
- Echo protocol: evolve to t/2, apply the (imperfect) reversal U†(t/2) with a controlled perturbation ε on the couplings (g_n → g_n(1+ε η_n), η_n ~ N(0,1)), evolve back, measure → p_bwd(r)
- O(t; ε) = ½[D_KL(p_fwd‖p_bwd) + D_KL(p_bwd‖p_fwd)]

In the baseline the echo is analytically computable (Loschmidt-echo structure of a pure-dephasing model). Extract t_O and include it in the co-transition analysis of §6. Purpose: (i) direct commensurability with QTI's proposed order parameter; (ii) tests whether the echo-based and trajectory-based (§5 Tier 2) irreversibility measures agree — they need not, and a discrepancy is itself reportable.

## 6. The co-transition test (the paper's result)

For each (g, N) in the grid:

    g ∈ {0.05, 0.1, 0.2, 0.5, 1.0} (units of J),  N ∈ {8, 12, 16, 20}

extract the three timescales t_R, t_P, t_σ and test:

1. **Coincidence:** are the ratios t_P/t_R and t_σ/t_R constant across the grid (within error bars)?
2. **Scaling:** do all three scale identically with g and N (expected ~ 1/(g√N) for the decoherence timescale in this model — verify against r(t) analytics)?
3. **Sharpness:** does the transition sharpen with N (finite-size scaling of the width of the P(t) rise)? If yes, extract exponents — the phase-transition language becomes defensible.

**Outcomes and their interpretation (decide before running):**
- Three timescales lock together → main claim supported; write the paper as "objectivity, T-odd fixation and irreversibility co-emerge"
- Redundancy lags polarization (expected from Riedel–Zurek: redundancy keeps growing after decoherence saturates) → the *fine structure* is the result: chirality fixes at decoherence time, becomes objective at redundancy time, and σ tracks one of the two. Either way there is a paper; the negative/split result is arguably more informative.

## 7. Variant B (robustness, after baseline works)

Add H_E = Σ_n ω_n σ_x^{(n)}, ω_n ~ Uniform[0, ω_max], ω_max ∈ {0.1g, g}. Records degrade; redundancy rises and falls (compare Riedel–Zurek–Zwolak "rise and fall"). Question: does branch polarization *revert* (recoherence of chirality) when records are erased, and does σ track the reversal? This directly probes "the arrow unfixes when redundancy is lost" — the strongest possible support for the framework if seen.

## 8. Implementation

- **Stack:** Python 3.11+, NumPy, SciPy (sparse where useful), QuTiP optional (prefer custom dense state-vector code — the structure is too special for generic solvers to be efficient). Matplotlib for figures. Numba only if profiling demands it.
- **State vector:** dim 3·2^N complex128; N = 20 → 3.1M amplitudes ≈ 50 MB. Fine on M1 Max. Evolution: exact via the conditional-rotation structure (no Trotterization needed in baseline); Krylov (expm_multiply) for variant B.
- **Entropies:** reduced density matrices by tensor reshaping + partial trace; eigh on ≤ 2^6·3 sized matrices (cap fragment size at m ≤ N/2, sufficient for R_δ).
- **Reproducibility:** single `config.yaml` for all parameters; fixed RNG seeds; every figure regenerable by one CLI command.
- **File layout:**

```
chiral-ring/
  config.yaml
  src/model.py        # H, Î, T operator, initial states, unit tests of §1
  src/exact.py        # analytic r(t), reduced ρ, closed-form I(S:F) (uniform g)
  src/evolve.py       # state-vector evolution (baseline + variant B)
  src/redundancy.py   # partial information plots, R_δ(t)
  src/trajectories.py # bath measurement sampling, P(t), A_T(t), histograms
  src/entropy_prod.py # Tier 1 + Tier 2 σ[ξ], fluctuation-theorem check
  src/scaling.py      # grid runs, timescale extraction, finite-size analysis
  notes/backward_process.md   # the one-page algebra of §5 Tier 2 — write FIRST
  figures/
```

- **Order of work:** unit tests of §1 → analytics vs brute force at N=8 → §3 redundancy (validate vs literature) → §4 trajectories → §5 Tier 2 → §6 grid → §7 variant B. Stop-loss: if §5 Tier 2's backward process cannot be defined cleanly in ~2 sessions, ship the proxy version and note it.

## 9. Honest framing (for the eventual paper — do not skip)

- This is a **spatial** ring with a **T-odd observable**, not a temporal loop. The model demonstrates decoherence-driven local breaking of time-reversal symmetry of the *state*, co-analyzed with irreversibility of the *process*. It does not model the emergence of time itself.
- The unconditional dynamics never breaks T; the breaking is branch-relative (spontaneous-symmetry-breaking structure). Say so explicitly; it preempts the obvious referee objection.
- Must-read/must-cite before writing: Riedel–Zurek–Zwolak (NJP 2012), Zwolak–Zurek (1703.10096), Zurek (Entropy 2022 review, §4.2.5), Maes–Netočný (cond-mat/0211252), Landi–Paternostro (RMP 2021), Mlodinow–Brun (PRE 2014), Guff–Shastry–Rocco (Sci Rep 2025), Chen, Entanglement Past Hypothesis (2405.03418).
- **Tlalpan diff (done, July 2026):** Frank's QTI (arXiv 2508.19301) is an interpretive essay: it proposes order parameters (χ, O, τ_RC) and postulates scaling laws with critical exponents, but derives nothing from microscopic dynamics — no QD mutual information, no redundancy, no entropy production, no solved Hamiltonian. The calculational niche is open. Positioning: this work is the first microscopic quantification of a QTI-style co-transition. Second, adversarial angle: QTI claims a *sharp threshold* discriminates it from standard decoherence (predicted "smooth decay"). If §6.3 finite-size sharpening is observed here — in a fully standard-QM model with no retrocausality — that discriminating prediction is neutralized: threshold-like behavior emerges from decoherence + redundancy alone. State both readings in the discussion.
- Physical cousins to mention in discussion: superconducting flux qubits (clockwise/counterclockwise superpositions), persistent-current rings — the chiral doublet is experimentally realizable hardware.
