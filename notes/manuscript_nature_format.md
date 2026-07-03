# Classical records make the local arrow of time permanent

**[Author name]**¹

¹ [Affiliation, City, Country]. e-mail: [address]

---

## Summary paragraph (Nature format, ~200 words)

The microscopic laws of physics are invariant under time reversal, yet every observer records a fixed arrow of time. Decoherence theory explains how quantum systems acquire classical properties¹⁻³, and quantum Darwinism explains how those properties become objective through redundant environmental records⁴⁻⁶, but whether classicality, objectivity and the local orientation of time emerge together — and in what order — has remained an open question, studied so far in disconnected literatures⁷⁻¹⁰. Here we show, in an exactly solvable model in which the pointer observable is itself odd under time reversal (the circulation direction of a chiral ring monitored by a qubit bath), that the arrow of time is spontaneously broken: the unconditional entropy production vanishes identically, while the branch-relative entropy production equals exactly the information content of the environmental records. Numerically, the orientation of the time-reversal-odd observable locks to the decoherence time across two decades of environment size, whereas Darwinian objectivity lags fourfold; branch polarization collapses onto a universal function of accumulated record distinguishability, with no hysteresis as coherent records degrade and revive. Only measurement converts this revocable arrow into a ratchet. Classicalization of records, not decoherence alone, is what makes the local arrow of time permanent.

---

## Main

Irreversibility is not written in the dynamical laws. Newtonian, quantum and relativistic dynamics are time-reversal symmetric, and the standard resolutions of this tension — the thermodynamic arrow from special initial conditions¹¹, the psychological arrow from the physics of memory⁷, the decoherent arrow from environmental monitoring³,⁹ — share a common intuition: the arrow of time is carried by *records*. Quantum Darwinism has made the record concept quantitative: a property of a system becomes objective when information about it proliferates redundantly into the environment, so that many observers can read it independently from small environmental fragments⁴⁻⁶. In parallel, stochastic thermodynamics has made irreversibility quantitative: the entropy production of a process equals the log-ratio of the probabilities of a trajectory and of its time reverse¹²,¹³. These two quantifications — of objectivity and of irreversibility — have developed with essentially no contact. Recent proposals that classicality and the arrow of time co-emerge at a threshold, as an amplification-driven symmetry-breaking transition¹⁰, have remained at the level of postulated order parameters. Here we connect the two frameworks in a model simple enough to be solved exactly, and specifically constructed so that the question is sharp: its pointer observable is odd under time reversal, so that fixing it *is* breaking time-reversal symmetry locally.

**An exactly solvable model with a T-odd pointer observable.** The system is a single particle on a three-site ring — the minimal system with a well-defined sense of circulation. Its degenerate excited doublet consists of two states of opposite chirality, |+⟩ and |−⟩, carrying opposite currents; the current operator Î commutes with the ring Hamiltonian and is odd under the antiunitary time-reversal operation, which exchanges the two chiralities. The environment is a bath of N qubits coupled to the current, so that each qubit's precession records the circulation direction: the chiral states are exact pointer states, and the pointer observable is T-odd (Fig. 1). Prepared in the symmetric superposition (|+⟩+|−⟩)/√2 — itself time-reversal invariant — the system decoheres without dissipation, and because the branch structure is exactly two-fold, every bipartite entropy has a closed form and the trajectory statistics of environmental measurements reduces to classical inference on a hidden fair coin: the branch. All numerical results below were validated against these closed forms at machine precision, and the redundancy curves reproduce the canonical partial-information plots of quantum Darwinism⁵ (Extended Data Fig. 1). The model has a direct experimental cousin: superconducting flux qubits, whose clockwise/counterclockwise persistent-current superpositions realize precisely such a T-odd doublet¹⁴.

**The arrow is spontaneously broken and equals record information.** With the correct time-reversal operation — which flips both the current and the bath spins, leaving the total dynamics invariant — we prove two exact results (Methods). First, the unconditional stochastic entropy production vanishes identically, trajectory by trajectory: σ[ξ] = 0 for every measurement record ξ. The ensemble has no arrow, exactly as a ferromagnet above threshold has no magnetization; the time asymmetry is spontaneously broken and exists only relative to a branch. Second, for an observer within a branch — for whom the chirality is a fact that cannot be conjugated away — the branch-relative entropy production is strictly positive and equals, in expectation, the accumulated Kullback–Leibler distinguishability of the records, ⟨σ⟩ = D(t) = Σₙ D_KL[p(rₙ|s) ‖ p(rₙ|−s)]: **the irreversibility of the process is identically the information content of the records about the direction of time-reversal-odd data**. The associated integral fluctuation theorem ⟨e^{−σ}⟩ = 1 holds exactly and was confirmed numerically to one per cent for up to ten records, beyond which the estimator variance grows exponentially, as is generic¹³. The identity was confirmed by trajectory sampling to a relative deviation of 10⁻⁴, limited only by sampling.

**Fixation precedes objectivity.** Tracking three quantities through the transition — the decoherence of the system (entropy of the reduced state), the branch polarization P (the mean magnitude of the conditional current, given the records), and the Darwinian redundancy R_δ (how many independent environmental fragments each supply the chirality) — reveals a strict temporal hierarchy (Fig. 2). The polarization locks to the decoherence time: across bath sizes N = 8 to 128 the ratio of their half-rise times is constant at 0.954 ± 0.012, with polarization marginally leading, since the posterior over branches concentrates as soon as any record discriminates them, while the system's coherence requires the full product of overlaps to die. Redundancy, in contrast, lags fourfold: the chirality is *fixed* long before it is *objective*, consistent with the known continuation of redundancy growth after decoherence saturates⁵. The distribution of conditional currents evolves from unimodal, centred on zero — the time-symmetric phase — to sharply bimodal at the two chiralities (Fig. 2b): the fingerprint of spontaneous T-symmetry breaking at the level of branches. All timescales scale as N^{−1/2}, and the transition does not sharpen relative to its own timescale: in the record variable there is no critical threshold, a point to which we return.

**A universal curve, ridden in both directions.** The deeper structure is that all these quantities are functionals of a single variable. Across bath sizes, coupling strengths and measurement protocols — including a collision protocol with entirely different dynamics — the branch polarization collapses onto one universal sigmoid P(D) of the accumulated record distinguishability (Fig. 3a). This is explained by the exact identity above: entropy production, polarization and redundancy are all controlled by the same per-record divergences. The collapse survives the harshest test we could construct: adding a bath self-Hamiltonian that degrades and revives the records (the model remains exactly solvable), D(t) becomes non-monotone, and the polarization retraces the *same* universal curve on the way down, with residuals below 0.04 in all regimes probed (Fig. 3b). Coherent records make the local arrow fully revocable: when the environment recoheres, the branch orientation unfixes along exactly the path by which it fixed, with no hysteresis.

**The classicalization ratchet.** The revocability disappears the moment records are classicalized. Running the identical degrading-record dynamics, but measuring each bath qubit once — converting its quantum record into a classical, permanent one — the polarization becomes monotone: each measured record is a click of a ratchet that no recoherence can undo (Fig. 4). The contrast isolates the precise role of the quantum-to-classical transition in the arrow of time. Decoherence fixes the orientation of a T-odd observable within a branch; redundant proliferation makes that orientation objective; but only the classicalization of the records makes it *permanent*. In the language of quantum Darwinism, the arrow of time inherits the modal status of the records that carry it.

**Discussion.** Three limitations bound the claims. The ring is spatial, and its T-odd observable is a current: the model demonstrates decoherence-driven local breaking of time-reversal symmetry of the state, co-analysed with the irreversibility of the process — it does not model the emergence of time itself, and, like all emergence-of-classicality arguments, it presupposes a dynamics with a locality structure. The branch-relative definition of entropy production is a choice; we defend it on the grounds that the alternative — conjugating the branch — is unavailable to any observer whose memories are themselves records within the branch⁷, and yields identically zero information. Finally, our universal sigmoid is smooth: standard decoherence plus redundancy produces no sharp threshold in the record variable. This sharpens, rather than undermines, the testability of recent proposals in which the quantum-classical transition is a genuine critical phenomenon¹⁰: any observed threshold-like disappearance of interference at a critical amplification would be a departure from the exactly solvable baseline established here. The flux-qubit realization of the chiral doublet¹⁴ suggests that measuring the polarization–distinguishability curve, and the ratchet contrast, is within reach of existing circuit-QED techniques.

---

## Methods

**Model.** H = H_S ⊗ 1 + (Î/I₀) ⊗ Σₙ gₙσ_z⁽ⁿ⁾ (+ Σₙ ωₙσ_x⁽ⁿ⁾ for degrading records), with H_S the three-site hopping Hamiltonian, Î the current operator, [H_S, Î] = 0, I₀ = √3J. Couplings gₙ drawn uniformly from [ḡ/2, 3ḡ/2]. Bath initialized in |+ₓ⟩^⊗N. The conserved chirality reduces the global state to two branches, each a product over bath qubits; all bipartite reduced states have rank ≤ 2, giving closed-form entropies (binary entropy of (1±|c|)/2 with c the appropriate product of single-qubit overlaps). Brute-force state-vector evolution at N = 8 agrees with the closed forms to 5×10⁻¹⁴.

**Time reversal.** Baseline: Θ = K_ring ⊗ (iσ_yK)^⊗N, under which Î → −Î, σ_z → −σ_z, hence ΘHΘ⁻¹ = H. With the bath self-Hamiltonian, Θ = K_ring ⊗ (σ_xK)^⊗N. Backward protocol: Θ-imaged preparation, motion-reversed evolution, Θ-imaged detectors retaining forward labels; the T-oddness of the record datum enters through the trajectory conjugation ξ̃ (reversed order, r̃ₙ = −rₙ).

**Entropy production.** Collision protocol: qubit n interacts during [tₙ₋₁, tₙ], measured in the σ_y (Helstrom-optimal) basis. P(rₙ|s) = (1 + rₙ s εₙ)/2 with εₙ = sin(2gₙΔt) (baseline) or √(1−|cₙ|²) (degrading records). σ[ξ] = ln P_fwd(ξ) − ln P_bwd(ξ̃) vanishes identically after branch averaging (double sign flip plus symmetric prior). Branch-relative: σ_s[ξ] = Σₙ ln[(1+rₙsεₙ)/(1−rₙsεₙ)], with ⟨σ_s⟩ = Σₙ εₙ ln[(1+εₙ)/(1−εₙ)] = D. All six formal steps were verified by an independent numerical suite (machine precision; repository).

**Observables.** Redundancy R_δ: smallest fragment fraction whose average mutual information with the system reaches (1−δ)S(ρ_S), δ = 0.1, 150 fragment samples per size. Polarization: posterior q(ξ) over branches from Bayes on the sampled records; the conditional current obeys ⟨Î⟩_ξ = I₀(2q−1) exactly (verified against exhaustive enumeration at N = 4); P = E|2q−1|, 3,000–20,000 trajectories. Timescales: interpolated 50% crossings; widths: 10–90% rise.

**Data and code availability.** All code (≈700 lines, Python/NumPy) and generation scripts for every figure are available at [repository]; every figure regenerates from fixed seeds by a single command.

---

## References

1. Zeh, H. D. *The Physical Basis of the Direction of Time* (Springer, 2007).
2. Joos, E. et al. *Decoherence and the Appearance of a Classical World in Quantum Theory* (Springer, 2003).
3. Zurek, W. H. Decoherence, einselection, and the quantum origins of the classical. *Rev. Mod. Phys.* **75**, 715 (2003).
4. Zurek, W. H. Quantum Darwinism. *Nat. Phys.* **5**, 181–188 (2009).
5. Riedel, C. J., Zurek, W. H. & Zwolak, M. The rise and fall of redundancy in decoherence and quantum Darwinism. *New J. Phys.* **14**, 083010 (2012).
6. Zurek, W. H. Quantum theory of the classical. *Entropy* **24**, 1520 (2022).
7. Mlodinow, L. & Brun, T. A. Relation between the psychological and thermodynamic arrows of time. *Phys. Rev. E* **89**, 052102 (2014).
8. Chen, E. K. The decoherent arrow of time and the entanglement past hypothesis. Preprint at arXiv:2405.03418 (2024).
9. Guff, T., Shastry, C. U. & Rocco, A. Emergence of opposing arrows of time in open quantum systems. *Sci. Rep.* **15** (2025).
10. Frank, A. Time symmetry, retrocausality, and emergent collapse. Preprint at arXiv:2508.19301 (2025).
11. Albert, D. Z. *Time and Chance* (Harvard Univ. Press, 2000).
12. Maes, C. & Netočný, K. Time-reversal and entropy. *J. Stat. Phys.* **110**, 269 (2003).
13. Landi, G. T. & Paternostro, M. Irreversible entropy production: from classical to quantum. *Rev. Mod. Phys.* **93**, 035008 (2021).
14. Friedman, J. R. et al. Quantum superposition of distinct macroscopic states. *Nature* **406**, 43–46 (2000).

---

## Figure legends

**Fig. 1 | A time-reversal-odd pointer observable.** [TO MAKE] Three-site chiral ring with degenerate circulation doublet |±⟩; qubit bath coupled to the current; time reversal exchanges the branches. Right: the two conditional bath-qubit trajectories on the Bloch sphere, whose growing distinguishability εₙ(t) constitutes the record.

**Fig. 2 | Fixation precedes objectivity.** **a**, Decoherence S(ρ_S), branch polarization P, record distinguishability D and redundancy R_δ versus time in the collision protocol (N = 120): P locks to S while R lags fourfold; the sampled entropy production ⟨σ⟩ (dots) falls exactly on D. **b**, Distribution of the conditional current across branches: unimodal and time-symmetric before the transition, bimodal at ±I₀ after — spontaneous T-symmetry breaking. **c**, Finite-size scaling: all timescales ∝ N^{−1/2}; the ratio t_P/t_S = 0.954 ± 0.012 across two decades.

**Fig. 3 | One universal curve, no hysteresis.** **a**, Branch polarization versus accumulated distinguishability for six parameter sets and two protocols: collapse onto a single sigmoid. **b**, With degrading and reviving records (bath self-Hamiltonian, three ratios ω/ḡ), the (D, P) trajectory rides the same master curve up and down; maximal residual 0.037.

**Fig. 4 | The classicalization ratchet.** Identical degrading-record dynamics, ω/ḡ = 2: with coherent records the polarization is revocable and oscillates; with each record measured once (classicalized), polarization is monotone — the arrow ratchets and cannot be undone by recoherence.

---

## Notes for revision (delete before submission)

- Fig. 1 schematic remains TO MAKE (vector graphics; Bloch-sphere inset from conditional_states()).
- Author, affiliation, ORCID, repository URL, data-availability DOI: placeholders.
- Error bars: rerun all figures with ≥10 coupling realizations and n_traj ≥ 5×10⁴ (checklist item 2).
- Honest fit assessment: Nature/Science accept ~none of this genre; realistic ladder is Nature Physics (long shot) → PRL (Letter rewrite, 4 pages) → PRA/Quantum. This draft's format converts to PRL with minor surgery: the summary paragraph becomes the abstract, subheads drop, Methods becomes Supplemental Material.
