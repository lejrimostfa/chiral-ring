# Paper Skeleton — Chiral Ring / Local Arrow Fixation

Working title (primary):
**"Fixation of a time-reversal-odd observable by decoherence: record distinguishability as the local arrow of time in an exactly solvable model"**

Alternates: "The arrow ratchet: classical records make local T-symmetry breaking permanent" (punchier, riskier); "Co-emergence of objectivity and irreversibility in quantum Darwinism" (safer, blander).

Target: Physical Review A (first choice — the content is now standard quantum-information methodology with exact results); Foundations of Physics (fallback); Entropy (second fallback). Length target: 12–16 pages two-column PRA format, 9 figures (trim to 7 for submission).

---

## Abstract (draft, ~180 words)

The microscopic laws are time-reversal invariant, yet local observers record a fixed arrow of time. We study how decoherence fixes the orientation of a time-reversal-odd observable in an exactly solvable model: a three-site chiral ring, whose degenerate current-carrying doublet forms a T-odd pointer basis, coupled to a bath of recording qubits. The model admits closed forms for all bipartite entropies, and its two-branch structure reduces trajectory statistics to classical inference. We prove that the unconditional stochastic entropy production vanishes identically — the arrow is spontaneously broken, existing only relative to a branch — and that the branch-relative entropy production equals exactly the accumulated Kullback–Leibler distinguishability D(t) of the environmental records. Numerically, chirality fixation locks to decoherence (ratio 0.95 across two decades of bath size), while Darwinian redundancy lags by a factor of four: the observable fixes before it becomes objective. Branch polarization collapses onto a universal function of D across protocols and couplings, with no hysteresis under record degradation — until records are measured, converting the revocable arrow into a ratchet. Classicalization, not decoherence alone, makes the local arrow permanent.

---

## Section plan

### 1. Introduction (1.5 pages)
Frame: two literatures that do not talk — quantum Darwinism (objectivity via redundancy: Zurek, Riedel–Zurek–Zwolak) and stochastic thermodynamics of irreversibility (entropy production as time-reversal breaking: Maes–Netočný, Landi–Paternostro). The question: do classicality, objectivity, and the local arrow emerge together, and in what order? Position against: Mlodinow–Brun (records and the psychological arrow), Chen's Entanglement Past Hypothesis, Guff–Shastry–Rocco (opposing arrows in open systems), and Frank's Tlalpan interpretation (co-emergence conjectured as a phase transition — our model provides the standard-QM quantitative baseline it must be compared to). One paragraph, no more, on the broader emergent-spacetime motivation (Van Raamsdonk, Cao–Carroll); keep it as motivation, not claim. State the four results at the end of the intro.

### 2. Model (1.5 pages) — from spec §1–2
Ring Hamiltonian, current operator, T-odd pointer doublet, recording bath, exact solvability, the two-branch structure. Emphasize why 3 sites (minimal chirality) and why a T-odd pointer observable is the right object for arrow-of-time questions (fixing it breaks T locally — unlike position or spin-z pointers in standard models). Figure: model schematic (TO MAKE: ring + bath cartoon with the two circulation branches). State the correct antiunitary Θ for baseline and variant B (note §1, §8) and microreversibility.

### 3. Entropy production for a recording bath (2 pages) — from note §2–5
Collision protocol; the classical hidden-branch reduction (Eq. 1–2); backward process (Eq. 3); Theorem 1: unconditional σ ≡ 0 (Eq. 4) with the spontaneous-symmetry-breaking reading; the branch-relative definition and its defense (observer-internal, Mlodinow–Brun); Theorem 2: ⟨σ_s⟩ = D(t) (Eq. 7). IFT as corollary. This section is the paper's formal contribution.

### 4. Numerical results (3 pages)
4.1 Validation: machine-precision agreement analytic/brute-force (one line + repo pointer). Partial information plots reproduce the canonical shape → **Fig: fig1_pip**.
4.2 Co-transition and hierarchy: decoherence, polarization, redundancy on one axis → **Fig: fig2 (snapshot) or fig4 (collision — prefer, monotone)**; histograms unimodal→bimodal → **Fig: fig3**. Timescale table: t_P/t_S = 0.95 ± 0.01 over N = 8–128; t_R/t_S ≈ 4.
4.3 Scaling: all timescales ∝ N^{−1/2}, relative width constant → **Fig: fig5**. Explicit statement: no relative sharpening; the transition is not critical in time.
4.4 Universality: data collapse P(D) across (N, g, protocol) → **Fig: fig6**. This is the central figure; the identity of Sec. 3 explains it (σ, P, R are all functionals of the record ensemble {ε_n}).

### 5. Degrading records: revocability and the ratchet (1.5 pages)
Variant B, still exactly solvable; rise-and-fall of redundancy (contact with Riedel–Zurek–Zwolak) → **Fig: fig7**. Hysteresis test: P rides the master curve both directions, residual < 0.04 → **Fig: fig8**. The ratchet: measured records freeze ε, P becomes monotone → **Fig: fig9**. Punchline paragraph: coherent records make the local arrow revocable; classicalization makes it permanent. Connect to quantum Darwinism's core claim: objectivity = redundant *classical* proliferation.

### 6. Discussion (1.5 pages)
(i) What is and is not shown: spatial ring, T-odd observable, branch-relative arrow — not emergence of time; the circularity caveat (dynamics presupposes locality). (ii) Relation to Tlalpan/QTI: our smooth universal sigmoid in D is the standard-QM baseline; a genuinely sharp threshold in amplification would be a discriminating signature — we sharpen its testability. (iii) Experimental cousin: superconducting flux qubits (clockwise/counterclockwise superpositions) — the chiral doublet is realizable hardware; sketch what measuring P vs D would require. (iv) Open: interacting environments (redundancy destruction), biased priors (unconditional σ ≠ 0 — one figure in appendix), and the graph/causal-order generalization (one sentence, future work).

### 7. Methods / Appendices
A: closed-form entropies (rank-2 structure). B: proof details of Eq. (4) and (7). C: Helstrom optimality of the readout basis; variant-B ε_n(t). D: numerical parameters, seeds, repo link.

---

## Figure inventory and trimming plan
Keep: schematic (new), fig4, fig3, fig5, fig6, fig8, fig9 (7 total). Move fig1, fig2, fig7 to appendix or supplement.

## Pre-submission checklist
1. Hand-check the five items of note §7 (STILL PENDING — do not skip).
2. Rerun everything with error bars: n_traj ≥ 50k, 10 coupling realizations, shaded bands. On the M1 Max via the repo.
3. Densify fig5 (add N = 256, 512 — analytic, cheap) and add the g-sweep.
4. Investigate the constant 0.037 hysteresis residual (suspected sampling floor at low D; raise n_traj and check it shrinks).
5. Literature lock: Scholar crawl on Riedel–Zurek 2012 forward-citations × "entropy production", and Maes–Netočný forward-citations × "Darwinism". One hour. If a collision, reposition before writing Sec. 3.
6. Co-author/endorsement outreach with figs 6, 8, 9 attached (Zwolak, Riedel, Paternostro; or NTNU quantum group for proximity).
7. Draft Sec. 3 first (the theorems), then 4–5 (results write themselves from figures), intro last.
