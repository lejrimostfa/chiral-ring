# Arrow definition and its reconciliation with Paper-2 C5 (Paper 3)

Companion to `src/lattice/`. Fixed before any physics run (sprint P3-A).

## The two candidate arrows

For a ring `i` monitored by a bath, two sign-valued quantities compete for the
name "arrow of time":

1. **Inference-ratio sign** (Paper 2, C5): `sign ⟨σ_i⟩`, where `σ_i` is the
   branch-relative entropy production (Option-A marginal log-likelihood ratio).
   `⟨σ_i⟩ > 0` means the records *fix* the chirality (positive distinguishability
   accumulated relative to the antipode); `⟨σ_i⟩ < 0` means the accessible
   records are net time-reversed (a mirror/consuming bath).

2. **Accumulation direction** (Paper 3, adopted): 
   `τ_i = sign( d⟨σ_i^local⟩/dt )` — is ring `i`'s accessible record information
   *growing* (forward) or *being consumed* (backward) right now?

`σ_i^local` uses only ring `i`'s accessible witnesses (its private bath plus
every shared link touching it), marginalizing the neighbour chiralities. Gate
G4 shows the local marginal agrees in sign with the global marginal (100% of
sites, mean correlation 0.97), so **the local marginal is the production
definition** — it is `O(2^{deg+1})` per site and therefore scales to a lattice,
whereas the global `2^L` marginal does not.

## Why accumulation, not the inference ratio

The inference ratio measures *how much has been fixed*, a state variable; the
arrow should be a *rate*. A bath can hold a large positive `⟨σ_i⟩` (much
distinguishability already recorded) while that stock is being consumed
(`d⟨σ_i⟩/dt < 0`). The physically contagious quantity — the thing shared
witnesses transmit between rings — is the *direction of accumulation*, so `τ_i`
is defined as its sign. This is the quantity swept in the contagion / percolation
studies (P3-B, P3-C).

## Reconciliation (they agree where it matters, and differ knowably)

At `f = 0` (no shared witnesses) the two definitions **coincide in sign** for
the canonical setups:
- fresh forward ring: `⟨σ_i⟩ = +D_i > 0` and `d⟨σ_i⟩/dt > 0` → both give **+1**;
- charged mirror ring: `⟨σ_i⟩ = −D_i < 0` and its records are consumed,
  `d⟨σ_i⟩/dt < 0` → both give **−1**.

Gate **G6** verifies this f=0 limit exactly against the Paper-2 C5 numbers
(`⟨σ_forward⟩ = +D`, `⟨σ_mirror⟩ = −D`, machine precision). So Paper 3 does not
contradict Paper 2; it *refines* the arrow from a state sign to a rate sign.

The two can only disagree in a **transient**: a bath whose stock `⟨σ_i⟩` is still
positive but already decreasing (just past its accumulation peak). Paper-2 C5
would still call it "forward" (positive fixation); Paper 3 calls it "turning
backward" (`τ_i = −1`). This window is exactly where the instantaneous protocol
produced spurious non-monotonicity (recurrences), which is why Paper 3 runs the
**collision protocol**.

## Collision protocol (kills the recurrences)

Production uses fixed-window collision records, `ε_n = sin(2 a_n(S) Δt)` with `Δt`
a constant (Paper-1 `src/` convention): each witness leaves a *permanent* record,
so accessible distinguishability accumulates **monotonically** and `τ_i` is
well-defined. The instantaneous protocol (bias `sin(2 a_n(S) t)`, oscillatory in
`t`) is kept only for the non-regression gates G1–G4. Consequently the sharing
fraction is swept **structurally**, via the number of shared witnesses `Ns`
(`f = Ns / (Np + 2 Ns)` per bond), never via the coupling strength — strong
shared coupling reintroduces folding of `|sin|` and pollutes `τ`.

## Status

Gates G1–G6 all PASS (`src/lattice/verify_lattice.py`). The engine reproduces:
the Paper-2 decomposition theorem on the lattice law (G1), Paper-1 Theorem 2 at
f=0 (G2), brute-force state-vector evolution (G3), local≡global marginals (G4),
exact collision additivity (G5), and the Paper-2 C5 f=0 limit (G6). No figure is
produced until this holds — it does.
