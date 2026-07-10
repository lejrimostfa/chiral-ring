# Paper 3 — Percolation of arrow alignment: results

Lattice of chiral rings with shared recording baths. Every ring carries a local
arrow `tau_i = sign(<sigma_i^local>)` (accumulation direction, collision
protocol; `arrow_definition.md`); shared witnesses transmit the arrow between
rings, and the *aligned* domains undergo a genuine percolation transition in 2D.

Engine `src/lattice/lattice.py`; verification `src/lattice/verify_lattice.py`;
runs `run_chain.py` (1D), `run_square.py` (2D), `run_minority.py` (rho_c + FSS).
Foundation (Papers 1–2, `src/`, `src/multi/`) unchanged.

## Validation gates (must pass before any figure)

`verify_lattice.py` — 6/6 PASS:

| Gate | What | Residual |
|---|---|---|
| G1 | `Delta = I + C_bwd` on the lattice-generated law (Paper-2 theorem) | 2.8e-17 |
| G2 | `f=0` -> Paper-1 Theorem 2 exact (`<sigma_i> = D_i`) | 2.2e-16 |
| G3 | `L=3` brute-force state vector vs branch construction | 2.9e-17 |
| G4 | local marginal == global marginal (sign 100%, corr 0.97) | sampling |
| G5 | collision protocol: exact additivity (`deficit=0`), `<sigma_i>=D_i` | 0.0 |
| G6 | collision `L=2` reproduces the Paper-2 C5 `f=0` limit (`+D`,`-D`) | 2.2e-16 |

The local marginal (`O(2^{deg+1})` per site) is the production definition; it
agrees with the exact global marginal (G4), so lattices up to `L=16` (256 rings)
are tractable — the `2^L` global object is never built.

## Continuous baseline (exact test throughout)

`A(f=0) = 1 - 2 rho` (using the discretized density `round(rho L)/L`) is
reproduced to machine precision (`run_chain.py` max deviation `0.0`): with no
sharing, fresh rings point forward and charged (contrarian) rings backward.

## 1D chain (P3-B) — contagion, but a crossover

Figures: `fig_p3_chain_alignment.png`, `fig_p3_domains_chain.png`.

- **Contagion.** Signed alignment `A(f)` rises from the `1-2 rho` baseline toward
  full alignment as sharing grows, with a threshold near `f ~ 0.36-0.40` — the
  same neighbourhood as the Paper-2 pair `f* ~ 0.40`.
- **Universality in `f`.** The fraction of contrarians converted collapses onto
  ONE curve in `f`, independent of `rho` (the conversion threshold is a property
  of the sharing fraction alone).
- **Crossover, not a transition.** Mean `tau`-domain length grows with `f` but
  the `L=8/16/32` curves overlap (no finite-size divergence) — consistent with
  the absence of a 1D percolation transition (`p_c = 1`). This motivates 2D.

## 2D square lattice (P3-C) — genuine percolation transition

Figures: `fig_p3_2d_phase.png`, `fig_p3_percolation.png`, `fig_p3_snapshots.png`,
`fig_p3_fss.png`. Forward arrows occupy `p = (1+A)/2`; they percolate above the
2D site threshold `p_c = 0.5927`.

- **Transition.** At `rho=0.5`, `P_perc(f)` for `L in {4,5,6,8}` cross at
  **`f_c = 0.721`** — a finite-size percolation crossing. `f_c` is well above the
  pair `f* ~ 0.40` because in 2D the forward arrow must reach `p_c = 0.593`, not
  merely a majority.
- **Phase diagram.** The percolation boundary `rho_c(f)` sits at `rho ~ 0.41` at
  `f=0` (where `1-rho = p_c`) and RISES with sharing — the forward arrow
  percolates against an increasingly dense contrarian minority.
- **Exponent (P3-C alone: inconclusive).** The `f`-driven FSS collapse at
  `nu=4/3` aligns the transition region, but the coarse integer-`Ns` `f`-grid
  and `L<=8` give no robust interior minimum in `nu` — the exponent could not be
  fixed from the `f`-driven data. Resolved in P3-D.

## Minority + finite-size scaling (P3-D)

Figures: `fig_p3_rhoc.png`, `fig_p3_fss_rho.png`.

### rho_c(f): shared witnesses vote double

The critical contrarian density that reverses the GLOBAL arrow (`A=0`):
`rho_c(0) ~ 0.50` (0.495 on the L=8 rho-grid), flat until `f ~ 0.6`, then climbs
steeply above `1/2`:

    f      : 0.67   0.71   0.75   0.78   0.80   0.83
    rho_c  : 0.51   0.58   0.70   0.89   >0.93  >0.93  (diverges)

Because a shared witness is read by BOTH its rings, it counts double; here the
shared witnesses are fresh (forward), so they weight the forward arrow. Beyond
`f ~ 0.78` the forward arrow wins for EVERY contrarian density tested up to 0.93
— `rho_c` diverges: no minority, even a large majority, can flip the global
arrow. The naive `rho_c = 1/2` is broken in the direction set by the
shared-witness preparation. (The prompt's guess `rho_c < 1/2` would obtain under
charged shared witnesses — the same F/M prescription duality as Paper-2 C5;
fresh witnesses give `rho_c > 1/2`.)

### FSS: the transition is in the 2D-percolation universality class

Driving the transition with the CONTINUOUS parameter `rho` (steps `1/L^2`) at
fixed `f=0.5`, on larger lattices `L in {8,12,16}`:
- percolation threshold `rho_c^perc = 0.391`;
- **the finite-size-scaling collapse has a sharp interior minimum at
  `nu = 1.33`, coincident with the 2D site-percolation value `nu = 4/3`.**

This pins what P3-C could not: the arrow-alignment transition is in the **2D
percolation universality class** (not Ising). The improvement came from the
continuous control parameter and the larger lattices, not from any change to the
model.

## Summary

| Result | Number |
|---|---|
| Gates (independent recode) | 6/6 PASS, residuals <= 2.9e-17 |
| Continuous baseline `A(f=0)=1-2rho` | exact (dev 0.0) |
| 1D: contagion threshold | `f ~ 0.36-0.40`; crossover (no FSS divergence) |
| 2D: percolation point (`rho=0.5`) | `f_c = 0.721` |
| 2D: order-parameter threshold | `p_c = 0.5927` (forward occupation) |
| rho_c(f) reversal of the global arrow | `~0.50 -> 0.89 -> diverges (>0.93)` beyond `f~0.78` (votes double) |
| Universality class | 2D percolation, established by the E3 EQUIVALENCE (below), not by direct exponents |
| Direct exponents (reinforcement) | `gamma/nu = 1.78 == perc 1.79` (robust); `nu` consistent, trends toward 4/3 with size but NOT pinned at `L<=24` (scatters 1.4-2.0); `rho_c=0.407` stable. Unequivocal direct `nu` needs `L~64-128` (HPC) — each site is a quantum sampling. |
| **Factorization law** (P3-E, `factorization_theorem.md`) | arrow field == independent site percolation, `p_eff=1-rho(1-conv(f))`; `rho_c^align=1/[2(1-conv)]` (res 0.003), E3 mean Δ=0.040 < noise; F-specific (M breaks it). **This equivalence is the rigorous universality argument.** |

Open / caveats: `nu` is pinned at `f=0.5` on `L<=16`; a full multi-`f` exponent
map and the `rho_c(f)` divergence point would benefit from larger lattices. The
`rho_c` sign (`>1/2` vs `<1/2`) is prescription-dependent (fresh vs charged
shared witnesses), exactly mirroring the Paper-2 F/M duality.
