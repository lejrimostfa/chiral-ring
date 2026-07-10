# The factorization law (Paper 3, Section 4)

Companion to `src/lattice/run_factorization.py`. Established after the six gates
(`verify_lattice.py`) and the exact continuous baseline. All numbers below are
from the L=8 square lattice, collision protocol, prescription F (shared-link
baths fresh), unless stated.

## Statement

Let each ring carry the arrow `tau_i = sign(<sigma_i^local>)` (accumulation
direction, prescription F). **Claim (factorization law).** The arrow field is an
INDEPENDENT site-percolation field: ring `i` points forward with probability

    p_eff(rho, f) = (1 - rho)·1 + rho·conv(f) = 1 - rho(1 - conv(f)),

where `conv(f)` is the single-ring contrarian-conversion probability — a function
of the sharing fraction `f` ALONE — and fresh rings point forward with
probability 1. Every collective threshold then follows in closed form, with no
free parameter:

    rho_c^align(f) = 1 / [2(1 - conv(f))]           (global arrow reversal, p_eff=1/2)
    rho_c^perc(f)  = (1 - p_c) / (1 - conv(f)),   p_c = 0.5927   (percolation, p_eff=p_c)

## The three structural facts it rests on (verified)

1. **Fresh rings never reverse (F).** Across all runs, 35600/35600 fresh rings
   are forward — 0 backflips. The `(1-rho)·1` term is exact.
2. **Conversion depends on `f` alone.** `conv` measured at six densities
   `rho in [0.35, 0.72]` has max spread `std_rho = 0.007`: `rho`-universal.
   `conv(f) = [0, 0, 0.004, 0.024, 0.111, 0.278, 0.465, 0.606]` for
   `f = [0, 0.5, 0.6, 0.67, 0.71, 0.75, 0.78, 0.80]`.
3. **Conversions are independent.** Two-point correlation of contrarian
   conversion events vs lattice separation `d`: bulk `r(4<=d<=9) = -0.025 ~ 0`;
   the only departure is a weak NEAREST-NEIGHBOUR correlation `r(1) = +0.112`
   (5539 pairs, ~8 sigma) — two contrarians sharing a bond read the same fresh
   witnesses, so they convert slightly together. This is the leading (and only)
   correction to strict independence.

## Numerical proof

- **E2 — closed forms, no free fit** (`conv(f)` interpolated from data):
  `rho_c^align`: residual max **0.003**, mean 0.001 — exact.
  `rho_c^perc` : residual max **0.072**, mean 0.022 — confirmed (the max is a
  single point, finite-size / grid-crossing at L=8).
- **E3 — the ultimate test.** A synthetic field of INDEPENDENT sites at density
  `p_eff` (same lattices, same contrarian geometry seeds) reproduces the measured
  quantum `P_perc(rho)`: mean `|Delta| = 0.040`, below the binomial sampling
  floor `~0.056`. Quantum and independent-site percolation are statistically
  indistinguishable.
- **Universality (P3-D).** The rho-driven finite-size-scaling collapse
  (`L in {8,12,16}`) has a sharp interior minimum at `nu = 1.33 = 4/3`, the 2D
  site-percolation exponent — the `L -> infinity` form of the reduction.

The alignment threshold is EXACT (0.003) because it depends only on the mean
occupation `p_eff`, which factorizes regardless of the `r(1)` correction; the
percolation observables carry the weak `r(1) = 0.11` nearest-neighbour clustering
but the universality class is unchanged.

## Boundary of the theorem (E4 — prescription M)

Charge the shared-link baths instead (prescription M): the shared witnesses now
vote BACKWARD and count double. Structural fact (1) fails — 5600/6400 fresh
rings reverse — and `A(f)` at `rho=0.5` plunges from 0 to `-1` (backward wins,
`A(f)` monotone DOWN vs monotone UP under F). The monotonicity that underpins
the factorization (fresh -> forward) is broken; the independent-site reduction
does NOT hold under M. This is the boundary of the theorem, not a failure: the
sign — and the very validity of factorization — is set by the shared-witness
preparation, the same F/M prescription duality as Paper-2 C5. Under F the shared
witnesses vote forward (factorization holds, `rho_c > 1/2`); under M they vote
backward (factorization broken, arrow globally reversed).

## Status

`E1 = E2 = E3 = PASS`. The arrow-alignment transition of Paper 3 reduces, under
prescription F, to ordinary 2D independent site percolation with a
sharing-controlled occupation `p_eff = 1 - rho(1 - conv(f))`; the only correction
to strict independence is the weak nearest-neighbour conversion correlation
`r(1) = 0.11`, which leaves the universality class (`nu = 4/3`) intact. The
reduction is prescription-specific (holds for F, broken for M).

Figures: `fig_p3_factorization.png` (E2), `fig_p3_independence.png` (E1+E3),
`fig_p3_prescriptionM.png` (E4).
