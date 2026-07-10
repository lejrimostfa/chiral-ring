# The sharing-fraction convention `f` (Paper 3) — traceability note

Companion to F0 (P3-F). The `f`-convention question is **CLOSED**; this note
records the two conventions actually used in the code and their exact relation so
the manuscript axes and cross-references stay consistent. Nothing here is an open
audit — it is documentation of a decision already made and already reflected in
the manuscript (Sec. 3.2 convention note).

## Two conventions, both live in the repo

The sharing fraction is "how much of a ring's witness content is shared with
neighbours". There are two natural ways to count it, and the 1D and 2D scripts
use different ones:

| Script (dim) | Code (`frac(Ns)`) | Axis label | Convention |
|---|---|---|---|
| `run_chain.py` (1D) | `Ns / (NP + 2*Ns)`, `NP=3` | `f = N_s/(N_p+2N_s)` | **per-bond** |
| `run_square.py` (2D) | `4*Ns / (NP + 4*Ns)`, `NP=8` | `f = 4N_s/(N_p+4N_s)` | **per-ring** |
| `run_factorization.py`, `run_universality.py`, `run_convergence.py` (2D) | `4*Ns/(NP+4*Ns)` | — | **per-ring** |

- **per-ring** (2D): a bulk ring of degree `z` carries `N_p` private witnesses and
  `z·N_s` shared witnesses (each of its `z` incident shared links holds `N_s`).
  The shared fraction of its witness content is `f_ring = z·N_s/(N_p+z·N_s)`; in
  2D bulk `z=4` → `4N_s/(N_p+4N_s)`. This is the physically-primary quantity: it
  is the actual fraction of a ring's records that another ring can also read.
- **per-bond** (1D): the label used on the 1D chain axis counts a single bond's
  worth, `f_label = N_s/(N_p+2N_s)` — the denominator `N_p+2N_s` is the per-ring
  witness count (degree `z=2` in the chain bulk), but the numerator is `N_s`, one
  link's share, not the ring's total `2N_s`.

## Exact relation (1D, bulk degree 2)

    f_ring = 2 N_s / (N_p + 2 N_s) = 2 · f_label .

So the 1D chain axis undercounts the per-ring sharing by exactly a factor 2. This
is the Eq.-(1) conversion referenced in the manuscript.

### Consequence for the 1D onset (already applied in the manuscript)

The contagion onset reported on the raw 1D axis (`f_label ≈ 0.36–0.40`) is a
**per-bond** number. Converted to the per-ring convention used everywhere in 2D:

    f_label ≈ 0.36–0.40  →  f_ring = 2·f_label ≈ 0.67–0.80  (near saturation,
    since f_ring is bounded by 2N_s/(N_p+2N_s) → 1).

More precisely, on `N_p=3` the per-bond values map as
`N_s: f_label → f_ring`:

    N_s=1: 0.200 → 0.400      N_s=3: 0.333 → 0.667
    N_s=2: 0.286 → 0.571      N_s=4: 0.364 → 0.727
                              N_s=6: 0.400 → 0.800

The manuscript (Sec. 3.2) therefore relocates the 1D onset to `f ≈ 0.67–0.80`
(per-ring) and compares it to the Paper-2 pair threshold at the collision
mid-window, `f*_F(t=T/2) ≈ 0.76`, **not** to the equal-effort point `f* ≈ 0.40`.
The two-paper comparison is like-for-like only in the per-ring convention.

## Rule going forward

- Report and cross-reference everything in the **per-ring** convention `f_ring`
  (the 2D convention). It is the physically-primary fraction and makes the 1D/2D
  and Paper-2/Paper-3 comparisons commensurate.
- Where a 1D figure keeps the per-bond axis for continuity with earlier drafts,
  state the conversion `f_ring = 2 f_label` in the caption (F3 relabels the 1D
  figures to per-ring so this footnote is no longer needed).

## Provenance

Definitions read directly from the code (F0, P3-F):
`run_chain.py:37`, `run_chain.py:101`; `run_square.py:40`, `run_square.py:100`.
No parameter was changed; this note only fixes the *reading* of `f`.
