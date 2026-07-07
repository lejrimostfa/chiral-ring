# chiral-ring

**An exactly solvable model of the arrow of time.**

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)
![Verification: machine precision](https://img.shields.io/badge/verification-machine%20precision-brightgreen.svg)

A three-site chiral ring — a time-reversal-odd pointer doublet — monitored by a
bath of recording qubits. The dynamics is pure dephasing and conserves the ring
chirality, so the whole model reduces to a two-branch (or, for several rings, a
`2^M`-branch) structure with closed forms for every observable: coherences,
record distinguishability, entropy production, redundancy, concurrence.

Nothing here is a fit or an approximation. **Every analytic claim is checked
against brute-force state-vector evolution and against an independently recoded
verification suite, at machine precision.**

---

## The two papers

This repository accompanies two papers.

- **Paper I** — *Fixation of a time-reversal-odd observable by decoherence:
  record distinguishability as the local arrow of time in an exactly solvable
  model* (manuscript under arXiv submission).
  Single ring. Code: [`src/`](src/) · verification:
  [`src/verify_note.py`](src/verify_note.py) · figures:
  [`figures/paper1/`](figures/paper1/) · notes: [`notes/paper1/`](notes/paper1/).

- **Paper II** — *Collective irreversibility and arrow competition in shared
  quantum environments* (in preparation).
  Multi-ring extension. Code: [`src/multi/`](src/multi/) · verification:
  [`src/multi/verify_multi.py`](src/multi/verify_multi.py) · figures:
  [`figures/paper2/`](figures/paper2/) · notes:
  [`notes/paper2/`](notes/paper2/).

Paper II builds on Paper I: `src/multi/` imports the Paper-I closed forms and
its unit tests confirm exact reduction to Paper I at `M = 1`.

---

## Key results

**Paper I.** Under a branch-relative definition of entropy production, the
irreversibility and the record content are the *same quantity*:
`⟨σ⟩(t) = D(t)` exactly (accumulated Kullback–Leibler record distinguishability).
Polarization, redundancy and entropy production co-transition on one master
curve; the pooled decoherence-to-polarization timescale ratio is
`t_P / t_S = 0.953 ± 0.005`, and the data collapse holds at RMS residual
`0.005`.

**Paper II** (verdicts, with the number that decides each):

| Conjecture | Verdict | Evidence |
|---|---|---|
| **C1** total current objectifies, individuals don't (DFS) | confirmed | `P_indiv = 0.50` (chance) vs `P_total = 1.0` |
| **C2** additivity deficit = a multi-information | reinterpreted | deficit is a **theorem**: `Δ = I(s₁;s₂\|ξ) + C_bwd`, verified to `1e-15`; extensive `Δ = 0.688·N` (`R² = 0.99995`); `I/Δ → 1/8` proven at leading order |
| **C3** entanglement lives off the records | confirmed | concurrence `C(t)` monotone in relative distinguishability; `C ≡ 1` in the decoherence-free subspace |
| **C4** detailed fluctuation theorem | confirmed | `ln[p(σ)/p(−σ)]` slope `= 1.00 ± 0.02` |
| **C5** arrow frustration under bath sharing | resolved | one arrow reverses at a disorder-stable threshold; at equal effort the F/M prescriptions are symmetric duals (`f* ≈ 0.40`) |

The Paper-I data collapse **survives bath sharing** (RMS `0.018`) and breaks
only in the decoherence-free subspace (RMS `0.40`) — the break is itself the
signature of where records fail to form.

Full verdicts and figures: [`notes/paper2/results_multi.md`](notes/paper2/results_multi.md).

---

## Quick start

```bash
pip install numpy scipy matplotlib sympy

python src/verify_note.py          # Paper I: re-derives every analytic claim (~1 min)
python src/multi/verify_multi.py   # Paper II: independent from-scratch suite
```

Both print `ALL CHECKS: PASS`. `verify_note.py` runs six exact checks
(time-reversal operator, forward/backward record laws, `σ = 0` identity,
integral fluctuation theorem, Bayes posterior); `verify_multi.py` verifies the
deficit decomposition theorem, the `1/8` constant (symbolically, via `sympy`),
and the sign survey — none of it importing the production code.

---

## Reproducing the figures

Every figure regenerates from a fixed seed. Output lands in
`figures/paper1/` and `figures/paper2/` automatically.

**Paper I** (from the repo root):

| Figures | Command |
|---|---|
| Figs. 1–3 (partial-info plots, co-transition, histograms) | `python src/run_baseline.py` |
| Figs. 4–6 (collision protocol, finite-size scaling, collapse) | `python src/run_next.py` |
| Figs. 7–9 (variant B, hysteresis, ratchet) | `python src/run_variantB.py` |
| Error bars / disorder statistics (Fig. 5b, ~20 min) | `python src/errorbars.py` |

**Paper II** (from the repo root):

| Figures | Command |
|---|---|
| Scenarios (a)–(f): DFS, deficit, frustration, fluctuation theorem | `python src/multi/run_scenarios.py all` |
| Phase diagram `f*(t)` (the paper's figure) | `python src/multi/run_phase_diagram.py` |
| Exact deficit decomposition | `python src/multi/analyze_c2.py` |
| Deficit scaling, collapse error bars, `f*` stability | `run_deficit_scaling.py`, `run_collapse_errors.py`, `run_fstar_stability.py` |

---

## Repository layout

```
src/            Paper I — model, exact closed forms, trajectory protocols
  verify_note.py    six machine-precision checks of the analytic claims
src/multi/      Paper II — multi-ring extension
  model_multi.py    closed forms (imports Paper I; M=1 reduction is a unit test)
  brute_multi.py    brute-force state-vector validation (N=8)
  verify_multi.py   independent from-scratch verification suite
figures/
  paper1/  paper2/   generated figures, one folder per paper
notes/
  paper1/           model spec, backward-process derivation, manuscripts
  paper2/           multi-ring definitions + theorem, results & verdicts
```

Validation gates: exact reduction to Paper I at `M = 1`; brute-force agreement
at `N = 8`; exact additivity (`Δ = 0`) for disjoint baths — all to machine
precision.

---

## Citation

```bibtex
@unpublished{lejri_chiralring_paper1,
  author = {Lejri, Mostfa},
  title  = {Fixation of a time-reversal-odd observable by decoherence:
            record distinguishability as the local arrow of time in an
            exactly solvable model},
  note   = {Manuscript under arXiv submission},
  year   = {2026}
}

@unpublished{lejri_chiralring_paper2,
  author = {Lejri, Mostfa},
  title  = {Collective irreversibility and arrow competition in shared
            quantum environments},
  note   = {In preparation},
  year   = {2026}
}
```

---

## Author

**Mostfa Lejri** — Independent researcher — ORCID
[0000-0002-0283-6069](https://orcid.org/0000-0002-0283-6069).

*Both papers made extensive use of Anthropic's Claude as a research assistant;
every analytic claim is checked by an independently recoded verification suite.
See the manuscripts' acknowledgments.*

Released under the [MIT License](LICENSE).
