# COMPTE RENDU — P1-R (repo = source unique de vérité, rangement + commits)

> **Règle gravée :** tout `.tex` compilé est committé au moment où son PDF est
> produit ; les PDF sans source sont interdits ; les copies hors repo sont
> supprimées (par l'auteur, pas par l'assistant — l'assistant LISTE).

Branche : **`paper1-v1.1`** (créée depuis `main`). Pas de merge, pas de push
(attente OK explicite). Interpréteur : `.venv/bin/python`.

Décisions auteur pour cette session : **papier 3 mis de côté** (traité en session
paper3 dédiée) ; **articles 1 & 2 enchaînés**.

---

## R0 — Branche et état

- P1-G est **déjà mergé dans `main`** (commits `d56b6a4`, `8c0f2b8`, `3ef77fd`,
  `468c486`). Aucun cherry-pick nécessaire : `paper1-v1.1` part de `main`, les
  commits P1-G en sont ancêtres.
- `git status` propre avant travaux (seul `trilogie_consolidation/` untracked).

## R1 — Arborescence cible

Politique de synchronisation appliquée fichier par fichier (bundle = version
certifiée du 22 juillet) :

| Fichier | État repo avant | Action |
|---|---|---|
| `article2/*` (paper2.tex, pdf, zip, 14 figs) | présents, **identiques** (cmp -s) | consolidés à plat (`git mv` hors de `arxiv_submission_paper2/`) |
| `letter/main.tex` | différent (33340→33693 o) | **bundle remplace** (annotation H2 post-P1-G — attendu) |
| `extended/manuscript_extended.tex` | différent (41297→42017 o) | **bundle remplace** (certifié) |
| `extended/manuscript_extended.pdf` | différent | **bundle remplace** |
| `letter/main.pdf`, `letter/*.zip` | absents | **ajoutés** depuis bundle |

### Consolidation (règle : jamais deux copies du même manuscrit)
- Lettre : `old/arxiv_submission_paper1/` → **`letter/`** (supprimé après).
- Extended : `files/` → **`extended/`** (supprimé après).
- `old/extended/` (figs + tex antérieurs) : **supprimé** (doublon).
- Article 2 : figures + `paper2.tex` sortis de `arxiv_submission_paper2/` vers le
  plat ; sous-dossier retiré. Le zip arXiv reste (`arxiv_submission_paper2.zip`).

### Conservés (repo-uniques, non-doublons — **décision auteur** pour suppression)
- `paper/article1/old/manuscript.pdf` — compile lettre **pré-v1.1** (11 juil.).
- `paper/article1/old/manuscript_extended.pdf` — compile extended **pré-v1.1**.
- `paper/article1/old/manuscript_extended_fr.pdf` — **version française**.

  Ces trois PDF n'ont plus de source `.tex` correspondante active (sources
  remplacées par les v1.1). Techniquement « PDF sans source » — mais repo-uniques,
  donc **non supprimés** par l'assistant, listés ici pour ta décision.

### Dates figées (`\today` → `July 7, 2026`, label préservé)
```
letter/main.tex               \date{\small Draft: July 7, 2026}
extended/manuscript_extended  \date{\small Extended version: July 7, 2026}
article2/paper2.tex           \date{\small Draft: July 7, 2026}
article3/paper3_lattice.tex   \date{\small Draft: July 2026}   (déjà figée — hors scope session)
```
Aucun `\today` restant dans les tex versionnés d'articles 1 & 2.

### Vérification figures (obligatoire) — **PASS**
Chaque `\includegraphics` / macro `\figorbox` / `\repofig` résout à côté de son
tex : **letter 8/8 · extended 12/12 · article2 14/14**. Aucune absence.

## R2 — Fixes du script de production (`src/run_baseline.py`)
1. Histogrammes pré-transition : `[0.5, 1.5, 3.0, 8.0]` → `[0.1, 0.5, 1.5, 8.0]` ;
   sortie `fig3_histograms.png` → `fig3_histograms_v2.png`.
2. Légende co-transition hors graphe : `loc="center right"` →
   `loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=4, frameon=False` ;
   sortie `fig2_cotransition.png` → `fig2_cotransition_v2.png`.

Relance (`.venv/bin/python src/run_baseline.py`) :
```
model.py unit tests: ALL PASSED
validation I(S:F): max |analytic - numeric| = 5.764e-14
validation P(t=6.0): analytic 0.9894 vs numeric 0.9897 (diff 0.0003, sampling noise ~0.005)
```
Committé : script + `figures/paper1/fig2_cotransition_v2.png` +
`figures/paper1/fig3_histograms_v2.png`. Anciennes figures conservées
(provenance) ; `fig1_pip.png` re-rendu par effet de bord **restauré** (non commité).

## R3 — .gitignore
`.gitignore` racine complété : ajout `*.aux *.log *.out *.synctex.gz *.toc *.fls
*.fdb_latexmk`. Bloc existant conservé (`__pycache__/`, `*.pyc`, `.DS_Store`,
`.venv/`, `all files of session 1/`). **PDF des manuscrits explicitement
versionnés** (non ignorés).

## R4 — Gates post-rangement — **PASS (bloquant)**
```
.venv/bin/python src/verify_note.py        -> ALL CHECKS: PASS  (7/7, CHECK 7 Branch Arrow inclus)
.venv/bin/python src/multi/verify_multi.py -> ALL CHECKS: PASS
```
Le rangement n'a cassé aucun chemin. (Gates de la branche `paper3-lattice` :
hors scope — non touchés cette session.)

## R5 — Doublons hors repo (LISTÉS, rien supprimé)
```
find ~ -maxdepth 4 \( -name "main.tex" -o -name "manuscript_extended*" \
  -o -name "paper3_lattice*" -o -name "arxiv_submission*" \) \
  -not -path "*/chiral-ring/*" -not -path "*/Library/*"
→ (aucun résultat)
```
Aucune copie de manuscrit hors du repo.

**Dans le repo (untracked)** : `trilogie_consolidation/` — le bundle d'échange
lui-même (livraison du 22 juillet, articles 1/2/3). Non committé. À supprimer par
l'auteur une fois la synchro validée.

## R6 — Commits (branche `paper1-v1.1`, non poussée)
```
3efe42b  P1-R/R1  arborescence letter+extended, consolidation article2, dates figées
bb8622e  P1-R/R2  fixes run_baseline.py (histogrammes + légende) + figures _v2
fc0152e  P1-R/R3  .gitignore artefacts LaTeX
```
Papier 3 : **non traité** (décision auteur — session paper3 dédiée). Rappel du
point signalé : le `paper3_lattice.tex` du bundle **diffère** de la branche
`paper3-lattice` (bundle 51983 o vs 43617 o, ~295 lignes) — version strictifiée
biblio 32 / ligne G7. À arbitrer en session paper3.

## Arborescence finale (tracked)
```
paper/article1/letter/     main.tex  main.pdf  arxiv_submission_paper1_v1_1.zip  + 8 figs
paper/article1/extended/   manuscript_extended.tex  .pdf  + 12 figs
paper/article1/old/        manuscript.pdf  manuscript_extended.pdf  manuscript_extended_fr.pdf  (repo-uniques, décision auteur)
paper/article2/            paper2.tex  paper2_final.pdf  arxiv_submission_paper2.zip  + 14 figs
```

## En attente d'OK explicite
- `git push origin paper1-v1.1` (autorisé par R6 mais non exécuté sans ton feu vert).
- **PAS** de merge vers `main`.
- Suppression par l'auteur : bundle `trilogie_consolidation/` + les 3 PDF `old/`.
