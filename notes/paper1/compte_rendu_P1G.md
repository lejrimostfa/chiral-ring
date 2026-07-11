# COMPTE RENDU — P1-G (Branch Arrow Theorem + intégration v1.1 papier 1)

Repo `chiral-ring`. Session **P1-G**. Branche courante : **`paper3-lattice`**
(NON poussée). Aucun push effectué (attente OK explicite).

Interpréteur : `/Users/mostfa/chiral-ring/.venv/bin/python` (Python 3.14.4,
numpy 3.14). `python` nu n'est **pas** sur le PATH.

---

## V0 — Sanity gates existants : **PASS**

- `.venv/bin/python src/verify_note.py`        → `ALL CHECKS: PASS` (6/6 avant P1-G).
- `.venv/bin/python src/multi/verify_multi.py` → `ALL CHECKS: PASS`.

Les deux verts avant toute modification. Rien n'a été committé avant ce point.

---

## V1 — Gate Branch Arrow Theorem : **PASS** (avec réserve fp cosmétique)

`src/verify_branch_arrow.py` créé (copie **verbatim** de la spec P1-G ; non
modifié). Sortie déterministe (`rng(1)`), **verbatim cette machine** :

```
(a) max |sigma[xi]| over 200 correlated instances : 0.00e+00
(b) max |<sigma_s> - D_KL(joint)|                 : 0.00e+00
    additivity gap |KL_joint - sum KL_n| (MAX over trials, correlated): 1.237 -> additive form needs cond. independence
    product law: |KL_joint - sum_n eps ln[(1+eps)/(1-eps)]| = 8.88e-16
(c) T-even control: <sigma_s> = 0.0 identically while D = 0.432 > 0
```

**Écart unique vs spec** : ligne « product law » = `8.88e-16` (spec P1-G :
`0.00e+00`). Bruit machine — la somme est analytiquement **exactement 0** ;
l'ordre de sommation des 32 termes (n=5) sur numpy 3.14 laisse un résidu
~1e-16. Les 4 autres lignes, **y compris les valeurs sensibles à la graine**
(`gap=1.237`, `D=0.432`), sont **bit-exactes** vs spec. Le script est la spec :
il n'a **pas** été « corrigé » pour forcer `0.00e+00`. Verdict : gate
substantiellement PASS ; l'écart est présentationnel (arrondi `%.2e` d'un zéro
machine), pas une violation du théorème.

**CHECK 7** ajouté dans `src/verify_note.py` : lance `verify_branch_arrow.py`
en **sous-process**, parse ses verdicts numériquement (tolérances : (a),(b) et
résidu produit `< 1e-9` ; `|gap-1.237| < 1e-3` ; T-even `< 1e-15` ; `D > 0`),
replie le résultat dans `ALL CHECKS`. Le parse tolérant absorbe le bruit fp
tout en asserant les conclusions du théorème. Le script **reste invocable
seul**. Suite complète re-lancée :

```
CHECK 7: Branch Arrow Theorem gate (src/verify_branch_arrow.py)
  (a) sigma[xi]=0 : True;  (b) <sigma_s>=D_KL : True;  product-law additivity exact : True
  correlated gap = 1.237 (>0, needs cond. indep.);  T-even <sigma_s>=0.0 while D=0.432>0 ... PASS: True
============================================================
ALL CHECKS: PASS
```

→ **7/7 PASS**. Commit `b03c0a1` (V1).

---

## V2 — Intégration des fichiers v1.1 : **PASS** (source ; PDF/zip à régénérer)

> Note historique : à l'ouverture de session ces fichiers étaient **absents**
> (disque + historique git toutes branches). L'auteur a ensuite ajouté les
> dossiers `paper/article1` & `paper/article2`, puis **remplacé `main.tex` et
> ajouté `fig2_sequence.png`** par la v1.1. État final vérifié ci-dessous.

Dossier cible confirmé : **`paper/article1/arxiv_submission_paper1/`**.

- `main.tex` = **v1.1** (vérifié) : théorème dans l'abstract (« …instances of a
  model-independent theorem (Methods)… » + « The laws carry no arrow; conditioned
  histories do. »), paragraphe **Branch Arrow Theorem** en Methods (H1-H3 +
  preuve 4 lignes + « scripts in the repository »), paragraphe decoherent
  histories avec refs `GellMann1993` (l.164, PRD 47, 3345 (1993)) et `Hartle2013`
  (l.165, arXiv:1301.2844).
- `fig2_sequence.png` : **présent**, référencé l.109.
- **Tous les 8 `\includegraphics` résolvent sur disque** à côté du tex :
  fig1_schematic, fig4, fig5b, fig3, fig6, fig8, fig9, fig2_sequence → OK.
- **Non recompilé** (pas de LaTeX ici) : `manuscript.pdf` est la compile
  **pré-v1.1** (7 juil. 16:28, antérieure aux éditions v1.1) ; l'ancien zip
  arXiv v1.0 a été **retiré** par l'auteur. PDF + zip → **à régénérer après
  validation auteur** (cf. RESTE À FAIRE).

Commit `d7c04f8` (V2) : ajoute `paper/article1` (lettre v1.1 + PDFs) et
`paper/article2` (artefacts fournis). Extended papier 1 : **intouchée**.

---

## V3 — Traçabilité notes : **PASS**

`notes/paper1/changelog_v1.1.md` créé : consigne la substance de v1.1
(théorème H1-H3 + preuve 4 lignes Methods ; paragraphe modèle-indépendant
essentiel/inessentiel ; `fig2_sequence` = séquence 4 étapes + boucle de
révocabilité ; phrase abstract ; Discussion decoherent histories + refs
**GellMann1993** & **Hartle2013**, vérifiées en session — consigné, à ne pas
re-vérifier), le gate CHECK 7, et l'échec T-even devenu corollaire. Le
document **marque explicitement** l'absence des fichiers manuscrit (V2 en
attente) et n'affirme pas leur intégration. Commit `1ace36c` (V3).

---

## RESTE À FAIRE (explicite)

1. ~~[BLOQUANT V2] Livrer main.tex v1.1 + fig2_sequence.png + dossier cible~~ →
   **FAIT** (auteur, en session) ; V2 vérifié PASS (source). ~~Vérifier que les
   figures résolvent~~ → **FAIT** (8/8).
2. **[LIVE] Recompiler la lettre v1.1** (hors de cette machine, pas de LaTeX
   ici) → produire le `manuscript.pdf` v1.1 (le PDF sur disque est la compile
   pré-v1.1 du 7 juil.) **et regénérer le zip arXiv** `arxiv_submission_paper1`
   (l'ancien v1.0 a été retiré ; son `main.tex` bundlé était pré-v1.1 et sans
   fig2_sequence), après validation auteur.
3. **Miroir v1.1 sur la version extended** du papier 1 — session dédiée, tex à
   fournir. Hors périmètre P1-G. (Les PDFs extended sur disque —
   `manuscript_extended{,_fr}.pdf` — sont committés tels quels, non modifiés.)
4. **Décision de branche** : les commits P1-G (V1-V4) sont sur `paper3-lattice`
   (HEAD courant). Le papier 1 « complet » vit sur `main`. À trancher : rebaser
   / cherry-pick les commits P1-G vers `main` ou une branche `paper1-v1.1`
   dédiée avant push.
5. **Push** : uniquement sur **OK explicite** de l'auteur. Rien n'est poussé.

---

## Commits de la session (branche `paper3-lattice`, non poussés)

- `b03c0a1` — V1 : gate `src/verify_branch_arrow.py` + CHECK 7 (7/7 PASS).
- `1ace36c` — V3 : `notes/paper1/changelog_v1.1.md` (traçabilité v1.1).
- `8bb69ed` — V4 : compte rendu P1-G (version initiale, V2 alors bloqué).
- `d7c04f8` — V2 : intégration lettre v1.1 (`paper/article1`) + `paper/article2`.
- *(ce commit)* — V3/V4 : notes remises à jour après livraison auteur (V2 PASS).
