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

## V2 — Intégration des fichiers v1.1 : **BLOQUÉ (fichiers jamais livrés)**

Les fichiers annoncés comme *fournis* — `main.tex` (v1.1), `fig2_sequence.png`,
`main.pdf` — sont **introuvables** :

- absents du disque (racine, `files/`, `paper/`, `figures/`, `~/Downloads`,
  `~/Desktop`, recherche système) ;
- **jamais** présents dans l'historique git (toutes branches : `main`,
  `paper3-lattice`, remotes) — vérifié via `git log --all --name-only`.

De plus, **aucun dossier d'artefacts « lettre » du papier 1 n'existe** dans ce
repo :

- seul tex/pdf présent : `paper/paper3_lattice.{tex,pdf}` (papier 3) ;
- `figures/paper1/` contient bien 11 PNG mais à **nomenclature différente**
  (`fig1_pip`, `fig2_cotransition`, `fig3_histograms`, `fig4_collision`,
  `fig5_scaling`, `fig5b_ratio_errorbars`, `fig6_collapse`, `fig7_variantB`,
  `fig8_hysteresis`, `fig9_ratchet`) — **pas** `fig1_schematic`, ni les
  `fig3/fig4/fig5b/fig6/fig8/fig9` référencés par le `main.tex` v1.1 attendu.

**Conséquence** : impossible de placer les fichiers, impossible de vérifier
que les `\includegraphics` de `main.tex` (fig1_schematic, fig3, fig4, fig5b,
fig6, fig8, fig9, fig2_sequence) résolvent sur disque. **Aucune fabrication**
de manuscrit / PDF / figure n'a été faite (ce serait inventer l'artefact).
Recompilation : de toute façon impossible ici (pas de LaTeX). V2 **en attente
de livraison** des 3 fichiers + confirmation du dossier cible.

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

1. **[BLOQUANT V2]** Livrer sur cette machine `main.tex` (v1.1),
   `fig2_sequence.png`, `main.pdf`, **et** confirmer/creer le dossier cible des
   artefacts « lettre » du papier 1 (il n'existe pas encore). Sans ces
   fichiers, V2 ne peut être ni exécuté ni vérifié.
2. Une fois V2 fait : vérifier que **toutes** les figures référencées par
   `\includegraphics` existent à côté du tex ; **ne pas recompiler** (pas de
   LaTeX ici ; le PDF fourni fait foi).
3. **Miroir v1.1 sur la version extended** du papier 1 — session dédiée, tex à
   fournir. Hors périmètre P1-G.
4. **Régénération des zips arXiv** après validation auteur.
5. **Décision de branche** : les commits P1-G (V1/V3) sont sur `paper3-lattice`
   (HEAD courant). Le papier 1 « complet » vit sur `main`. À trancher : rebaser
   / cherry-pick les commits P1-G vers `main` ou une branche `paper1-v1.1`
   dédiée avant push.
6. **Push** : uniquement sur **OK explicite** de l'auteur. Rien n'est poussé.

---

## Commits de la session (branche `paper3-lattice`, non poussés)

- `b03c0a1` — V1 : gate `src/verify_branch_arrow.py` + CHECK 7 (7/7 PASS).
- `1ace36c` — V3 : `notes/paper1/changelog_v1.1.md` (traçabilité v1.1).
- *(ce fichier)* — V4 : compte rendu P1-G.

V2 : **aucun commit** (bloqué, fichiers absents).
