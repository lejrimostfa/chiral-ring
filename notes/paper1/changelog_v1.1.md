# Paper 1 (lettre) — changelog v1.1 + Branch Arrow Theorem

Session **P1-G**. Repo `chiral-ring`. Traçabilité de la révision v1.1 du papier 1
(version lettre) et du gate modèle-indépendant qui la certifie.

> ⚠️ **État des fichiers manuscrit v1.1** : `main.tex` (v1.1), `fig2_sequence.png`
> et `main.pdf` (compilé côté Claude, 9 pages) étaient annoncés comme *fournis*
> mais **ne sont présents ni sur le disque ni dans l'historique git** (toutes
> branches). Aucun dossier d'artefacts « lettre » du papier 1 n'existe dans ce
> repo (seul `paper/paper3_lattice.{tex,pdf}` est présent ; `figures/paper1/`
> contient des PNG à nomenclature différente — `fig1_pip`, `fig2_cotransition`,
> …, **pas** `fig1_schematic`/`fig3`/`fig5b`…). **V2 est donc en attente de
> livraison des fichiers** (voir RESTE À FAIRE dans le compte rendu P1-G). Ce
> présent document consigne la *substance* de v1.1 et le gate déjà livré ; il
> n'affirme pas que le tex/pdf/figure sont intégrés.

## 1. Branch Arrow Theorem (cœur de v1.1)

Théorème **modèle-indépendant** énoncé dans l'abstract, le Main et les Methods.

**Hypothèses** (aucune ne mentionne l'anneau, la géométrie, N, ni les couplages) :
- **(H1)** deux branches `s = ±1` échangées par l'opérateur antiunitaire `Theta`.
- **(H2)** covariance antipodale des records `P(xi|-s) = P(-xi|s)` **et** protocole
  arrière Theta-imagé `P_bwd(xi~|s) = P_fwd(xi|-s)`, avec `xi~ = reverse(-xi)`.
- **(H3)** prior symétrique sur les branches.

**Conclusions** (preuve 4 lignes en Methods) :
- **(a)** `sigma[xi] = 0` pour *chaque* record `xi` — les lois ne portent aucune
  flèche.
- **(b)** `<sigma_s> = D_KL( P(.|s) || P(.|-s) ) ≥ 0` — la flèche vit dans les
  histoires *conditionnées* ; **additive sur les records SSI ceux-ci sont
  conditionnellement indépendants**.
- **(c)** contrôle **T-pair** (corollaire) : si la conjugaison *fixe* la branche au
  lieu de l'échanger, alors `sigma_s = 0` identiquement **alors même que** la
  distinguabilité marginale `D` des deux branches croît. L'ancien « échec T-even »
  devient ce corollaire.

## 2. Gate : `src/verify_branch_arrow.py` (= CHECK 7)

Script autonome (numpy seul, `rng(1)` figé). Vérifie (a)-(c) sur des lois de
records **aléatoires CORRÉLÉES** satisfaisant (H1)-(H3), **sans aucun input du
modèle en anneau** — le théorème est modèle-indépendant, son gate aussi. Tient
la promesse « scripts in the repository » des Methods.

Sorties déterministes (verbatim, cette machine — numpy 3.14 / Python 3.14) :

```
(a) max |sigma[xi]| over 200 correlated instances : 0.00e+00
(b) max |<sigma_s> - D_KL(joint)|                 : 0.00e+00
    additivity gap |KL_joint - sum KL_n| (MAX over trials, correlated): 1.237 -> additive form needs cond. independence
    product law: |KL_joint - sum_n eps ln[(1+eps)/(1-eps)]| = 8.88e-16
(c) T-even control: <sigma_s> = 0.0 identically while D = 0.432 > 0
```

**Écart vs spec** : la seule ligne « product law » affiche `8.88e-16` là où la
spec P1-G prédit `0.00e+00`. C'est du **bruit machine** (la somme est
analytiquement exactement nulle ; l'ordre de sommation des 32 termes sur numpy
3.14 laisse un résidu ~1e-16). Les 4 autres lignes — dont les valeurs sensibles
à la graine `gap=1.237`, `D=0.432` — sont **bit-exactes** vs spec. Le script
n'a **pas** été modifié (copie verbatim de la spec).

**CHECK 7** dans `src/verify_note.py` : lance le script en sous-process, parse
ses verdicts numériquement (tolérances : (a),(b),résidu produit `< 1e-9` ;
`gap ≈ 1.237` ; T-even `< 1e-15` ; `D > 0`) et replie le résultat dans
`ALL CHECKS: PASS`. Suite désormais **7/7 PASS**. Le script reste invocable seul.

## 3. Paragraphe « The statement is model-independent » (Main)

- **Essentiel** = { pointeur T-odd binaire ; covariance antipodale des records ;
  condition initiale record-free }.
- **Inessentiel** = { les 3 sites ; la géométrie ; les couplages ; N }.

## 4. `fig2_sequence.png` (figure conceptuelle, **en attente de livraison**)

Séquence en **4 étapes**, chaque étape certifiée par une quantité *mesurée* du
papier ; boucle de **révocabilité** dessinée sous les étages cohérents.
À déposer à côté du `main.tex` lettre (même dossier que les autres `fig*` de la
lettre), une fois les fichiers fournis.

## 5. Abstract & Discussion

- **Abstract** : phrase du théorème + « The laws carry no arrow; conditioned
  histories do. »
- **Discussion** : paragraphe *decoherent histories* + 2 références
  **Gell-Mann–Hartle, PRD 47, 3345 (1993)** et **Hartle, arXiv:1301.2844**.
  Les deux ont été **vérifiées en session contre les pages APS/arXiv** ;
  **ne pas re-vérifier, ne pas modifier**.

## 6. Hors périmètre (P1-G)

- Version **extended** du papier 1 : recevra son propre miroir de v1.1 dans une
  session dédiée (tex à fournir). **Ne pas toucher** ici.
