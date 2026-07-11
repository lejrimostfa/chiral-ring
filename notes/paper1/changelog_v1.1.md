# Paper 1 (lettre) — changelog v1.1 + Branch Arrow Theorem

Session **P1-G**. Repo `chiral-ring`. Traçabilité de la révision v1.1 du papier 1
(version lettre) et du gate modèle-indépendant qui la certifie.

> ✅ **État des fichiers manuscrit v1.1** (mis à jour après livraison auteur) :
> dossier `paper/article1/arxiv_submission_paper1/` — `main.tex` **v1.1**
> intégré (Branch Arrow Theorem dans abstract/Main/Methods ; abstract clos par
> « The laws carry no arrow; conditioned histories do. » ; paragraphe
> decoherent histories + refs `GellMann1993` & `Hartle2013`), `fig2_sequence.png`
> présent, et les 8 `\includegraphics` (fig1_schematic, fig3, fig4, fig5b, fig6,
> fig8, fig9, fig2_sequence) **résolvent tous sur disque**. Le Methods promet
> « scripts in the repository » → tenu par CHECK 7 / `src/verify_branch_arrow.py`.
> **Restent à régénérer** (pas de LaTeX ici, non recompilé) : `manuscript.pdf`
> (compile pré-v1.1 du 7 juil.) et le zip arXiv (l'ancien v1.0 a été retiré) —
> après validation auteur. Version **extended** du papier 1 : intouchée
> (miroir v1.1 en session dédiée).

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

## 4. `fig2_sequence.png` (figure conceptuelle, **intégrée**)

Séquence en **4 étapes**, chaque étape certifiée par une quantité *mesurée* du
papier ; boucle de **révocabilité** dessinée sous les étages cohérents.
Présente à côté du `main.tex` lettre (`paper/article1/arxiv_submission_paper1/`),
référencée à la ligne `\includegraphics{fig2_sequence.png}` du tex.

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
