# COMPTE RENDU — Paper 3, sprint P3-F (périodique + gel figures) pour reprise Claude

Repo `chiral-ring`, branche **`paper3-lattice`**. À lire avant de continuer.
Ce sprint a fait UNE grande chose non prévue par le prompt : rendre le moteur 2D
**périodique (tore)** sur décision utilisateur, ce qui a recalculé TOUS les
nombres 2D et modifié le manuscrit en conséquence. Le manuscrit `.tex` est
maintenant la **source de vérité** (les notes `results_lattice.md` /
`factorization_theorem.md` gardent des chiffres PRÉ-périodiques — à régénérer
depuis le `.tex` si on les réutilise).

## 0. État du repo

- Manuscrit : `paper/paper3_lattice.tex` (+ `paper3_lattice.pdf` fourni par
  l'utilisateur, compilé AILLEURS). **Pas de LaTeX sur cette machine**
  (`pdflatex` absent) → l'utilisateur recompile (`pdflatex ×2`).
- Commits P3-F (dans l'ordre) :
  - `351ddbd` F0 : note de convention f (`notes/paper3/f_convention.md`)
  - `e2b9197` moteur périodique + `torus_percolation` + wiring 5 scripts
  - `21d559e` F4 : manuscrit aux nombres périodiques + narratif exposants
  - `9acc8ca` F2 : test d'indépendance affûté (torus, barres bootstrap)
  - `b669629` F3 : 13 figures régénérées périodiques + captions
  - `69e5fa4` F5 : ce compte rendu
  - `4028b8d` F1/G7 : `verify_C1.py` + `exact_reduction.md` + déclaration G7
- Gates : **7/7 PASS** (G1-G6 `verify_lattice.py` sur `build_chain` resté OUVERT ;
  G7 `verify_C1.py` factorisation par bond + positivité d'incrément).

## 1. LA décision : bords périodiques (2D uniquement)

Le prompt F2 supposait des bords périodiques ; le moteur était en réalité
**OUVERT** (bug de cohérence : `build_square` n'ajoutait pas de liens
wraparound ; le manuscrit disait pourtant "periodic BC" ×3). Fork proposé à
l'utilisateur → choix : **rendre le moteur périodique**.

Implémenté :
- `build_square(..., periodic=True)` : **tore**, coordination 4 partout, liens
  wraparound `((x+1)%Lx, (y+1)%Ly)`, garde anti-double-lien pour Lx=2.
- `build_chain` reste OUVERT (le manuscrit ne revendique le périodique QUE pour
  le 2D ; les 6 gates utilisent `build_chain` → intacts).
- Nouveau `torus_percolation(mask)` dans `lattice.py` : clustering torique +
  critère de **wrapping** (BFS à coordonnées "déroulées" ; un cluster qui se
  reconnecte à son image périodique en x OU y percole). Retourne
  `(perc, Pinf, S)`. **Testé** unitairement (ligne/colonne pleine → wrap ;
  bloc 2×2 → non ; ligne trouée → non). Les 5 scripts 2D
  (`run_square/minority/factorization/universality/convergence`) l'utilisent
  (remplace le critère haut-bas ouvert, incohérent avec un partage périodique).

## 2. Verdicts F0–F4

- **F0** : gates 6/6 (résidus ≤ 2.9e-17). `f_convention.md` : 1D per-bond
  `Ns/(Np+2Ns)`, 2D per-ring `4Ns/(Np+4Ns)`, relation exacte 1D `f_ring=2·f_bond`.
- **F1 / G7** : **FAIT** (`4028b8d`). L'utilisateur a fourni `verify_C1.py` +
  `exact_reduction.md`. Script copié dans `src/lattice/` (adapté à la signature
  production `local_marginal_sigma(R,G,Tc,param,S_true,i)`), lancé, **PASS** :
  LEMME factorisation par bond vs moteur = **1.2e-15** ; C1 K=0→1 =
  `b·ln[(2+b)/(2-b)]` (sympy) ; grille K=1→2 min +4.4e-9 ; 20000 configs
  min +4.8e-7 ; séquences min +6.7e-5. **Prop. C = identité vérifiée.**
  Suite de gates maintenant **7/7** (G1-G6 `verify_lattice.py` + G7
  `verify_C1.py`). Manuscrit : ligne G7 dans le tableau, six→sept, `verify_C1`
  en annexe, résidu fait(1) synchronisé à 1.2e-15. `exact_reduction.md` (Prop.
  A/B prouvées, C vérifiée) dans `notes/paper3/`.
- **F2** : `run_independence_v2.py` (tore, `d_eff` min-image périodique,
  barres bootstrap). **PASS** :
  `r(1)=+0.029±0.012` (2.5σ, 6344 paires), `r(2)=-0.004±0.008`,
  `r(3)=-0.002±0.008` → prédiction `r(d≥2)=0` confirmée AVEC barres d'erreur.
  E3 mean|Δ|=0.022. Figure `fig_p3_independence.png` (copiée depuis `_v2`).
- **F4** : manuscrit mis à jour (voir §3). `verify_C1` numbers intacts.

## 3. Nombres périodiques (AVANT → APRÈS) — tous dans le `.tex`

| Quantité | Ouvert (avant) | **Périodique (après)** |
|---|---|---|
| f_c (crossing ρ=0.5) | 0.721 | **0.692** |
| conv(f) [f=.6,.67,.71,.75,.78,.80] | .004,.024,.111,.278,.465,.606 | **.005,.039,.199,.42,.66,.84** |
| std_ρ(conv) | 0.007 | **0.011** |
| r(1) | +0.112 (8σ, sans barre) | **+0.029±0.012** (2.5σ) |
| r(bulk d≥2) | −0.025 | **−0.004** (et r(2),r(3)≈0 avec barres) |
| ρ_c^align résidu | 0.003 / 0.001 | **0.003 / 0.001** (EXACT, inchangé) |
| ρ_c^perc résidu | 0.072 / 0.022 | **0.078 / 0.048** |
| E3 mean|Δ| | 0.040 | **0.022** (plancher 0.056) |
| ρ_c(f) table (.67→.80) | .51,.58,.70,.89,>.93 | **.52,.63,.87,>.93,>.93** |
| divergence ρ_c | f≈0.78 | **f≈0.77** (conv passe ½) |
| ρ_c^perc (intercept) | 0.407 | **≈0.43** fini-L (offset wrapping ; 0.407 = L→∞) |
| β/ν (convergence, robuste) | 0.27 (mauvais) | **0.107** (perc 0.104 — quasi EXACT) |
| γ/ν | 1.78 | **1.57–1.86** (perc 1.79 ; window-sensible) |
| ν | 2.04→drift 1.43 | **1.15 (slope), 1.27 (collapse), 1.49 (paire 16→24)** — bracket 4/3, non fixé |
| univ β/ν, γ/ν, ν (f=.5/.71) | 0.170/0.126, 1.70/2.00, 1.74/1.38 | **0.220/0.269, 1.77/2.04, 1.90/1.46** |
| no-backflip / M backflips | 35600/35600, 5600/6400 | **inchangés** |

**Lecture** : la loi de factorisation SURVIT et se renforce (ρ_c^align exact,
E3 0.022, β/ν quasi exact, ρ_c(f)=forme close via conv). Renversement narratif
exposants : β/ν devient le check robuste (était γ/ν). ν toujours non fixé
(→ HPC L~64-128). Universalité = **équivalence E3**, pas l'exposant.

## 4. F3 — figures (régénérées EN PLACE, mêmes noms, pas de suffixe _v2)

Décision : comme TOUTES les données ont changé (périodique), on régénère en
place (l'historique git garde les anciennes) → aucun re-câblage des
`\includegraphics`. 10 figures du manuscrit, toutes en données périodiques :
- `fig_p3_independence` : version F2 (barres + `d_eff` tore + E3).
- `fig_p3_convergence` : replot `replot_convergence.py` depuis
  `convergence_data.npz`, **titre neutre** (pas d'exposant en titre/légende).
- `fig_p3_chain_alignment`, `fig_p3_domains_chain` : **axe per-ring** (×2,
  Eq. 1), titres neutres (`run_chain.py` modifié).
- `fig_p3_2d_phase` : **overlay Eq. (7)** `(1-p_c)/(1-conv)` en pointillé blanc
  (`run_square.py` modifié) + titres neutres.
- `fig_p3_factorization`, `fig_p3_percolation`, `fig_p3_snapshots`,
  `fig_p3_prescriptionM`, `fig_p3_universality` : données périodiques, titres OK.
- Les 4 TODO(P3-F) du `.tex` sont retirés ; captions cohérentes.

## 5. Pièges rencontrés (P3-F)

1. **Bords ouverts vs périodiques** : `build_square` était ouvert (deg 2/3/4
   selon la position) alors que le manuscrit disait périodique. Vérifié
   empiriquement (24 vs 32 liens sur 4×4). → moteur rendu périodique.
2. **Percolation sur tore** : le critère haut-bas ouvert est incohérent sur un
   tore (bords adjacents). Il FAUT le critère de wrapping + labelling torique
   (sinon P_perc trivial ou clusters coupés à la couture).
3. **Offset de wrapping isotrope** : le crossing P=½ du wrapping "x OU y" est
   SOUS p_c à L fini (valeur universelle R_e(p_c)≈0.69) → ρ_c^perc mesuré ~0.43
   au lieu de 0.407 (forme close L→∞). C'est la source du résidu perc plus gros
   (0.048), pas un échec de la loi.
4. **Fenêtre de densité décentrée** : `run_convergence` avait `[0.34,0.46]`
   (ancien ρ_c 0.407) ; le crossing périodique 0.43 tombait au bord → slope
   biaisé → ν=0.88. Recentré `[0.37,0.49]` N_DENS=11 → ν=1.15 (paire 16→24 :
   1.49). TOUJOURS relancer avec la fenêtre centrée sur le vrai crossing.
5. **conv(f) monte** (tore = deg 4 partout, pas de dilution de bord) → tous les
   seuils bougent ; ρ_c(f) reste = forme close via la nouvelle conv.
6. **stdout Python bufferisé** vers fichier → pas de progression live ; les
   logs ne flushent qu'à la fin (les tracebacks, eux, sortent via stderr).
7. `RuntimeWarning nanstd empty slice` dans les métriques de collapse = inoffensif.

## 6. État du `.tex` et reste à faire

- `.tex` : nombres + narratif + captions + G7 à jour, 0 TODO. **NON compilé ici**
  (pas de LaTeX) → l'utilisateur lance `pdflatex paper3_lattice.tex` ×2, exiger
  0 erreur. Les figures sont dans `figures/paper3/` (le `\repofig` les inclut).
  **C'est le seul reste bloquant côté machine.**
- **F1/G7** : FAIT (§2) — 7/7 gates.
- **Notes stale** : `results_lattice.md`, `factorization_theorem.md`,
  `files/compte_rendu_paper3.md` contiennent des chiffres pré-périodiques.
  À régénérer depuis le `.tex` si réutilisés (le `.tex` fait foi).
- **ν sans équivoque** : HPC `run_convergence.py` (`LS=[16,24,32,48]`, R≥200) —
  fenêtre déjà recentrée [0.37,0.49].
- **NE PAS pousser sans accord utilisateur** (règle maintenue ; la branche EST
  déjà sur origin d'avant, mais ne pas pousser les nouveaux commits sans OK).

## 7. Commandes

    cd /Users/mostfa/chiral-ring/src/lattice
    ../../.venv/bin/python verify_lattice.py        # 6 gates (obligatoire)
    ../../.venv/bin/python run_independence_v2.py   # F2 r(d_eff) + barres (~6 min)
    ../../.venv/bin/python run_square.py            # phase(overlay)+percolation (~15 min)
    ../../.venv/bin/python run_minority.py          # rho_c(f) + FSS (~40 min)
    ../../.venv/bin/python run_factorization.py     # E1-E4 (~15 min)
    ../../.venv/bin/python run_universality.py      # nu(f), beta/gamma (~30 min)
    ../../.venv/bin/python run_convergence.py       # nu convergence (~15 min, fenêtre [0.37,0.49])
    ../../.venv/bin/python replot_convergence.py    # figure convergence titre neutre (instantané)
    ../../.venv/bin/python run_chain.py             # 1D per-ring (~3 min)

État : P3-F ≈ complet sauf F1/G7 (fichiers utilisateur manquants) et la
recompilation LaTeX (côté utilisateur). Moteur périodique, loi de
factorisation validée plus fort. Branche `paper3-lattice`.
