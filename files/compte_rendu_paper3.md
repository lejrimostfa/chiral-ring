# COMPTE RENDU — Paper 3 (lattice) pour reprise Claude Code

Repo `chiral-ring`. Branche **`paper3-lattice`** (NON poussée). État à la fin de
la session P3-A..P3-D. À lire avant de continuer ; ne pas redécouvrir ce qui suit.

## 1. État du repo

- `main` (poussé sur github.com/lejrimostfa/chiral-ring) = Papers 1 & 2 complets
  + README + LICENSE. Branche `paper2-multiring` supprimée (fusionnée).
- `paper3-lattice` (local seulement) — commits :
  - `c9828fd` prototype `files/` (ajouté par l'utilisateur)
  - `78fd013` P3-A moteur + collision + 6 gates
  - `e94148e` P3-B chaîne 1D
  - `d1c9054` P3-C percolation 2D
  - `919884e` P3-D minorité + FSS + rapport
- Nouveau code : `src/lattice/` (lattice.py, verify_lattice.py, run_chain.py,
  run_square.py, run_minority.py). Notes : `notes/paper3/` (arrow_definition.md,
  results_lattice.md). Figures : `figures/paper3/` (8 png).
- Env : `.venv` (numpy 2.4, scipy 1.17, matplotlib 3.11, sympy 1.14). Lancer les
  scripts avec `/Users/mostfa/chiral-ring/.venv/bin/python`.
- Papers 1 & 2 (`src/`, `src/multi/`) : INTOUCHABLES (fondation validée).

## 2. Validé — NE PAS redécouvrir (verify_lattice.py, 6/6 PASS)

- G1 `Delta=I+C_bwd` sur la loi réseau (théorème papier 2) — 2.8e-17
- G2 `f=0` -> Th.2 papier 1 (`<sigma_i>=D_i`) — 2.2e-16
- G3 force brute L=3 (state vector) vs branches — 2.9e-17
- G4 **marginal local == global** (signe 100%, corr 0.97) -> le marginal local
  EST la définition de production, coût `O(2^(deg+1))/site`, jamais `2^L`.
  C'est ce qui rend L=16 (256 anneaux) traitable.
- G5 collision : additivité exacte (deficit=0), `<sigma_i>=D_i` — 0.0
- G6 collision L=2 reproduit la limite C5 papier 2 à f=0 (`+D`,`-D`) — 2.2e-16

## 3. Conventions & décisions (fixées, à respecter)

- **Flèche** : `tau_i = sign(<sigma_i^local>)` = direction d'ACCUMULATION des
  records accessibles (PAS le signe du ratio d'inférence). Réconcilié avec C5
  papier 2 dans `notes/paper3/arrow_definition.md` (coïncident à f=0, diffèrent
  dans le transitoire).
- **Protocole collision** (production) : fenêtre fixe `dt`, biais
  `beta_n(S)=sin(2 a_n(S)(dt-Tc_n))`, `a_n(S)=sum_i S_i G[i,n]`. Frais `Tc=0`
  tiré forward `(1+r beta)/2` ; chargé (contrarian) `Tc=T` tiré backward
  `(1-r beta)/2`. Records permanents/monotones -> tue les récurrences.
- **Balayer f par Ns** (structurel), JAMAIS par la force de couplage
  (récurrences). `f=4Ns/(Np+4Ns)` en 2D (bulk deg 4), `Ns/(Np+2Ns)` en 1D.
- **Params** : Np=8 en 2D (résolution f fine), Np=3 en 1D. dt=0.3, T_charge=1.0,
  g_scale=0.25. Branche vraie `S_true = tous +1`.
- **Contrarians** : fraction `rho` d'anneaux, bain PRIVÉ chargé (Tc=T), liens
  partagés FRAIS.

## 4. Résultats acquis (chiffres)

- Baseline continue exacte : `A(f=0)=1-2*round(rho*L)/L` (dev 0.0).
- 1D : contagion seuil `f~0.36-0.40` (= f* de la paire papier 2) ; fraction
  convertie UNIVERSELLE en f (indép. de rho) ; **crossover** (domaines L=8/16/32
  se superposent, pas de divergence -> pas de transition en 1D).
- 2D : **vraie transition de percolation**, `f_c=0.721` (rho=0.5) ; seuil d'occup.
  forward `p_c=0.5927`. Phase map : frontière `rho_c(f)` monte avec f.
- `rho_c(f)` (renversement flèche GLOBALE, A=0) : ~0.50 (f<=0.6) -> 0.51,0.58,
  0.70,0.89 (f=0.67..0.78) -> **DIVERGE (>0.93)** au-delà de f~0.78. « Témoins
  partagés votent double » : frais=forward -> rho_c>1/2. SIGNE dépend de la
  prescription (frais vs chargés) = dualité F/M du papier 2.
- **Universalité : nu=4/3 (percolation 2D)** — FSS pilotée par rho (continu) sur
  L={8,12,16}, minimum INTÉRIEUR net à nu=1.33 ; rho_c^perc=0.391 à f=0.5.

## 5. Pièges rencontrés (vrais gotchas de la session)

1. Baseline : comparer à `1-2*round(rho*L)/L`, PAS `1-2*rho` (discrétisation du
   nombre de contrarians ; sinon fausse « déviation » ~0.025).
2. Percolation visible seulement pour `rho > ~0.41` : sinon forward
   (fraction 1-rho) est déjà > p_c=0.593 à f=0, ça percole avant tout partage.
3. Ns entier -> grille f GROSSIÈRE -> **impossible de fixer nu par une FSS
   pilotée par f**. Solution : piloter par **rho continu** + **L plus grand**.
4. Métrique de collapse naïve RAIL vers grand nu (L^(1/nu)->1 collapse tout
   trivialement). Restreindre à la fenêtre de transition + plage nu physique.
5. Longueur de corrélation intégrée fragile en 1D -> utiliser la longueur de
   domaine (échelle de corrélation robuste du champ de signes).
6. `rho_c` peut sortir de la plage testée -> nan ; étendre la plage (jusqu'à
   0.93) et marquer la divergence.
7. Coût : L=16 ~2s/réalisation (n_traj=2500). Goulot = boucle Python sur les
   anneaux (local_marginal_sigma par site). n_traj=1500-2500 suffit pour le SIGNE.
8. cwd du shell : lancer avec chemins ABSOLUS (le cwd dérive entre calls).

## 6. MISE À JOUR — P3-E (factorisation) + exposants (fait après P3-D)

**P3-E — Loi de factorisation (`8b36d64`, `factorization_theorem.md` = Section 4).**
Le champ de flèches sous prescription F EST une percolation de sites
indépendants d'occupation `p_eff = 1 - rho(1-conv(f))`. Prouvé (run_factorization.py, L=8) :
- E1 indépendance : bulk `r(4<=d<=9)=-0.025`~0 ; seule correction voisine
  `r(1)=+0.112` (lien partagé). No-backflip 35600/35600. conv(f) universelle en
  rho (std 0.007).
- E2 formes closes SANS fit : `rho_c^align=1/[2(1-conv)]` résidu **0.003 (exact)** ;
  `rho_c^perc=(1-p_c)/(1-conv)` résidu max 0.072.
- E3 test ultime : champ de sites indépendants reproduit P_perc, mean|Δ|=0.040 <
  plancher de bruit 0.056 → **indistinguable**.
- E4 frontière : sous prescription **M** (liens partagés chargés), la
  factorisation CASSE (5600/6400 anneaux frais retournés, A→-1). Signe fixé par
  la prépa des témoins = dualité F/M du papier 2.

**Exposants — statut HONNÊTE (`0eecea6`, `bbacaca`).** À L<=24 ils NE sont PAS
pinnés proprement :
- `run_universality.py` (métrique de collapse, fragile) : nu=1.74 (f=0.5) vs
  1.38 (f=0.71) → scatter, PAS f-indépendant. Le nu=4/3 "propre" de P3-D était
  en partie une grille chanceuse.
- `run_convergence.py` (estimateur robuste = pente au seuil ; densité exacte
  `n_c/L^2` ; L={12,16,24}) : `rho_c`=0.407 stable ; **`gamma/nu=1.78` == perc
  1.79 (quasi exact, robuste)** ; nu : la pente se courbe vers la valeur perc
  aux grandes tailles (paire 16->24 -> nu_eff=1.43->4/3) mais bruit R=50 empêche
  de clouer ; beta/nu=0.27 (le plus dur). Données brutes : `convergence_data.npz`.
- **Conclusion** : cohérent avec percolation 2D (gamma/nu exact, nu qui converge),
  mais la mesure DIRECTE ne le prouve pas à ces tailles. **La preuve rigoureuse
  d'universalité = E3 (équivalence), PAS l'exposant en force.** Un ν sans
  équivoque en direct demande L~64-128 (HPC) car chaque site coûte un
  échantillonnage quantique.

## 7. Reste à faire / ouvertures

- **Rédaction du manuscrit Paper 3** (pas commencée). Squelette = `results_lattice.md`
  + `factorization_theorem.md` (Section 4). Angle à vendre : "au lieu de forcer
  l'exposant sur des réseaux géants, on PROUVE la réduction à une percolation
  indépendante (E3) et l'universalité suit".
- **Exposant sans équivoque (optionnel, HPC)** : `run_convergence.py` est prêt et
  configurable (changer `LS=[16,24,32,48]`, R>=200, NTRAJ>=4000 en tête de
  fichier) → lancer sur cluster pour extrapoler nu->4/3.
- Fil des 3 papiers : théorème `Delta=I+C_bwd` (P2) -> réseau (G1) ; dualité F/M
  (P2 C5) -> signe de rho_c (P3).
- Extensions : autres réseaux (triangulaire, 3D) ; contagion dynamique ; frontière
  `rho_c(f)=1/[2(1-conv(f))]` déjà analytique (P3-E).

## 8. Commandes

    cd /Users/mostfa/chiral-ring/src/lattice
    ../../.venv/bin/python verify_lattice.py     # 6 gates (obligatoire avant figs)
    ../../.venv/bin/python run_chain.py           # 1D            (~2-3 min)
    ../../.venv/bin/python run_square.py          # 2D percolation (~5-8 min)
    ../../.venv/bin/python run_minority.py        # rho_c + FSS    (~20 min)
    ../../.venv/bin/python run_factorization.py   # loi de factorisation (~14 min)
    ../../.venv/bin/python run_universality.py    # exposants nu(f) (~30 min)
    ../../.venv/bin/python run_convergence.py      # convergence nu/beta/gamma (~30 min)

État : Paper 3 complet côté RÉSULTATS (P3-A..E + exposants). Branche
`paper3-lattice`, NON poussée. Prochaine grande étape = le manuscrit.
Règle : AUCUNE figure tant que les 6 gates ne passent pas. NE PAS pousser sans
accord de l'utilisateur.
