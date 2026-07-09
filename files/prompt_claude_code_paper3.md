# PROMPT CLAUDE CODE — PAPIER 3 : Percolation de l'alignement des flèches
# Repo chiral-ring. Nouvelle branche: paper3-lattice.
# Le moteur de base est FOURNI et DÉJÀ VALIDÉ (4 gates passés en session Claude) :
# lattice_p3.py + preview_p3.py + fig_p3_preview.png. Pars de lui, ne le réécris pas.

## CONTEXTE ACQUIS (ne pas redécouvrir)
1. Moteur : branches jointes S∈{±1}^L, β_n(S)=sin(2·a_n(S)·(t−Tc_n)),
   topologies chaîne (build_chain) ; contrariants = bains privés chargés (Tc=T).
2. GATES DÉJÀ PASSÉS (à porter en tests unitaires, tolérances incluses) :
   G1 Δ=I+C_bwd sur la loi réseau à L=2 (2.8e-17) ; G2 f=0 → Th.2 papier 1 (2.2e-16) ;
   G3 force brute L=3 vs branches (2.2e-17) ; G4 marginal LOCAL ≡ global en signe
   (100%, corr 0.97) → LE MARGINAL LOCAL EST LA DÉFINITION DE PRODUCTION (scalable).
3. DÉCISION CONCEPTUELLE (issue de la preview, cruciale) :
   la flèche d'un anneau est τ_i = sign(d⟨σ_i^local⟩/dt) — la DIRECTION
   D'ACCUMULATION de ses records accessibles — PAS le signe du ratio d'inférence
   (qui mesure la fixation, positive même pour un bain qui se consume).
   Réconcilier explicitement avec la convention C5 du papier 2 (src/multi) et
   documenter dans notes/arrow_definition.md.
4. LEÇON DE PROTOCOLE (preview) : en protocole instantané, fort couplage partagé
   → récurrences (qualité non-monotone) qui polluent τ. Pour la production :
   PROTOCOLE DE COLLISION (records permanents, monotones), et balayer f
   structurellement via Ns (nombre de témoins partagés), pas via la force.

## SPRINTS
### P3-A — Intégration (1 session)
- Porter lattice_p3.py dans src/lattice/ ; adapter au protocole de collision
  (ε_n = sin(2 a_n Δt) par fenêtre, comme src/ du papier 1).
- Ajouter build_square(Lx, Ly, ...) (réseau carré, liens 4-voisins).
- Tests: G1-G4 + G5 (collision: additivité exacte des D_i pour bains disjoints)
  + G6 (L=2 collision reproduit les nombres C5 du papier 2 à f donné).
- Commit. AUCUNE figure avant que les 6 gates passent.

### P3-B — Chaîne 1D (1-2 sessions)
- Sweep: L∈{8,16,32}, f=Ns/(Np+2Ns) via Ns∈{0..6}, ρ∈{0.1,0.2,0.3,0.4,0.5},
  10 réalisations de désordre, collision, 2e4 trajectoires.
- Sorties: A_signed(f,ρ), fraction convertie, tailles de domaines de τ,
  longueur de corrélation ξ(f). Baseline exacte A(f=0)=1−2ρ comme test continu.
- Figures: fig_p3_chain_alignment.png, fig_p3_domains_chain.png.

### P3-C — Réseau 2D et percolation (2 sessions)
- Carrés 4×4→8×8 (marginal local: coût O(2^5) par site, OK).
- Percolation des domaines alignés: probabilité de cluster percolant vs f,
  scaling en taille finie, estimation f_c(ρ) et exposants (comparaison
  percolation standard / Ising 2D).
- Figures: fig_p3_2d_phase.png (carte A(f,ρ)), fig_p3_percolation.png,
  fig_p3_snapshots.png (cartes de τ sous/au seuil), fig_p3_fss.png.

### P3-D — La minorité pondérée + consolidation (1 session)
- ρ_c(f): densité critique de contrariants qui renverse la flèche GLOBALE ;
  tester la prédiction "témoins partagés votent double" → ρ_c < 1/2 aux grands f,
  lien quantitatif avec le f*≈0.40 de la paire.
- Barres d'erreur partout (std sur désordre), graines fixes, une commande
  par figure. Rapport final: valeurs, R², écarts aux baselines exactes.

## GARDE-FOUS
- Chaque quantité nouvelle : d'abord un test contre un cas exactement calculable.
- Toute non-monotonie inattendue : suspecter les récurrences avant la physique.
- Ne pas pousser sans mon accord ; commits par sprint ; log des seeds.
