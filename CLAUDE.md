# CLAUDE.md — CMFSN

## Règles de collaboration (prioritaires)

- **Rôle principal : revue de code et suggestions.** Explique, critique, propose.
- **N'implémente rien sans mon autorisation explicite.** Aucune création, modification ou suppression de fichier sans mon accord. En cas de doute, demande.
- **Avançons pas à pas** : une question ou un sujet à la fois. Propose un plan court avant tout changement non trivial et attends ma validation.
- Quand tu proposes un changement, montre le diff (ou l'extrait de code) et explique **pourquoi**, en particulier les raisons numériques ou physiques.
- Classe tes remarques par priorité : 🔴 bloquant / bug, 🟠 important, 🟡 suggestion.
- Ne lance aucune commande destructive (`rm`, `git reset --hard`, `git push --force`…), et ne fais ni commit ni push sans mon accord.
- Réponds en français.

## Le projet

CMFSN (Mécanique des Fluides et Solides Numériques) est un logiciel de simulation **3D d'écoulements visqueux subsoniques et transsoniques** autour d'avions, avec **déformations structurelles**, c'est-à-dire du couplage fluide-structure et de l'aéroélasticité.

Points d'attention spécifiques au domaine :
- Écoulements compressibles et transsoniques : chocs, schémas de flux, limiteurs, stabilité (CFL).
- Viscosité, couche limite et éventuellement modèles de turbulence.
- Couplage fluide-structure : transfert des charges et des déplacements à l'interface, déformation du maillage fluide, conservation à l'interface.
- Performance : parallélisme (OpenMP / Kokkos), localité mémoire, absence d'allocations dans les boucles chaudes.

## Stack technique

- **C++20**, **CMake**.
- **Eigen** (algèbre linéaire) : trouvé sur le système, sinon téléchargé via FetchContent.
- **Kokkos** (parallélisme portable) : trouvé sur le système, sinon téléchargé via FetchContent. Il exige CMake ≥ 3.22.
- **OpenMP**.
- IDE : VS Code et Visual Studio Community.

## Structure du dépôt

```
CMakeLists.txt        # CMake racine
Main/main.cpp         # point d'entrée (exécutable CMFSN)
CAD/                  # géométrie : InCAD / OutOfCAD   (non encore ajouté au build)
Mesh/                 # maillage : MeshInCAD / MeshOutOfCAD
Solver/
  SolverFluid/        # solveur fluide
  SolverStructure/    # solveur structure
  SolverHandler/      # orchestration / couplage
UI/                   # interface utilisateur
tests/                # tests (non encore ajouté au build)
```

Chaque module suit la convention `include/` (en-têtes) et `src/` (implémentation).

> À compléter : la signification précise de `InCAD` / `OutOfCAD`.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build
```