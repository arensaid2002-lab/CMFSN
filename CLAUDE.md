# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.


# CMFSN — Solveur fluide (branche SolveurFluid---Said)

Logiciel C++20 de mécanique des fluides numérique pour avions transsoniques (à terme couplé à la structure).
Cette branche contient le solveur FVM 2D structuré : Euler2D et Navier-Stokes laminaire (`RANS2D`), parallélisé avec Kokkos (OpenMP ou CUDA).

## Façon de travailler (prioritaire)
- Avancer **pas à pas** : proposer, expliquer, attendre la validation.
- **Ne rien implémenter ni modifier sans autorisation explicite.** Une revue ou une suggestion n'est pas une autorisation.
- Pour chaque problème signalé : fichier:ligne, scénario concret qui le déclenche, correction proposée.
- Répondre en français.

## Build et tests
Les presets sont dans `CMakePresets.json` : `openmp`, `cuda`, `openmp-debug`, `cuda-debug`.
```
cmake --preset openmp -DBUILD_TESTING=ON
cmake --build build/openmp-clang --target CMFSN CMFSNTests
ctest --preset openmp                 # 32 tests (16 cas × float/double)
ctest --preset openmp -L gradients    # filtrer par famille ou par précision
```
Build CUDA (c'est le build `build/cuda` qui existe sur cette machine) :
```
cmake --preset cuda -DBUILD_TESTING=ON
cmake --build build/cuda --target CMFSN CMFSNTests
ctest --test-dir build/cuda --output-on-failure
```
Lancer un seul test :
```
ctest --test-dir build/cuda -R fvm_rk_ck_double
build/cuda/tests/CMFSNTests.exe rk_ck double      # <cas> <float|double>
```
- `cmake`, `ninja` et `cl` ne sont pas dans le PATH : charger d'abord `vcvars64.bat` de Visual Studio 18 (`VC\Auxiliary\Build`), et ajouter le CMake et le Ninja fournis par Visual Studio (`Common7\IDE\CommonExtensions\Microsoft\CMake`).
- Ajouter un cas de test se fait à deux endroits : la liste `foreach` de `tests/CMakeLists.txt` et `dispatch` dans `tests/FVMTests.cpp`.
- `near` (tests) est une tolérance absolue en dessous de 1. Les grandeurs de l'ordre de `mu` (flux visqueux) se comparent en unités de `mu`, sinon le test ne voit rien en `float`.
- Kokkos 5.2.2 (version exacte) et Eigen 5.0.1 sont téléchargés par FetchContent s'ils sont absents.
- Tout changement numérique doit passer en `float` **et** en `double`.
- Les paramètres de simulation sont codés en dur dans `Main/main.cpp`. Les résultats vont dans `results/run_XXXX/` (`residuals.csv`, `.vts`, `.pvd`).

## Carte du code
- `Solver/include/SolverFluid/SolverFVM/`
  - `SolverFVM.h/.tpp` : boucle `run`/`run_until`, diagnostics, résidu RMS
  - `CellsHandler/CellsData.h` : toutes les `View` Kokkos (SoA), `FlowConditions` (conditions choisies par l'utilisateur) et `initialization_data`, qui en hérite et ajoute l'état amont dérivé
  - `CellsHandler/CellsHandler.h/.tpp` : validation des conditions, géométrie, état initial et fantômes
  - `CellOps/CellOps.h` : géométrie, conditions aux limites, flux, JST, viscosité, RK5 Jameson, RK Carpenter-Kennedy
  - `CellOps/CellCopy.h` : copie de cellules (périodicité)
  - `NumOps/` : opérateurs sur `Kokkos::Array`, thermodynamique, τ
- `IO/` : import du maillage, export VTS/PVD. `Mesh/` : `MeshData_Structured`. `tests/FVMTests.cpp` : tests CTest.
- `UI/visualize.py` : visualisation (pyvista + PySide6). Champs, convergence, distribution de Cp à la paroi, Cl, Cd, Cm intégrés depuis la rangée j = 0 et leur évolution (colonnes `cl, cd, cm` de `residuals.csv` si présentes, sinon recalculés sur les états de la série), lecture au survol, thèmes sombre/clair (couleurs centralisées dans `THEMES`), export des graphes (PNG, SVG, PDF, CSV) et GIF de la série (`--gif`, `--save-curve`).
  ```
  python -m pip install -r UI/requirements.txt
  python UI/visualize.py            # ouvre le results/run_*/solution le plus récent
  ```
  Il n'a pas de tests automatiques.
- Le solveur et `CellsHandler` ont un seul constructeur : `(FlowConditions, maillage)`. `Main/RunSimulation.h` relie maillage, solveur et export ; `SimulationParameters` contient `flow` et les réglages du calcul.

## Conventions à respecter
- **Indexation** : `k = i + Ni_total_cells * j`. Il y a 3 couches fantômes de chaque côté (`CellsData<T>::N_ghost_layers`) : elles sont nécessaires au stencil JST à 6 points. Une cellule physique a `i, j ∈ [3, 3 + N_physical)`.
- **Faces** : 0 = bas (j−), 1 = droite (i+), 2 = haut (j+), 3 = gauche (i−). La face opposée est `(d+2)%4`. Les normales sont unitaires et sortantes. `Normals(k,d,2)` contient la longueur d'arête.
- **Topologie** : maillage en O. j = 0 est la paroi, j = Nj−1 le champ lointain (invariants de Riemann), et la direction i est périodique.
- **Ordre des conditions aux limites** : paroi → champ lointain → périodicité, pour que les coins soient remplis.
- **Adimensionnement** : ρ∞ = p∞ = 1, c∞ = √γ, θ = p/ρ (R = 1), μ = M∞√γ/Re, k = γ/(γ−1)·μ/Pr. `alpha` est en radians.
- **État conservatif** : `*_old` est l'état courant, `*_new` est la sortie de l'étage. On échange les handles avec `std::swap` après chaque étage. `D_history` sert aussi d'accumulateur `G` pour Carpenter-Kennedy.
- **Reynolds** : `Euler2D` ignore le Reynolds fourni (`mu_inf = k_inf = 0`). L'export n'écrit `reynolds_inf` dans le `.vts` que s'il est fini ; le visualiseur s'en sert pour la part de frottement du Cd.
- **Métriques** : tous les `DJacobien` sont calculés une seule fois dans `eval_GeometricalData`, y compris celui de la face commune dans les rangées fantômes adjacentes à la paroi et au champ lointain (opposé à celui de la cellule physique). Les conditions aux limites n'y touchent pas. Sans ces métriques fantômes, le flux visqueux de bord est divisé par deux.
- **Fichiers du maillage** : en-tête `Ni Nj` en ASCII, puis tous les X, puis tous les Y. Ce n'est pas du vrai PLOT3D multibloc.

## Règles Kokkos / GPU
- Le code dans `KOKKOS_LAMBDA` / `KOKKOS_FUNCTION` doit compiler pour le device : utiliser `Kokkos::sqrt`, `Kokkos::abs`, `Kokkos::pow`, jamais les versions `std::` ou non qualifiées. Le solveur et l'export utilisent `Kokkos::` partout, code hôte compris.
- Écrire les littéraux avec `T(0.5)` plutôt que `0.5` (sinon calcul FP64 sur GPU en `float`).
- Les lambdas capturent par valeur. `CellsData` se copie à moindre coût (ce sont des `View`). Ne jamais capturer de référence vers un objet hôte.
- Un `parallel_for` ne doit jamais écrire là où un autre thread lit dans le même kernel : sources et destinations disjointes.
- Mettre un `fence` avant toute lecture côté hôte (`create_mirror_view_and_copy`).

## Points sensibles connus (issus de la revue)
- `NumOps/NumFunctions.h` contient encore des littéraux `0.5` et `1.0` (calcul FP64 sur GPU en `float`), contrairement à la règle `T(0.5)`.
- Le frottement pariétal n'est validé sur aucun cas physique (plaque plane, par exemple), et l'ordre de convergence des schémas n'est pas mesuré.
- Le critère d'arrêt en stationnaire est absolu (`rms <= 1e-6`).
- Paroi isotherme : `2·θw − θint` peut devenir ≤ 0, et JST crée un flux de masse à travers la paroi.
- Aucune vérification d'aire positive ni d'orientation du maillage à l'import.
- Chaque flux de face est calculé deux fois. La périodicité copie toutes les données (`All`) à chaque étage.

## À éviter
- Ne pas ouvrir les maillages volumineux (`NACA_0012_mesh_PLOT3D_format/x.1`–`x.3`, jusqu'à 80 Mo). Pour tester, utiliser `x.9` (9×9) ou `Oplot/NACA0012_nc8.x`.
- Ne pas ajouter de `using namespace` dans les en-têtes (il y en a déjà trop).
- Ne pas affaiblir un test pour le faire passer : expliquer pourquoi il échoue.
- Les fichiers sont en CRLF dans la copie de travail (`core.autocrlf=true`). Un `sed -i` sous Git Bash les convertit en LF : repasser `unix2dos` ensuite.
