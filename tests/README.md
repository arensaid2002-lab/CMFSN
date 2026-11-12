# Tests du solveur

Les tests sont actifs avec `BUILD_TESTING=ON` (valeur par defaut de CTest).
Ils utilisent le backend Kokkos de la configuration choisie et testent
chaque cas en `float` et en `double`.

Depuis un terminal de developpement Visual Studio, a la racine du depot :

```powershell
cmake --preset cuda -DBUILD_TESTING=ON
cmake --build build/cuda --target CMFSN CMFSNTests
ctest --preset cuda
```

Pour le CPU, remplacer `cuda` par `openmp` et le repertoire de compilation
par `build/openmp-clang`. Le preset de test CPU reprend le `PATH` du preset
de configuration pour trouver la bibliotheque d'execution OpenMP.
Les presets de test `cuda-debug` et `openmp-debug` sont egalement disponibles.

Filtrer une famille ou une precision :

```powershell
ctest --preset cuda -L gradients
ctest --preset cuda -L double
```

## Couverture

Les 32 tests couvrent l'indexation, la thermodynamique, la geometrie et les
metriques des quatre faces, les gradients de champs affines sur deux grilles,
les contraintes visqueuses, la periodicite et ses coins, les etats de paroi
adiabatique/isotherme, les quatre regimes du
farfield, la conservation des faces interieures, le pas de temps et les
constructeurs. Deux essais courts verifient l'avancement instationnaire
Jameson et Carpenter-Kennedy. Un export VTS verifie les indices et les valeurs
des cellules physiques.

Ces tests verifient des contrats locaux et des parcours courts. Ils ne
remplacent pas une etude de convergence spatiale/temporelle, un benchmark
laminaire ni une validation du cisaillement pres de la paroi sur une grille
fortement etiree.

## Configuration visqueuse

Le constructeur avec Reynolds et Prandtl accepte a la fin la condition
thermique puis sa temperature. Le nombre de couches fantomes est fixe a trois
par `CellsData<T>::N_ghost_layers` et ne fait plus partie des arguments des
constructeurs. La paroi par defaut est adiabatique.

```cpp
Solver::SolverFVM<T, Solver::SolverTypes::RANS2D, TimeScheme, RKScheme> solver(
    Mach, gamma, alpha, Re, Pr, mesh,
    Solver::BoundaryCondition::isothermalWall, wall_temperature);
```

La temperature utilise la convention adimensionnelle `theta = p/rho`, avec
`theta_inf = 1`. Ce parametre n'est pas automatiquement en kelvins.
La temperature fantome est `2*theta_wall - theta_inside` pour le prolongement
symetrique utilise ici. Elle doit rester positive ; ce traitement ne contient
pas de limiteur de positivite et ne garantit pas l'ordre deux de la fermeture
de paroi sur toute grille deformee.

Le nom historique `RANS2D` active actuellement les flux laminaires a viscosite
constante. Aucune equation de turbulence Spalart-Allmaras ou SST n'est ajoutee.
Le type de solveur et le schema temporel restent explicites dans `Main/main.cpp`.
Les valeurs physiques, le CFL, le temps final, les iterations, les chemins
et les cadences sont directement modifiables dans le bloc `parameters` de
`Main/main.cpp`. Le programme ne lit plus de parametres de simulation dans
les arguments de ligne de commande.

Le lancement utilise exclusivement les valeurs du main. En mode visqueux,
Reynolds et Prandtl doivent etre positifs et finis. La valeur infinie de
Reynolds conserve l'ancien defaut Euler ; il faut la remplacer par la valeur
du cas laminaire. Le selecteur `wall_condition` choisit le type de paroi et
`wall_temperature` donne sa temperature adimensionnelle si elle est isotherme.

## Dimensions

`Ni_nodes`, `Nj_nodes` comptent les sommets physiques ;
`Ni_physical_cells`, `Nj_physical_cells` comptent les cellules physiques ;
`Ni_total_cells`, `Nj_total_cells` incluent les fantomes. `row_stride` est le
saut d'une rangee, donc `Ni_total_cells`.

```cpp
Cell_index(i, j, Ni_total_cells) == i + Ni_total_cells*j;
```

L'ancien parametre de `Cell_index` comptait une dimension de sommets augmentee ;
il fallait lui retrancher un. Le nouveau parametre compte directement les
cellules stockees. Le changement de nom et la suppression des `-1` vont ensemble.

## JST et echecs physiques

Le calcul JST conserve sa formule d'origine sur toutes les faces, y compris
la paroi. Le test `jst` compare les quatre composantes a une reference scalaire
independante, en interieur et en frontiere. Il ne prescrit pas une annulation
de la dissipation au mur.

Les tests ne prescrivent pas une fermeture numerique particuliere du flux
au mur. Le controle d'assemblage dans `viscous_flux` porte sur les faces
interieures ou les deux etats, metriques, normales et flux sont prepares.
Les gradients et contraintes restent compares a des champs analytiques sur
les quatre orientations. La validation du flux complet de paroi reste a faire.
Un echec futur doit etre examine sans modifier automatiquement l'algorithme.
