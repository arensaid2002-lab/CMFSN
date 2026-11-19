# InCAD : noyau CAD Gmsh + import des sections NACA

Ce module prend les **sections NACA de la Phase 2** (out-of-CAD), les importe dans le noyau CAD
OpenCASCADE de Gmsh, construit l'aile (loft) et la maille.

```
Phase 2 : sectionLoop + sectionAt   ->  13 sections (corde, flèche, vrillage déjà appliqués)
InCAD   : addPoint + addSpline      ->  1 B-spline par section (passe par les 201 points)
          addThruSections           ->  solide de l'aile (demi-envergure)
Gmsh    : generate(2)               ->  maillage de surface (triangles)
```

Les stations d'envergure sont celles du VLM et de la poutre (`vlmNj + 1 = 13`, `SineTip`).

## Contenu (aucun fichier d'Eddie n'est modifié)

| Fichier | Rôle |
|---|---|
| `include/incad/InCAD.h`, `src/InCAD.cpp` | la bibliothèque `incad_core` |
| `src/main_incad.cpp` | programme `cmfsn_incad` (mêmes entrées que `Main/main.cpp`) |
| `tests/test_incad.cpp` | test `ctest` |
| `cmake/FindGmsh.cmake` | trouve `gmsh.h` + `libgmsh` |
| `CMakeLists.txt` | construit le tout |

## Installer Gmsh pour le C++

Il faut l'en-tête `gmsh.h` et la bibliothèque `libgmsh`, avec OpenCASCADE inclus.
Le plus simple : le **SDK** de gmsh.info (dossier contenant `include/` et `lib/`).

```
cmake -S . -B build -DGMSH_ROOT="C:/chemin/vers/gmsh-4.xx.x-Windows64-sdk"
```

Si Gmsh est introuvable, CMake affiche un avertissement et **saute** InCAD : le reste du projet compile.

## Compiler et lancer

**(FAIT) Dans le projet complet** : ajouter UNE ligne dans le `CMakeLists.txt` racine, avec les autres
`add_subdirectory`, avant `include(CTest)` :

```cmake
add_subdirectory(InCAD)
```

puis les commandes habituelles :

```
cmake -S . -B build
cmake --build build
build/cmfsn_incad          (Windows : build\Debug\cmfsn_incad.exe)
ctest --test-dir build -C Debug
```

**Seul** (sans toucher à la racine) : `cmake -S InCAD -B build-incad`, puis les mêmes commandes.

## Fichiers produits

- `Results/Geometry/wing_NACA0012.step` : la géométrie CAD
- `Results/Mesh/InCAD/wing_NACA0012_incad.msh` : le maillage (groupes `Skin`, `Root`, `Tip`, `Wing`)
- `Results/Mesh/InCAD/incad_report.txt` : les chiffres ci-dessous

## Résultats attendus (Gmsh 4.15.2, Linux)

```
sections importées (Phase 2)      : 13 x 201 points
surfaces de peau                  : 2
écart max points importés/surface : 3.0e-07 m
volume (à partir du maillage)     : 0.19118 m3      (formule : 0.19179)
aire de la face d'emplanture      : 0.08184 m2      (formule : 0.08220)
maillage : 45033 noeuds, 83804 triangles
```

L'écart, le volume attendu et les deux faces planes sont reproductibles. Le **nombre de noeuds et de
triangles peut changer** d'une version de Gmsh à l'autre (le mailleur évolue) : compare plutôt
l'écart max et le volume, qui doivent rester dans les tolérances du test (1e-4 m et 0,5 %).

## Limites connues

- Le volume vient du maillage, pas de `gmsh::model::occ::getMass` : sur cette aile, `getMass`
  donne un volume faux de 19 % (la géométrie, elle, est correcte : voir `InCAD.h`).
- Le bord de fuite doit être épais (NACA standard). Un bord de fuite fermé lance une exception.
- Seul le maillage de surface est testé (`MeshOptions::dim = 2`). Le maillage de volume (`dim = 3`) existe mais n'est pas validé.
- Testé uniquement sous Linux. Windows (copie de la DLL de Gmsh à côté des `.exe`, dans le CMake) n'est pas testé.
