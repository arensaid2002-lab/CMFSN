# ---------------------------------------------------------------------------
# InCAD : noyau CAD (OpenCASCADE via Gmsh). Il importe les sections NACA calculees
# par la Phase 2 (CAD/ : sectionLoop + sectionAt), en fait des B-splines, les relie
# en une aile (loft), puis la maille avec Gmsh.
#
# Deux facons de l'utiliser :
#   1. Dans le projet complet : la racine fait add_subdirectory(InCAD).
#   2. Seul :   cmake -S InCAD -B build-incad   (il va chercher ../CAD et ../Mesh)
#
# Si Gmsh est introuvable, le module est simplement saute (le reste du projet compile).
# ---------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.16)
project(CMFSN_InCAD VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---- 1. Gmsh (gmsh.h + libgmsh) ------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
find_package(Gmsh)
if(NOT Gmsh_FOUND)
    message(WARNING "InCAD : Gmsh introuvable (il faut gmsh.h et libgmsh). Module InCAD ignore.\n"
                    "         Donne le SDK avec -DGMSH_ROOT=<dossier du SDK Gmsh>.")
    return()
endif()
message(STATUS "InCAD : Gmsh trouve (${GMSH_LIBRARY})")

# ---- 2. Les bibliotheques (en mode autonome seulement, si le CMake du projet complet ne marche pas, celui s'active ) -------------
if(NOT TARGET nacawing_core)
    add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../CAD  ${CMAKE_CURRENT_BINARY_DIR}/_cad)
endif()
if(NOT TARGET mesh_core)
    add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../Mesh ${CMAKE_CURRENT_BINARY_DIR}/_mesh)
endif()

# Dossier Results/ (le meme que celui de Main/ : <depot>/Results)
set(CMFSN_RESULTS_DIR "${CMAKE_CURRENT_LIST_DIR}/../Results" CACHE PATH "Dossier Results")

# ---- 3. La bibliotheque incad_core ----------------------------------------
add_library(incad_core src/InCAD.cpp)
target_include_directories(incad_core PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(incad_core PUBLIC nacawing_core Gmsh::Gmsh)

# ---- 4. Le programme : sections de la Phase 2 -> aile CAD -> maillage ------
add_executable(cmfsn_incad src/main_incad.cpp)
target_link_libraries(cmfsn_incad PRIVATE incad_core mesh_core)
target_compile_definitions(cmfsn_incad PRIVATE CMFSN_RESULTS_DIR="${CMFSN_RESULTS_DIR}")

# ---- 5. Le test (ctest) -----------------------------------------------------
add_executable(test_incad tests/test_incad.cpp)
target_link_libraries(test_incad PRIVATE incad_core mesh_core)
enable_testing()
add_test(NAME incad COMMAND test_incad)

# ---- 6. Windows : la DLL de Gmsh doit etre a cote des .exe dans meme dossier -----------------
# (si la DLL n'est pas copier, copie-la a la main)
if(WIN32)
    file(GLOB GMSH_DLLS "${GMSH_ROOT}/bin/*.dll" "${GMSH_ROOT}/lib/*.dll")
    if(GMSH_DLLS)
        foreach(tgt cmfsn_incad test_incad)
            add_custom_command(TARGET ${tgt} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${GMSH_DLLS} $<TARGET_FILE_DIR:${tgt}>)
        endforeach()
    endif()
endif()
