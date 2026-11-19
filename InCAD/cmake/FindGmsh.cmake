# FindGmsh.cmake
# Cherche l'API C++ de Gmsh : l'en-tete gmsh.h et la bibliotheque libgmsh.
# Cree la cible importee  Gmsh::Gmsh  (a lier avec target_link_libraries).
#
# Si CMake ne trouve pas Gmsh tout seul, donne le dossier du SDK :
#     cmake -S . -B build -DGMSH_ROOT="C:/chemin/vers/gmsh-4.xx.x-Windows64-sdk"
# (ce dossier contient include/gmsh.h et lib/...). Le SDK se telecharge sur gmsh.info.

find_path(GMSH_INCLUDE_DIR
    NAMES gmsh.h
    HINTS ${GMSH_ROOT} ENV GMSH_ROOT
    PATH_SUFFIXES include)

find_library(GMSH_LIBRARY
    NAMES gmsh
    HINTS ${GMSH_ROOT} ENV GMSH_ROOT
    PATH_SUFFIXES lib lib64 bin)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Gmsh REQUIRED_VARS GMSH_INCLUDE_DIR GMSH_LIBRARY)

if(Gmsh_FOUND AND NOT TARGET Gmsh::Gmsh)
    add_library(Gmsh::Gmsh UNKNOWN IMPORTED)
    set_target_properties(Gmsh::Gmsh PROPERTIES IMPORTED_LOCATION "${GMSH_LIBRARY}")
    # SYSTEM : les avertissements (-Wpedantic...) de gmsh.h ne polluent pas la compilation.
    target_include_directories(Gmsh::Gmsh SYSTEM INTERFACE "${GMSH_INCLUDE_DIR}")
endif()

mark_as_advanced(GMSH_INCLUDE_DIR GMSH_LIBRARY)
