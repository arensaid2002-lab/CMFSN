# ---------------------------------------------------------------------------
# Library: everything about meshes (structured grid, PLOT3D reader, later the
# O-grid generator). Programs and tests link to it, like nacawing_core in CAD.
# ---------------------------------------------------------------------------
add_library(mesh_core
    src/Plot3D.cpp
    src/OGrid.cpp
)

# Anyone who links to mesh_core can #include "StructuredMesh2D.h" and "Plot3D.h".
target_include_directories(mesh_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)


# The mesh is built from the geometry, so Mesh depends on CAD (never the other way).
# PUBLIC because StructuredMesh2D.h itself includes "nacawing/NACASurface.hpp".
target_link_libraries(mesh_core PUBLIC nacawing_core)

