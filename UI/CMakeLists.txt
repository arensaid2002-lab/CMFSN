# ---------------------------------------------------------------------------
# UI: the graphical programs. They sit here, above CAD and Mesh, so a viewer
# can link to both libraries (Mesh depends on CAD, never the other way).
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# FETCHCONTENT: raylib, a simple graphics library for the GUI.
# ---------------------------------------------------------------------------
include(FetchContent)
FetchContent_Declare(raylib
    GIT_REPOSITORY https://github.com/raysan5/raylib.git
    GIT_TAG 5.5            # a fixed version, so the build never changes by surprise
    GIT_SHALLOW TRUE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)   # don't build raylib's own examples
FetchContent_MakeAvailable(raylib)


FetchContent_Declare(raygui
    GIT_REPOSITORY https://github.com/raysan5/raygui.git
    GIT_TAG 4.0            # a fixed version, so the build never changes by surprise
    GIT_SHALLOW TRUE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)   # don't build raylib's own examples
FetchContent_MakeAvailable(raygui)

# NACA viewer: the airfoil surfaces only, so it needs CAD.
add_executable(naca_viewer src/naca_viewer.cpp)
target_include_directories(naca_viewer SYSTEM PRIVATE ${raygui_SOURCE_DIR}/src)
target_link_libraries(naca_viewer PRIVATE nacawing_core raylib)

# Mesh viewer: the O-grid, so it needs Mesh (which brings CAD with it).
add_executable(mesh_viewer src/mesh_viewer.cpp)
target_include_directories(mesh_viewer SYSTEM PRIVATE ${raygui_SOURCE_DIR}/src)   # raygui: the buttons
target_link_libraries(mesh_viewer PRIVATE mesh_core raylib)
# Where the reference grids x.9 ... x.5 (Vassberg-Jameson PLOT3D files) live, baked into the program.
target_compile_definitions(mesh_viewer PRIVATE REF_DIR="${CMAKE_SOURCE_DIR}/tests/data")
