# ---------------------------------------------------------------------------
# cmfsn: the main program. It reads the inputs written at the top of main.cpp,
# builds the geometry and the meshes, and saves them in Results/.
# It sits above CAD and Mesh, so it can use both (Mesh already links CAD).
# ---------------------------------------------------------------------------
add_executable(cmfsn main.cpp)
target_link_libraries(cmfsn PRIVATE mesh_core)   # mesh_core brings nacawing_core with it (PUBLIC link)

# Full path of the Results folder, fixed at build time, so the program writes to
# the same place whatever folder you start it from.
target_compile_definitions(cmfsn PRIVATE CMFSN_RESULTS_DIR="${CMAKE_SOURCE_DIR}/Results")
