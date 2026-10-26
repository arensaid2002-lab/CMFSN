# ---------------------------------------------------------------------------
# Tests: one executable per test, registered with CTest.
# Run them from the VS Code Testing tab or with `ctest --test-dir build`.
# Paths to .cpp files are relative to this folder (tests/).
# ---------------------------------------------------------------------------

# --- CAD ---
add_executable(test_version test_version.cpp)
target_link_libraries(test_version PRIVATE nacawing_core)
add_test(NAME version COMMAND test_version)

add_executable(test_naca_equation test_naca_equation.cpp)
target_link_libraries(test_naca_equation PRIVATE nacawing_core)
add_test(NAME naca_equation COMMAND test_naca_equation)

add_executable(test_wing test_wing.cpp)
target_link_libraries(test_wing PRIVATE nacawing_core)
add_test(NAME wing COMMAND test_wing)

# --- Mesh ---
# The mesh file path is passed as an argument (argv[1]) to the test program.
add_executable(test_plot3d test_plot3d.cpp)
target_link_libraries(test_plot3d PRIVATE mesh_core)
add_test(NAME plot3d COMMAND test_plot3d ${CMAKE_CURRENT_SOURCE_DIR}/data/x.9)

add_executable(test_ogrid test_ogrid.cpp)
target_link_libraries(test_ogrid PRIVATE mesh_core)
add_test(NAME ogrid COMMAND test_ogrid ${CMAKE_CURRENT_SOURCE_DIR}/data/x.9)

add_executable(test_ogrid_wall test_ogrid_wall.cpp)
target_link_libraries(test_ogrid_wall PRIVATE mesh_core)
add_test(NAME ogrid_wall COMMAND test_ogrid_wall)
