add_library(Solver src/Solverfluid/SolverEuler/SolverEuler.cpp src/Solverfluid/SolverEuler/Cell_Euler.cpp)
target_include_directories(Solver PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/SolverFluid")

target_link_libraries(Solver PUBLIC Eigen3::Eigen Mesh)