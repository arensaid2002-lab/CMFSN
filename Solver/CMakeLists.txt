add_library(Solver src/Solverfluid/SolverEuler.cpp)
target_include_directories(Solver PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/SolverFluid")