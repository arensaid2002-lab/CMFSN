add_library(Solver STATIC src/SolverHandler/SolverHandler.cpp src/SolverHandler/SolverStatus.cpp)

target_include_directories(Solver PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/SolverHandler")

target_link_libraries(Solver PUBLIC Eigen3::Eigen Kokkos::kokkos)