add_library(Solver STATIC src/SolverHandler.cpp src/SolverStatus.cpp)
target_include_directories(Solver PUBLIC include)

target_link_libraries(Solver PUBLIC Eigen3::Eigen Kokkos::kokkos)