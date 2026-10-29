add_library(Solver src/Solverfluid/SolverEuler/SolverEuler.cpp src/Solverfluid/SolverEuler/Cells_Euler.cpp)

if(Kokkos_ENABLE_CUDA)
    set_source_files_properties(
        SolverFluid/SolverEuler/src/SolverEuler.cpp
        SolverFluid/SolverEuler/src/Cells_Euler.cpp
        PROPERTIES LANGUAGE CUDA
    )
endif()

target_include_directories(Solver PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/SolverFluid")

target_link_libraries(Solver PUBLIC Eigen3::Eigen Mesh)