add_library(Solver src/SolverFluid/SolverEuler/SolverEuler.cpp)

target_sources(Solver PRIVATE
    include/SolverFluid/SolverEuler/Cells_Euler.h
    include/SolverFluid/SolverEuler/CellsData.h
    include/SolverFluid/SolverEuler/CellOps/CellOps.h
    include/SolverFluid/SolverEuler/CellOps/CellCopy.h
    include/SolverFluid/SolverEuler/Cells_Euler.tpp
)

set_source_files_properties(
    include/SolverFluid/SolverEuler/Cells_Euler.tpp
    PROPERTIES HEADER_FILE_ONLY TRUE
)

target_include_directories(Solver PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/SolverFluid/SolverEuler"
)

target_link_libraries(Solver PUBLIC Eigen3::Eigen Mesh Kokkos::kokkos)

if(Kokkos_ENABLE_CUDA)
    set_source_files_properties(
        src/SolverFluid/SolverEuler/SolverEuler.cpp
        PROPERTIES LANGUAGE CUDA
    )
endif()

if(Kokkos_ENABLE_CUDA AND MSVC)
    # Keep the MSVC host and standalone-header IntelliSense in C++20.
    target_compile_options(Solver PRIVATE
        "$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcompiler /std:c++20>"
    )
endif()
