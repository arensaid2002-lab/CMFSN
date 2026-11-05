add_library(Solver src/SolverFluid/SolverFVM/SolverFVM.cpp)

target_sources(Solver PRIVATE
    include/SolverFluid/SolverFVM/SolverFVM.h
    include/SolverFluid/SolverFVM/SolverFVM.tpp
    include/SolverFluid/SolverFVM/CellsHandler/CellsHandler.h
    include/SolverFluid/SolverFVM/CellsHandler/CellsData.h
    include/SolverFluid/SolverFVM/CellOps/CellOps.h
    include/SolverFluid/SolverFVM/CellOps/CellCopy.h
    include/SolverFluid/SolverFVM/CellsHandler/CellsHandler.tpp
    include/SolverFluid/SolverFVM/NumOps/ArrayOps.h
    include/SolverFluid/SolverFVM/NumOps/NumFunctions.h
)

set_source_files_properties(
    include/SolverFluid/SolverFVM/SolverFVM.tpp
    include/SolverFluid/SolverFVM/CellsHandler/CellsHandler.tpp
    PROPERTIES HEADER_FILE_ONLY TRUE
)

target_include_directories(Solver PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/SolverFluid/SolverFVM"

)

target_link_libraries(Solver PUBLIC Eigen3::Eigen Mesh Kokkos::kokkos)

if(Kokkos_ENABLE_CUDA)
    set_source_files_properties(
        src/SolverFluid/SolverFVM/SolverFVM.cpp
        PROPERTIES LANGUAGE CUDA
    )
endif()

if(Kokkos_ENABLE_CUDA AND MSVC)
    # Keep the MSVC host and standalone-header IntelliSense in C++20.
    target_compile_options(Solver PRIVATE
        "$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcompiler /std:c++20>"
    )
endif()
