add_library(NumOps INTERFACE)

target_sources(NumOps INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/include/KokkosArrayOps.h"
)

target_include_directories(NumOps INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include/KokkosArrayOps.h")

target_link_libraries(NumOps INTERFACE Kokkos::kokkos)