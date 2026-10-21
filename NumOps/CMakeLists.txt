add_library(NumOps INTERFACE)
target_include_directories(NumOps INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(NumOps INTERFACE Kokkos::kokkos)