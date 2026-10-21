add_library(Mesh src/MeshHandler/MeshHandler.cpp)

target_include_directories(Mesh PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/MeshHandler")
target_link_libraries(Mesh PUBLIC Kokkos::kokkos)