add_library(IO src/FileIO.cpp)
target_include_directories(IO PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(IO PUBLIC Mesh)


if(Kokkos_ENABLE_CUDA)
    set_source_files_properties(
        src/FileIO.cpp
        PROPERTIES LANGUAGE CUDA
    )
endif()