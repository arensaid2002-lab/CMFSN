add_library(IO src/FileIO.cpp)
target_include_directories(IO PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(IO PUBLIC Mesh)
