cmake_minimum_required(VERSION 3.16)
project(CMFSN VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
##===========Dependencies===========##

find_package(Eigen3 CONFIG QUIET)
find_package(Kokkos CONFIG QUIET)
find_package(OpenMP REQUIRED COMPONENTS CXX)


IF (NOT Eigen3_FOUND)
    message(STATUS "EIgen not found, downloading and building Eigen")
    include(FetchContent)
    FetchContent_Declare(        
        Eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 5.0.1)
    FetchContent_MakeAvailable(Eigen)
ENDIF()
IF(NOT Kokkos_FOUND)
    message(STATUS "Kokkos not found, downloading and building Kokkos")
    include(FetchContent)
    FetchContent_Declare(
        Kokkos
        GIT_REPOSITORY https://github.com/kokkos/kokkos.git
        GIT_TAG 5.2.2
    )
    FetchContent_MakeAvailable(Kokkos)
ENDIF()

##===========subdirectories===========##
add_subdirectory(CAD)
add_subdirectory(Mesh)
add_subdirectory(Solver)
add_subdirectory(UI)


##===========Tests===========##
include(CTest)
enable_testing()
add_subdirectory(tests)   # every test is defined in tests/CMakeLists.txt
