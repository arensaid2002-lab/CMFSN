cmake_minimum_required(VERSION 3.10.0)
project(CMFSN VERSION 0.1.0 LANGUAGES C CXX)

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
add_subdirectory(Mesh)
add_subdirectory(Solver)
add_subdirectory(UI)


##===========Executable===========##
add_executable(CMFSN Main/main.cpp)


##===========Tests===========##
include(CTest)
enable_testing()

