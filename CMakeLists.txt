cmake_minimum_required(VERSION 3.10.0)
project(CMFSN VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(Kokkos_ENABLE_CUDA)
    enable_language(CUDA)

    set(CMAKE_CUDA_STANDARD 20)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

    set(
        Kokkos_ENABLE_COMPILE_AS_CMAKE_LANGUAGE
        ON CACHE BOOL "" FORCE
    )
endif()
##===========Dependencies===========##

find_package(Eigen3 CONFIG QUIET)
find_package(Kokkos 5.2.2 EXACT CONFIG QUIET)
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
        SYSTEM
    )
    FetchContent_MakeAvailable(Kokkos)
ENDIF()

if(MSVC)
    foreach(kokkos_target IN ITEMS kokkoscore kokkoscontainers)
        target_compile_definitions(${kokkos_target} PRIVATE
            _CRT_SECURE_NO_WARNINGS
        )

        target_compile_options(${kokkos_target} PRIVATE
            /wd4244
            /wd4267
        )
    endforeach()
endif()


##===========subdirectories===========##
add_subdirectory(IO)
add_subdirectory(Mesh)
add_subdirectory(Solver)
add_subdirectory(NumOps)
add_subdirectory(UI)

##===========Executable===========##
add_executable(CMFSN Main/main.cpp)

target_link_libraries(CMFSN PRIVATE IO Mesh NumOps Solver)
##===========Tests===========##
include(CTest)
enable_testing()


