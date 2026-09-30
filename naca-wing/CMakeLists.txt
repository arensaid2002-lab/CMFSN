# Minimum CMake version. 3.16 is available on every current OS package manager.
cmake_minimum_required(VERSION 3.16)

project(NacaWing
    VERSION 0.1.0
    DESCRIPTION "NACA airfoil and wing generator"
    LANGUAGES CXX)

# Use C++17 everywhere, and fail rather than silently fall back to an older standard.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to a Release build when the user doesn't pick one (single-config generators only).
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ---------------------------------------------------------------------------
# Library: all the geometry code lives here, so the program and the tests
# (and later a GUI or solver) can share it.
# ---------------------------------------------------------------------------
add_library(nacawing_core
    src/version.cpp
    src/naca.cpp
    src/upperSurface.cpp
)

target_include_directories(nacawing_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Turn on useful warnings. MSVC and GCC/Clang spell the flags differently.
# /Zc:__cplusplus makes MSVC report the real C++ standard in __cplusplus.
if(MSVC)
    target_compile_options(nacawing_core PUBLIC /W4 /permissive- /Zc:__cplusplus)
else()
    target_compile_options(nacawing_core PUBLIC -Wall -Wextra -Wpedantic)
endif()

# ---------------------------------------------------------------------------
# Executable: the command-line program.
# ---------------------------------------------------------------------------
add_executable(naca_wing src/main.cpp)
target_link_libraries(naca_wing PRIVATE nacawing_core)

# ---------------------------------------------------------------------------
# Tests: run with `ctest --test-dir build`.
# ---------------------------------------------------------------------------
enable_testing()
add_executable(test_version tests/test_version.cpp)
target_link_libraries(test_version PRIVATE nacawing_core)
add_test(NAME version COMMAND test_version)

add_executable(test_naca_equation tests/test_naca_equation.cpp)
target_link_libraries(test_naca_equation PRIVATE nacawing_core)
add_test(NAME naca_equation COMMAND test_naca_equation)

