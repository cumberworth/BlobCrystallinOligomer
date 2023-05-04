cmake_minimum_required(VERSION 3.15)

project(
  BlobCrystallinOligomer
  VERSION 0.1.0
  LANGUAGES CXX)

# Build type
set(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(
    STATUS
      "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "${DEFAULT_BUILD_TYPE}"
      CACHE STRING "Choose the type of build." FORCE)

  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(
  BlobCrystallinOligomer_lib
  src/config.cpp
  src/energy.cpp
  src/ifile.cpp
  src/monomer.cpp
  src/movetype.cpp
  src/ofile.cpp
  src/param.cpp
  src/particle.cpp
  src/potential.cpp
  src/random_gens.cpp
  src/simulation.cpp
  src/space.cpp)

target_compile_features(BlobCrystallinOligomer_lib PRIVATE cxx_std_17)
target_include_directories(BlobCrystallinOligomer_lib PUBLIC include)

# Boost
find_package(Boost 1.69 REQUIRED COMPONENTS program_options)

message(STATUS "Boost version: ${Boost_VERSION}")
target_link_libraries(BlobCrystallinOligomer_lib PUBLIC Boost::program_options)

# Interprocedular optimization
include(CheckIPOSupported)
check_ipo_supported(RESULT RESULT)
if(RESULT)
  message(STATUS "Using interprocedural optimization.")
  set_target_properties(BlobCrystallinOligomer_lib
                        PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# CCache (increased compiliation speed)
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()

# Compiler warnings

# target_compile_options(BlobCrystallinOligomer_lib PUBLIC -Wall -Wextra
# -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused
# -Woverloaded-virtual -Wpedantic -Wconversion -Wsign-conversion
# -Wnull-dereference -Wdouble-promotion -Werror)

# target_compile_options(BlobCrystallinOligomer_lib PUBLIC -Wall -Wextra
# -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Woverloaded-virtual
# -Wpedantic -Wconversion -Wsign-conversion -Wnull-dereference
# -Wdouble-promotion)

# Main simulation program
add_executable(blobCrystallinOligomer apps/main.cpp)
target_link_libraries(blobCrystallinOligomer PUBLIC BlobCrystallinOligomer_lib)

# Testing
find_package(Catch2 REQUIRED)
add_executable(tests test/test_main.cpp test/test_config.cpp
                     test/test_particle.cpp test/test_potential.cpp)
target_link_libraries(tests BlobCrystallinOligomer_lib Catch2::Catch2)
#include(CTest)
#include(Catch)
#catch_discover_tests(tests)
