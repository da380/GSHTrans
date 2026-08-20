cmake_minimum_required(VERSION 3.20)

project(GSHTrans
        VERSION 1.0.0
        DESCRIPTION "Generalized spherical harmonic transforms and a tensor field algebra on the sphere"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

include(FetchContent)
include(GNUInstallDirs)

# Dependencies.
#
# Each is looked for on the system first and fetched only if it is not there,
# so that a distribution package or an existing checkout is used when one
# exists and nothing is downloaded twice. FetchContent_MakeAvailable does the
# find_package half itself when given FIND_PACKAGE_ARGS.
#
# All URLs are https. GaussQuad was previously fetched over SSH, which meant
# the project could not be configured at all without GitHub keys.
FetchContent_Declare(
  NumericConcepts
  GIT_REPOSITORY https://github.com/da380/NumericConcepts.git
  GIT_TAG main
  FIND_PACKAGE_ARGS)

FetchContent_Declare(
  GaussQuad
  GIT_REPOSITORY https://github.com/da380/GaussQuad.git
  GIT_TAG main
  FIND_PACKAGE_ARGS)

FetchContent_Declare(
  FFTWpp
  GIT_REPOSITORY https://github.com/da380/FFTWpp.git
  GIT_TAG main
  FIND_PACKAGE_ARGS)

FetchContent_MakeAvailable(NumericConcepts GaussQuad FFTWpp)

find_package(OpenMP REQUIRED)

# The library itself: headers only.
add_library(GSHTrans INTERFACE)
add_library(GSHTrans::GSHTrans ALIAS GSHTrans)

target_link_libraries(GSHTrans INTERFACE
  NumericConcepts
  FFTWpp
  GaussQuad
  OpenMP::OpenMP_CXX)

target_compile_features(GSHTrans INTERFACE cxx_std_23)

target_include_directories(GSHTrans INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Installation.
#
# The headers install; the CMake package does not, and cannot yet. An
# INTERFACE target can only be exported if everything it links is itself
# exported or imported, and two of the three dependencies -- FFTWpp and
# GaussQuad -- carry no install rules at all, so there is nothing for a
# consumer's find_package to find. Shipping a config file that calls
# find_dependency on them would produce a package that fails at configure
# time, which is worse than not shipping one.
#
# So: install the headers, which is useful on its own, and leave the package
# for when the dependencies can be installed. Adding install and export rules
# to FFTWpp and GaussQuad is a small change in each and is the whole of what
# blocks this; docs/thoughts.md section 5 records it as a hand-over item.
# Until then, consumption is by add_subdirectory or FetchContent, which works
# today.
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/GSHTrans
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING
          PATTERN "*.h"
          PATTERN "All"
          PATTERN "Core"
          PATTERN "Field"
          PATTERN "Tensor"
          PATTERN "Expansion")

# Developer targets, only when this is the top-level project.
if(PROJECT_IS_TOP_LEVEL)
  option(GSHTRANS_BUILD_EXAMPLES   "Build the examples"   ON)
  option(GSHTRANS_BUILD_TESTS      "Build the test suite" ON)
  option(GSHTRANS_BUILD_BENCHMARKS "Build the benchmarks" ON)

  if(GSHTRANS_BUILD_EXAMPLES)
    add_subdirectory(examples)
  endif()
  if(GSHTRANS_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
  endif()
  if(GSHTRANS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
  endif()
endif()
