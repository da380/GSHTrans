cmake_minimum_required(VERSION 3.20)

project(GSHTrans
        VERSION 1.0.0
        DESCRIPTION "Generalized spherical harmonic transforms and a tensor field algebra on the sphere"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# The version in GSHTrans/src/Version.h is written by hand, so that copying
# the GSHTrans directory into another project is a complete way to vendor the
# library. Check it against the version declared above, so that the two cannot
# drift apart.
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/GSHTrans/src/Version.h _version_header)
set(_header_version)
foreach(_part MAJOR MINOR PATCH)
  string(REGEX MATCH "define GSHTRANS_VERSION_${_part} ([0-9]+)"
         _match "${_version_header}")
  list(APPEND _header_version ${CMAKE_MATCH_1})
endforeach()
list(JOIN _header_version "." _header_version)
if(NOT _header_version STREQUAL PROJECT_VERSION)
  message(FATAL_ERROR
    "Version mismatch: CMakeLists.txt declares ${PROJECT_VERSION} but "
    "GSHTrans/src/Version.h declares ${_header_version}. Update both, "
    "including GSHTRANS_VERSION_STRING.")
endif()
if(NOT _version_header MATCHES "GSHTRANS_VERSION_STRING \"${PROJECT_VERSION}\"")
  message(FATAL_ERROR
    "GSHTRANS_VERSION_STRING does not match ${PROJECT_VERSION}.")
endif()
unset(_version_header)
unset(_header_version)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# No FindFFTW.cmake here, and no CMAKE_MODULE_PATH entry for one.
#
# FFTWpp ships its own, written around imported targets, and asks for
# FFTW::Double. This project used to prepend its own module directory holding
# an older variable-setting FindFFTW.cmake, which won the search and set
# FFTW_LIBRARIES instead -- so configuring against a current FFTWpp failed
# outright with "target FFTW::Double not found". GSHTrans uses no FFTW symbol
# directly, so it has no business finding FFTW at all.

include(FetchContent)
include(GNUInstallDirs)

option(GSHTRANS_INSTALL "Generate the install and export rules"
       ${PROJECT_IS_TOP_LEVEL})

# A dependency that was fetched rather than found has to install itself into
# the same prefix, or the package written below refers to packages that are
# not there and find_package(GSHTrans) fails at the first find_dependency.
# GaussQuad and NumericConcepts install by default; FFTWpp does not when it is
# not the top-level project, so ask it to.
#
# Set without FORCE, so that an explicit -DFFTWPP_INSTALL=OFF still wins, and
# harmless when the dependency is found on the system rather than fetched --
# in that case nothing reads the variable.
if(GSHTRANS_INSTALL)
  set(FFTWPP_INSTALL ON CACHE BOOL "Install FFTWpp alongside GSHTrans")
endif()

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

# Interpolation, which is optional and is the one dependency only part of the
# library uses.
#
# Nothing angular wants a one-dimensional interpolator: the transform, the
# index algebra and the surface gradient are complete without it, and so are
# all three ready-made radial derivatives, which are self-contained by
# design (field-algebra-plan.md section 19.5 [R4]). What it buys is
# resampling a layered field onto a different radial grid, where the whole
# menu of schemes comes for the price of a policy argument.
#
# On by default because it is header-only, depends on nothing but
# NumericConcepts -- which is already here -- and is small. Off is a supported
# configuration and the suite runs in it.
#
# Tracking main, like the three above. Interpolation's own README asks
# consumers to pin a reviewed commit, which is right for a stranger and wrong
# for four libraries with one author that move together (section 19.5 [R2]).
option(GSHTRANS_WITH_INTERPOLATION
       "Radial resampling, built on da380/Interpolation" ON)

if(GSHTRANS_WITH_INTERPOLATION)
  FetchContent_Declare(
    Interpolation
    GIT_REPOSITORY https://github.com/da380/Interpolation.git
    GIT_TAG main
    FIND_PACKAGE_ARGS)
  FetchContent_MakeAvailable(Interpolation)
endif()

# If a dependency was fetched rather than found, it is an ordinary target in
# this build and must be in an export set of its own for GSHTrans to be
# exportable at all. CMake does notice, but only at the generate step and with
# a message that says nothing about how to fix it -- and the commonest way to
# arrive here is a build tree configured before this file asked for the
# install, whose cache still holds the old answer. So say it plainly.
if(GSHTRANS_INSTALL)
  # Interpolation is not in this list, and that is not an omission: its
  # install rules are unconditional, so a fetched copy installs itself with no
  # option to set and there is nothing here to check.
  foreach(_dependency NumericConcepts GaussQuad FFTWpp)
    if(TARGET ${_dependency})
      get_target_property(_imported ${_dependency} IMPORTED)
      string(TOUPPER ${_dependency} _upper)
      if(NOT _imported AND NOT ${_upper}_INSTALL)
        message(FATAL_ERROR
          "GSHTRANS_INSTALL is on, but ${_dependency} was fetched into this "
          "build with ${_upper}_INSTALL off, so it will not install itself "
          "and the exported GSHTrans package would refer to a package that "
          "is not there. Either configure with -D${_upper}_INSTALL=ON, or "
          "turn off -DGSHTRANS_INSTALL=OFF. A build tree configured before "
          "GSHTrans asked for this will have the old value cached; deleting "
          "CMakeCache.txt is the quickest cure.")
      endif()
    endif()
  endforeach()
endif()

find_package(OpenMP REQUIRED)

# The library itself: headers only.
add_library(GSHTrans INTERFACE)
add_library(GSHTrans::GSHTrans ALIAS GSHTrans)

# Namespaced names throughout. An exported target file names its
# dependencies literally, and a consumer resolving them through find_package
# only ever gets the namespaced form -- so linking the bare names here would
# build fine and produce a package that cannot be used. All three dependencies
# provide the alias in add_subdirectory mode as well, so the one spelling
# works in both.
target_link_libraries(GSHTrans INTERFACE
  NumericConcepts::NumericConcepts
  FFTWpp::FFTWpp
  GaussQuad::GaussQuad
  OpenMP::OpenMP_CXX)

# The optional half, and the macro that guards the headers using it. Carried
# on the target rather than generated into a header, because Version.h is
# written by hand so that copying the GSHTrans directory remains a complete
# way to vendor the library -- and a build-tree-generated header would break
# exactly that. An INTERFACE definition is exported with the target, so an
# installed consumer sees the same answer this build did.
if(GSHTRANS_WITH_INTERPOLATION)
  target_link_libraries(GSHTrans INTERFACE Interpolation::Interpolation)
  target_compile_definitions(GSHTrans INTERFACE GSHTRANS_HAVE_INTERPOLATION)
endif()

target_compile_features(GSHTrans INTERFACE cxx_std_23)

target_include_directories(GSHTrans INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Installation, and the CMake package.
#
# This could not be done until recently. An INTERFACE target can only be
# exported if everything it links is itself exported or imported, and two of
# the three dependencies -- FFTWpp and GaussQuad -- carried no install rules
# at all, so a config file calling find_dependency on them would have produced
# a package that fails at configure time. Both now ship packages of their own,
# which is what unblocked this; see docs/FFTWpp-for-GSHTrans.md and
# docs/GaussQuad-for-GSHTrans.md.
#
# Consumption by add_subdirectory or FetchContent still works and is
# unaffected.
if(GSHTRANS_INSTALL)
  include(CMakePackageConfigHelpers)
  set(GSHTRANS_CMAKE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/GSHTrans")

  install(TARGETS GSHTrans EXPORT GSHTransTargets)

  # The public headers are included as <GSHTrans/Core>, so the whole GSHTrans
  # directory goes to the include prefix as it stands. The extensionless
  # umbrella headers have to be named one by one: FILES_MATCHING takes
  # patterns, and there is no extension to match on.
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/GSHTrans
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
          FILES_MATCHING
            PATTERN "*.h"
            PATTERN "All"
            PATTERN "Core"
            PATTERN "Field"
            PATTERN "Tensor"
            PATTERN "Expansion"
            PATTERN "Layered")

  install(EXPORT GSHTransTargets
          FILE GSHTransTargets.cmake
          NAMESPACE GSHTrans::
          DESTINATION "${GSHTRANS_CMAKE_INSTALL_DIR}")

  configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/GSHTransConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/GSHTransConfig.cmake"
    INSTALL_DESTINATION "${GSHTRANS_CMAKE_INSTALL_DIR}")

  # SameMajorVersion rather than AnyNewerVersion: the library is pre-1.0 in
  # spirit and its interfaces have been changing, so promising a consumer that
  # any newer version will do would be a promise it cannot keep.
  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/GSHTransConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
    ARCH_INDEPENDENT)

  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/GSHTransConfig.cmake"
                "${CMAKE_CURRENT_BINARY_DIR}/GSHTransConfigVersion.cmake"
          DESTINATION "${GSHTRANS_CMAKE_INSTALL_DIR}")
endif()

# Developer targets, only when this is the top-level project.
if(PROJECT_IS_TOP_LEVEL)
  option(GSHTRANS_BUILD_EXAMPLES   "Build the examples"   ON)
  option(GSHTRANS_BUILD_TESTS      "Build the test suite" ON)
  option(GSHTRANS_BUILD_BENCHMARKS "Build the benchmarks" ON)

  # Before the subdirectories, because the examples register tests too.
  enable_testing()

  if(GSHTRANS_BUILD_EXAMPLES)
    add_subdirectory(examples)
  endif()
  if(GSHTRANS_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
  endif()
  if(GSHTRANS_BUILD_TESTS)
    add_subdirectory(tests)
  endif()
endif()
