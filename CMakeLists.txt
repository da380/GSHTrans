# 3.24, not the 3.20 this declared until now. Every dependency below is
# declared with FetchContent_Declare(... FIND_PACKAGE_ARGS), which is a 3.24
# feature, so a 3.20 configure never worked -- it simply failed later and less
# clearly. Raising the floor does not break anyone who was working, because
# nobody between 3.20 and 3.23 was.
cmake_minimum_required(VERSION 3.24)

project(GSHTrans
        VERSION 1.2.0
        DESCRIPTION "Generalized spherical harmonic transforms and a tensor field algebra on the sphere"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# The version in GSHTrans/Version.hpp is written by hand, so that copying
# the GSHTrans directory into another project is a complete way to vendor the
# library. Check it against the version declared above, so that the two cannot
# drift apart.
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/GSHTrans/Version.hpp _version_header)
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
    "GSHTrans/Version.hpp declares ${_header_version}. Update both, "
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
#
# Which revision of each is fetched is a cache variable, and the default is
# main. The four are one author's and move together, so on a development
# branch tracking main is the point. It is wrong for a *release*: a tag that
# fetches whatever main is that day does not name a state anyone can rebuild,
# and an upstream change has already broken a tagged tree once. So the commit
# that tags a release sets these four defaults to the commits it was tested
# against, and the commit after it sets them back. They are cache variables
# rather than literals so that the same thing can be done from outside, by
# anyone who needs a build to stay still:
#
#   cmake -DGSHTRANS_FFTWPP_TAG=<commit> ...
set(GSHTRANS_NUMERICCONCEPTS_TAG main CACHE STRING
    "The NumericConcepts revision to fetch if none is installed")
set(GSHTRANS_GAUSSQUAD_TAG main CACHE STRING
    "The GaussQuad revision to fetch if none is installed")
set(GSHTRANS_FFTWPP_TAG main CACHE STRING
    "The FFTWpp revision to fetch if none is installed")
set(GSHTRANS_INTERPOLATION_TAG main CACHE STRING
    "The Interpolation revision to fetch if none is installed")
mark_as_advanced(GSHTRANS_NUMERICCONCEPTS_TAG GSHTRANS_GAUSSQUAD_TAG
                 GSHTRANS_FFTWPP_TAG GSHTRANS_INTERPOLATION_TAG)
FetchContent_Declare(
  NumericConcepts
  GIT_REPOSITORY https://github.com/da380/NumericConcepts.git
  GIT_TAG ${GSHTRANS_NUMERICCONCEPTS_TAG}
  FIND_PACKAGE_ARGS)

FetchContent_Declare(
  GaussQuad
  GIT_REPOSITORY https://github.com/da380/GaussQuad.git
  GIT_TAG ${GSHTRANS_GAUSSQUAD_TAG}
  FIND_PACKAGE_ARGS)

FetchContent_Declare(
  FFTWpp
  GIT_REPOSITORY https://github.com/da380/FFTWpp.git
  GIT_TAG ${GSHTRANS_FFTWPP_TAG}
  FIND_PACKAGE_ARGS)

FetchContent_MakeAvailable(NumericConcepts GaussQuad FFTWpp)

# Interpolation, which is optional and is the one dependency only part of the
# library uses.
#
# Nothing angular wants a one-dimensional interpolator: the transform, the
# index algebra and the surface gradient are complete without it, and so are
# all three ready-made radial derivatives, which are self-contained by
# design. What it buys is
# resampling a layered field onto a different radial grid, where the whole
# menu of schemes comes for the price of a policy argument.
#
# On by default because it is header-only, depends on nothing but
# NumericConcepts -- which is already here -- and is small. Off is a supported
# configuration and the suite runs in it.
#
# Tracking main by default, like the three above. Interpolation's own README
# asks consumers to pin a reviewed commit, which is right for a stranger and
# wrong for four libraries with one author that move together -- except at a
# release, which is what GSHTRANS_INTERPOLATION_TAG is for.
option(GSHTRANS_WITH_INTERPOLATION
       "Radial resampling, built on da380/Interpolation" ON)

if(GSHTRANS_WITH_INTERPOLATION)
  FetchContent_Declare(
    Interpolation
    GIT_REPOSITORY https://github.com/da380/Interpolation.git
    GIT_TAG ${GSHTRANS_INTERPOLATION_TAG}
    FIND_PACKAGE_ARGS)
  FetchContent_MakeAvailable(Interpolation)
endif()

# The GEMM behind the matrix transform kernel.
#
# Found and not fetched, unlike the four above, because BLAS is neither
# header-only nor ours: a system OpenBLAS, MKL or Accelerate is what a caller
# on a real machine already has, and building one here would be this library
# taking on a responsibility it has declined for every other dependency. It is
# therefore the first non-header-only dependency, which is why it is an option
# and not a requirement -- everything the library did before this exists
# without it, and a build with it off simply does not offer
# TransformKernel::Matrix().
#
# Three states rather than two, because ON has to mean two different things
# depending on who said it. AUTO -- the default -- uses BLAS if it is there
# and says so if it is not. ON demands it, and fails the configure if it is
# missing, which is what a caller who typed it deserves rather than a silent
# downgrade to the loop kernel. OFF never looks.
#
# A note for whoever links this against a threaded BLAS: the matrix kernel
# threads over orders itself, so the BLAS must be single-threaded when it
# does. A BLAS on the same OpenMP runtime is held to one thread by the library
# itself; one with its own pthread pool is out of its reach and needs
# OPENBLAS_NUM_THREADS=1 or its equivalent.
set(GSHTRANS_WITH_BLAS AUTO CACHE STRING
    "Matrix transform kernel, built on a system BLAS (AUTO, ON, OFF)")
set_property(CACHE GSHTRANS_WITH_BLAS PROPERTY STRINGS AUTO ON OFF)

# The three states are compared as strings, so the other spellings CMake
# accepts for a boolean have to be folded onto them first. Unfolded, "0" and
# "FALSE" were neither "OFF" nor "ON" and so meant AUTO -- BLAS linked for a
# caller who had asked for none -- and "1" and "TRUE" lost the hard failure
# that ON exists to give. Anything else is refused rather than guessed at.
string(TOUPPER "${GSHTRANS_WITH_BLAS}" _with_blas)
if(_with_blas MATCHES "^(ON|1|TRUE|YES|Y)$")
  set(_with_blas ON)
elseif(_with_blas MATCHES "^(OFF|0|FALSE|NO|N)$")
  set(_with_blas OFF)
elseif(NOT _with_blas STREQUAL "AUTO")
  message(FATAL_ERROR
    "GSHTRANS_WITH_BLAS is \"${GSHTRANS_WITH_BLAS}\", "
    "which is none of AUTO, ON and OFF.")
endif()

set(GSHTRANS_BLAS_FOUND OFF)
if(NOT _with_blas STREQUAL "OFF")
  if(_with_blas STREQUAL "ON")
    find_package(BLAS REQUIRED)
  else()
    find_package(BLAS)
  endif()
  if(BLAS_FOUND)
    set(GSHTRANS_BLAS_FOUND ON)
  else()
    message(STATUS
      "GSHTrans: no BLAS found, so the matrix transform kernel is not built. "
      "The loop kernel is unaffected and is what the library has always used. "
      "Configure with -DGSHTRANS_WITH_BLAS=ON to make this an error instead.")
  endif()
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

# OpenMP, which the library threads with and does not need.
#
# On by default, and then required: a build that asked for threads and did not
# get them should say so at configure time and not by running on one core.
# Off, every parallel region is a serial loop and every policy that asks for
# threads gets one, with the same answers. That is a supported configuration,
# and CI runs the suite in it: it is what makes a bare copy of the GSHTrans
# directory usable with nothing but a compiler and FFTW, it builds where
# OpenMP is awkward to come by, and it is the build to take to a debugger.
#
# Nothing is defined for the headers either way. They ask the compiler, through
# _OPENMP, which is set exactly when the flag that enables OpenMP is given --
# so there is no macro of ours to get out of step with it.
option(GSHTRANS_WITH_OPENMP "Thread with OpenMP" ON)
if(GSHTRANS_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

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
  GaussQuad::GaussQuad)

if(GSHTRANS_WITH_OPENMP)
  target_link_libraries(GSHTrans INTERFACE OpenMP::OpenMP_CXX)
endif()

# The optional half, and the macro that guards the headers using it. Carried
# on the target rather than generated into a header, because Version.hpp is
# written by hand so that copying the GSHTrans directory remains a complete
# way to vendor the library -- and a build-tree-generated header would break
# exactly that. An INTERFACE definition is exported with the target, so an
# installed consumer sees the same answer this build did.
if(GSHTRANS_WITH_INTERPOLATION)
  target_link_libraries(GSHTrans INTERFACE Interpolation::Interpolation)
  target_compile_definitions(GSHTrans INTERFACE GSHTRANS_HAVE_INTERPOLATION)
endif()

if(GSHTRANS_BLAS_FOUND)
  target_link_libraries(GSHTrans INTERFACE BLAS::BLAS)
  target_compile_definitions(GSHTrans INTERFACE GSHTRANS_HAVE_BLAS)
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
# which is what unblocked this.
#
# Consumption by add_subdirectory or FetchContent still works and is
# unaffected.
if(GSHTRANS_INSTALL)
  include(CMakePackageConfigHelpers)
  set(GSHTRANS_CMAKE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/GSHTrans")

  install(TARGETS GSHTrans EXPORT GSHTransTargets)

  # The public headers are included as <GSHTrans/GSHTrans.hpp>, so the whole
  # GSHTrans directory goes to the include prefix as it stands. All and Core
  # are the two forwarding headers kept for code written against the old
  # layout; they have no extension, so they are named.
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/GSHTrans
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
          FILES_MATCHING
            PATTERN "*.hpp"
            PATTERN "All"
            PATTERN "Core")

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
  # The benchmark measures threading and calls the OpenMP runtime directly, so
  # without OpenMP there is nothing for it to be.
  if(GSHTRANS_BUILD_BENCHMARKS AND GSHTRANS_WITH_OPENMP)
    add_subdirectory(benchmarks)
  endif()
  if(GSHTRANS_BUILD_TESTS)
    add_subdirectory(tests)
  endif()

  # The API documentation, built on request rather than by default: doxygen is
  # not needed to build or use the library, and a missing one should not fail
  # a configure. `cmake --build build --target docs`.
  find_package(Doxygen QUIET)
  if(Doxygen_FOUND)
    set(DOXYGEN_INPUT_DIR  "${PROJECT_SOURCE_DIR}/GSHTrans")
    set(DOXYGEN_OUTPUT_DIR "${PROJECT_BINARY_DIR}/docs")
    configure_file("${PROJECT_SOURCE_DIR}/docs/Doxyfile.in"
                   "${PROJECT_BINARY_DIR}/Doxyfile" @ONLY)
    add_custom_target(docs
      COMMAND ${CMAKE_COMMAND} -E make_directory "${DOXYGEN_OUTPUT_DIR}"
      COMMAND ${DOXYGEN_EXECUTABLE} "${PROJECT_BINARY_DIR}/Doxyfile"
      WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
      COMMENT "Building the API documentation with Doxygen"
      VERBATIM)
  endif()
endif()
