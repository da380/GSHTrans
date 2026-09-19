include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  # A release and not main. GoogleTest is someone else's, its main raises its
  # minimum compiler and CMake versions without notice, and nothing here needs
  # anything newer than this.
  GIT_TAG v1.17.0
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
# Do not let a test dependency install itself alongside the library: without
# this, "cmake --install" ships six hundred files of GoogleTest.
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(TestConcepts
               TestConcepts.cpp)
target_link_libraries(TestConcepts PRIVATE GSHTrans gtest_main)

add_executable(TestSpinField
               TestSpinField.cpp)
target_link_libraries(TestSpinField PRIVATE GSHTrans gtest_main)

add_executable(TestTensorIndices
               TestTensorIndices.cpp)
target_link_libraries(TestTensorIndices PRIVATE GSHTrans gtest_main)

add_executable(TestTensorField
               TestTensorField.cpp)
target_link_libraries(TestTensorField PRIVATE GSHTrans gtest_main)

add_executable(TestTensorAlgebra
               TestTensorAlgebra.cpp)
target_link_libraries(TestTensorAlgebra PRIVATE GSHTrans gtest_main)

add_executable(TestTensorReality
               TestTensorReality.cpp)
target_link_libraries(TestTensorReality PRIVATE GSHTrans gtest_main)

add_executable(TestSinglePrecision
               TestSinglePrecision.cpp)
target_link_libraries(TestSinglePrecision PRIVATE GSHTrans gtest_main)

add_executable(TestLargeDegree
               TestLargeDegree.cpp)
target_link_libraries(TestLargeDegree PRIVATE GSHTrans gtest_main)

add_executable(TestPreconditions
               TestPreconditions.cpp)
target_link_libraries(TestPreconditions PRIVATE GSHTrans gtest_main)

add_executable(TestTensorOrbitValues
               TestTensorOrbitValues.cpp)
target_link_libraries(TestTensorOrbitValues PRIVATE GSHTrans gtest_main)

add_executable(TestExpansion
               TestExpansion.cpp)
target_link_libraries(TestExpansion PRIVATE GSHTrans gtest_main)

add_executable(TestInterpolate
               TestInterpolate.cpp)
target_link_libraries(TestInterpolate PRIVATE GSHTrans gtest_main)

add_executable(TestThreeJ
               TestThreeJ.cpp)
target_link_libraries(TestThreeJ PRIVATE GSHTrans gtest_main)

add_executable(TestTuning
               TestTuning.cpp)
target_link_libraries(TestTuning PRIVATE GSHTrans gtest_main)

add_executable(TestConventions
               TestConventions.cpp)
target_link_libraries(TestConventions PRIVATE GSHTrans gtest_main)

add_executable(TestLayered
               TestLayered.cpp)
target_link_libraries(TestLayered PRIVATE GSHTrans gtest_main)

add_executable(TestWigner
  TestWigner.cpp)
target_link_libraries(TestWigner PRIVATE GSHTrans gtest_main)


add_executable(TestGaussLegendreGrid
               TestGaussLegendreGrid.cpp)
target_link_libraries(TestGaussLegendreGrid PRIVATE GSHTrans gtest_main)

add_executable(TestRealFieldSymmetry
               TestRealFieldSymmetry.cpp)
target_link_libraries(TestRealFieldSymmetry PRIVATE GSHTrans gtest_main)

# Every core header must compile on its own, and twice.
# One translation unit per header, generated at configure time.
#
# The list is every header the library ships, plus the umbrella module
# headers, and the check below is what makes that sentence true: the list is
# compared with what is on disk, and a header that is there and not here stops
# the configure. It is written out and not globbed because adding a header
# should be a deliberate act -- but "deliberate" had quietly become
# "forgotten" for four of them, all in Layered/, and the claim above stood in
# the README regardless.
set(GSHTRANS_SELF_SUFFICIENT_HEADERS
  GSHTrans/src/3j.h
  GSHTrans/src/Blas.h
  GSHTrans/src/Concepts.h
  GSHTrans/src/GaussLegendreGrid.h
  GSHTrans/src/SphericalGrid.h
  GSHTrans/src/Indexing.h
  GSHTrans/src/Policies.h
  GSHTrans/src/Tuning.h
  GSHTrans/src/Utility.h
  GSHTrans/src/Version.h
  GSHTrans/src/Views.h
  GSHTrans/src/Wigner.h
  GSHTrans/src/WignerMatrices.h
  GSHTrans/src/SpinField/SpinWeighted.h
  GSHTrans/src/SpinField/SpinField.h
  GSHTrans/src/SpinField/SpinFieldNodes.h
  GSHTrans/src/SpinField/SpinFieldOverloads.h
  GSHTrans/src/SpinField/SpinFieldView.h
  GSHTrans/src/Tensor/MultiIndex.h
  GSHTrans/src/Tensor/Orbits.h
  GSHTrans/src/Tensor/TensorField.h
  GSHTrans/src/Tensor/TensorExpr.h
  GSHTrans/src/Tensor/BundleMaps.h
  GSHTrans/src/Expansion/SpinExpansion.h
  GSHTrans/src/Expansion/Eth.h
  GSHTrans/src/Expansion/TensorExpansion.h
  GSHTrans/src/Expansion/ContravariantDerivative.h
  GSHTrans/src/Expansion/IntrinsicDerivative.h
  GSHTrans/src/Expansion/BundleMaps.h
  GSHTrans/src/Expansion/Interpolate.h
  GSHTrans/src/Layered/RadialGrid.h
  GSHTrans/src/Layered/LayeredSpinField.h
  GSHTrans/src/Layered/LayeredTensorField.h
  GSHTrans/src/Layered/LayeredGradient.h
  GSHTrans/src/Layered/RadialOperator.h
  GSHTrans/src/Layered/RadialMajor.h
  GSHTrans/src/Layered/RadialDerivatives.h
  GSHTrans/src/Layered/RadialSplineDerivative.h
  GSHTrans/src/Layered/RadialResample.h
  GSHTrans/Core
  GSHTrans/Field
  GSHTrans/Tensor
  GSHTrans/Expansion
  GSHTrans/Layered
  GSHTrans/All
)

file(GLOB_RECURSE _shipped_headers
     RELATIVE "${PROJECT_SOURCE_DIR}"
     CONFIGURE_DEPENDS
     "${PROJECT_SOURCE_DIR}/GSHTrans/src/*.h")
foreach(header IN LISTS _shipped_headers)
  if(NOT header IN_LIST GSHTRANS_SELF_SUFFICIENT_HEADERS)
    message(FATAL_ERROR
      "${header} is shipped but is not in GSHTRANS_SELF_SUFFICIENT_HEADERS, so "
      "nothing checks that it compiles on its own. Add it to the list in "
      "tests/CMakeLists.txt.")
  endif()
endforeach()

set(GSHTRANS_SELF_SUFFICIENCY_SOURCES)
foreach(header IN LISTS GSHTRANS_SELF_SUFFICIENT_HEADERS)
  string(REGEX REPLACE "[/.]" "_" stem "${header}")
  set(generated "${CMAKE_CURRENT_BINARY_DIR}/self_sufficiency/${stem}.cpp")
  file(WRITE "${generated}"
    "// Generated by tests/CMakeLists.txt. Do not edit.\n"
    "#include <${header}>\n"
    "#include <${header}>\n")
  list(APPEND GSHTRANS_SELF_SUFFICIENCY_SOURCES "${generated}")
endforeach()

add_library(HeaderSelfSufficiency OBJECT ${GSHTRANS_SELF_SUFFICIENCY_SOURCES})
target_link_libraries(HeaderSelfSufficiency PRIVATE GSHTrans)

# Warnings, on the tests and nowhere else. The library is header-only, so its
# warnings are raised in whatever includes it -- which for a consumer building
# with -Werror is a build they cannot do. The tests include every header, so
# compiling them with the usual three is how the library is kept quiet for
# everyone else. PRIVATE, and on these targets: nothing here reaches the
# GSHTrans interface target, which must not impose flags on anyone.
get_property(_test_targets DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
foreach(target IN LISTS _test_targets)
  if(target MATCHES "^Test" OR target STREQUAL "HeaderSelfSufficiency")
    target_compile_options(${target} PRIVATE
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>)
  endif()
endforeach()

include(GoogleTest)
gtest_discover_tests(TestConcepts)
gtest_discover_tests(TestSpinField)
gtest_discover_tests(TestTensorIndices)
gtest_discover_tests(TestTensorField)
gtest_discover_tests(TestTensorAlgebra)
gtest_discover_tests(TestTensorReality)
gtest_discover_tests(TestTensorOrbitValues)
gtest_discover_tests(TestPreconditions)
gtest_discover_tests(TestLargeDegree)
gtest_discover_tests(TestSinglePrecision)
gtest_discover_tests(TestExpansion)
gtest_discover_tests(TestInterpolate)
gtest_discover_tests(TestTuning)
gtest_discover_tests(TestThreeJ)
gtest_discover_tests(TestConventions)
gtest_discover_tests(TestLayered)
gtest_discover_tests(TestWigner)
gtest_discover_tests(TestGaussLegendreGrid)
gtest_discover_tests(TestRealFieldSymmetry)
