# A standalone project that consumes an installed GSHTrans through
# find_package. It exists so that CI notices if the install or export rules
# stop producing a usable package -- which is a failure mode the main build
# cannot see, since there the targets are all present as real ones. It is
# deliberately not part of the main build.
cmake_minimum_required(VERSION 3.20)
project(GSHTransPackageTest LANGUAGES CXX)

find_package(GSHTrans REQUIRED)

add_executable(PackageTest main.cpp)
target_link_libraries(PackageTest PRIVATE GSHTrans::GSHTrans)

enable_testing()
add_test(NAME PackageTest COMMAND PackageTest)
