# Not run by ctest: these are measurements, not assertions.
add_executable(TransformBenchmark TransformBenchmark.cpp)
target_link_libraries(TransformBenchmark PRIVATE GSHTrans)
