# Not run by ctest: these are measurements, not assertions.
add_executable(gshtrans_transform_benchmark TransformBenchmark.cpp)
set_target_properties(gshtrans_transform_benchmark PROPERTIES
                      OUTPUT_NAME TransformBenchmark)
target_link_libraries(gshtrans_transform_benchmark PRIVATE GSHTrans)
