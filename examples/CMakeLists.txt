add_executable(QuadExample QuadExample.cpp)
target_link_libraries(QuadExample
  Eigen3::Eigen GaussQuad)

add_executable(FFTWExample FFTWExample.cpp)
target_link_libraries(FFTWExample FFTWpp)

