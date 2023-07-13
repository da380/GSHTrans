add_library (GSHT INTERFACE GSHT.h)

target_link_libraries(GSHT INTERFACE Eigen3::Eigen FFTWpp GaussQuad)

target_include_directories (GaussQuad INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>
)
