add_library(GSHT INTERFACE
  Legendre.h
)

target_include_directories(GSHT INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/GSHT>
  )
