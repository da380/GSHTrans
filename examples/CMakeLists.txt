# The numbered series is meant to be read in order: each one introduces a
# single aspect and assumes the ones before it.
set(GSHTRANS_EXAMPLES
  01-grid
  02-scalar-field
  03-spin-weight
  04-lazy-evaluation
  05-views
  06-transform
  07-batches-and-threads
  08-tensor-fields
  09-tensor-algebra
  10-real-tensors
  11-elasticity
  12-expansions
  13-raising-and-lowering
  14-wigner-functions
  15-surface-gradient
  16-layered-fields
)

foreach(example IN LISTS GSHTRANS_EXAMPLES)
  add_executable(${example} ${example}.cpp)
  target_link_libraries(${example} GSHTrans)
endforeach()

# Standalone Wigner 3-j experiments, unrelated to the library.
add_executable(wig wig.cpp)
target_link_libraries(wig GSHTrans)

add_executable(wigner3j_tests wigner3j_tests.cpp)
