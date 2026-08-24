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
  17-tangential-tensors
  18-radial-operators
  19-layered-models
  20-transform-kernels
  21-interpolation
  22-wigner-3j
)

foreach(example IN LISTS GSHTRANS_EXAMPLES)
  add_executable(${example} ${example}.cpp)
  target_link_libraries(${example} PRIVATE GSHTrans)
  # Registered as tests, so that they are run rather than merely compiled.
  #
  # They were compile-only until now, and that hid a real breakage: every
  # example ended with FFTWpp::CleanUp(), which upstream has since made throw
  # when any plan is still alive -- so several of them aborted, and the build
  # stayed green. An example nobody runs documents nothing.
  #
  # They are not assertions and mostly just print, so this catches crashes,
  # throws and non-zero exits rather than wrong numbers. Where an example
  # states an identity it checks it itself.
  add_test(NAME example.${example} COMMAND ${example})
endforeach()
