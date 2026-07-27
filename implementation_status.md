# Implementation Status

Branch: `fix/correctness-and-safety`

## Phase 1: Core numerical safety

- [x] Use the configured upper index in `Wigner` single-index convenience
  accessors.
- [x] Give degree-zero grids one longitude with the correct quadrature weight.
- [x] Correct degree-zero forward and inverse normalization.
- [x] Add focused Wigner and degree-zero regression tests.
- [x] Cover the valid upper-index boundaries `n = ±lMax`.
- [x] Build and run the Phase 1 test suite.

### Verification

- Clean Debug build completed in `/tmp/gshtrans-phase1-mwIqx2`.
- All 14 tests passed, including five repeated runs of every test.
- The six focused regressions passed with AddressSanitizer and UBSan in
  `/tmp/gshtrans-phase1-asan-KpyRbq`; leak detection was disabled because it is
  unsupported in the ptraced execution environment.
- `git diff --check` passed.

## Phase 2: Canonical field correctness

- [x] Use indexed access for every compound field operation.
- [x] Materialize lazy expressions through the field's grid-only constructor.
- [x] Own unary and scalar callable state by value and initialize it with
  forwarding.
- [x] Delete the unusable field default constructor.
- [x] Make copy and move assignment copy values without rebinding the
  destination grid.
- [x] Test real and complex arithmetic, unary and scalar expressions,
  materialization, callable ownership, and assignment.

### Verification

- The Debug build completed in `/tmp/gshtrans-phase1-mwIqx2`, including all
  examples and test executables.
- All 18 tests passed, including five repeated runs of every test.
- All 18 tests passed with AddressSanitizer and UBSan in
  `/tmp/gshtrans-phase1-asan-KpyRbq`; leak detection was disabled because it is
  unsupported in the ptraced execution environment.
- The four focused canonical-field regressions passed in both builds.
- `git diff --check` passed.

## Phase 3: Canonical expansion repair

- [x] Replace public `GSHIndices` inheritance with a composed index member.
- [x] Use complex coefficient storage and `Scalar` for real and complex
  expansions.
- [x] Expose only nonnegative orders for `RealValued` and every order for
  `ComplexValued`.
- [x] Repair degree/order iteration, coefficient indexing, and size queries.
- [x] Provide mutable and const coefficient data views.
- [x] Repair copy/move assignment, addition/subtraction assignment, and scalar
  multiplication/division while preserving the destination grid.
- [x] Restore the canonical component and scalar expansion aliases.
- [x] Remove obsolete commented implementations from the two repaired
  expansion headers.
- [x] Test real and complex traits, indexing, iteration, storage, assignment,
  and compound arithmetic.

### Verification

- The Debug build completed in `/tmp/gshtrans-phase1-mwIqx2`, including all
  examples and test executables.
- All 22 tests passed, including five repeated runs of every test.
- All 22 tests passed with AddressSanitizer and UBSan in
  `/tmp/gshtrans-phase1-asan-KpyRbq`; leak detection was disabled because it is
  unsupported in the ptraced execution environment.
- The four focused canonical-expansion regressions passed in both builds.
- `git diff --check` passed.

## Phase 4: Real-field symmetry validation

- [x] Document `RealValued` as a reduced `m >= 0` representation of real grid
  samples whose stored coefficients remain complex.
- [x] Verify reduced coefficients at `n = ±2` agree with the nonduplicated
  nonnegative-order portion of full complex transforms of the same samples.
- [x] Verify the `m = 0` and Nyquist coefficients satisfy the required real
  Fourier constraints.
- [x] Treat the Nyquist corner explicitly: the full complex transform stores
  the shared FFT bin at negative order and zeros the duplicate positive-order
  corner.
- [x] Verify
  `a(l,-m,-n) = (-1)^(m-n) conjugate(a(l,m,n))` between full transforms at
  `+n` and `-n`, with the duplicated Nyquist corner checked separately.
- [x] Verify reduced inverse synthesis equals synthesis from the explicit
  positive-`n` and negative-`n` Hermitian pair.

### Verification

- All three focused symmetry regressions passed in Debug and with
  AddressSanitizer/UBSan.

## Phase 5: Integrated verification

- [x] Configure a fresh Debug build in `/tmp`.
- [x] Build all examples and test executables.
- [x] Run the complete test suite repeatedly.
- [x] Run focused AddressSanitizer/UBSan probes for Wigner access, degree zero,
  field expressions, expansions, and real-field symmetry.
- [x] Check whitespace and confirm the changed-file scope.

### Verification

- Fresh Debug configure and build completed in
  `/tmp/gshtrans-phase5-debug-Z6EAn1`, reusing only previously downloaded
  dependency source trees.
- All 26 tests passed, including five repeated runs of every test.
- All 26 tests passed with AddressSanitizer and UBSan in
  `/tmp/gshtrans-phase1-asan-KpyRbq`; leak detection was disabled because it is
  unsupported in the ptraced execution environment.
- `git diff --check` passed.
- Only the approved source headers, focused tests, test registration, and this
  status file changed.

## Pull request review follow-up

- [x] Add the direct `<type_traits>` dependency used by callable decay.
- [x] Decay lvalue callable types during class template argument deduction and
  forward them into expression-owned state.
- [x] Add a regression proving unary and scalar expressions copy lvalue
  callables rather than retaining references.
- [x] Remove assumptions about moved-from values from field and expansion
  assignment tests.

### Verification

- All 26 tests passed, including five repeated Debug runs of every test.
- All 26 tests passed with AddressSanitizer and UBSan; leak detection was
  disabled because it is unsupported in the ptraced execution environment.
- The new lvalue-callable regression passed in both builds.

All approved phases are complete.
