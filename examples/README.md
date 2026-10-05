# Examples

A tutorial series, meant to be read in order: each example introduces one
aspect of the library and assumes the ones before it. Each is a short,
self-contained `main` whose header comment says what it demonstrates, what to
read first, which types and functions it introduces, and what its output
shows. All of them print something that shows the point rather than only
asserting it, and where an example states an identity it checks it and exits
non-zero if it fails.

The reader assumed is one who knows C++ and some spherical harmonics, but not
this library. The mathematics behind each step — conventions, normalisation,
the operators — is in `docs/gshtrans-reference.tex`, which the examples cite
by section name. The top-level `README.md` is the overview.

## Building and running

The examples are built with the library (`GSHTRANS_BUILD_EXAMPLES`, on by
default) and registered as tests, so `ctest` runs every one of them:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
cd build && ctest -R example
```

Each executable is named after its source file and lands in `build/bin`, so
`build/bin/12-expansions` runs one directly. Example 14 writes a data file
(`wigner-l6-m2.dat`) to the current directory, and a plot if `gnuplot` is on
the path.

A few examples have parts that depend on optional components. Those parts are
compiled only when the component is present, and the example says so when it
is not:

* the `Interpolation` dependency (`GSHTRANS_WITH_INTERPOLATION`): the spline
  derivative in 18, splines and resampling in 19, and the bilinear and
  bicubic schemes in 21;
* a BLAS (`GSHTRANS_WITH_BLAS`): the matrix transform kernel, which is the
  whole of 20.

## The series

### Fields on the sphere

| | |
|---|---|
| `01-grid` | the grid, `ForBand` and oversampling, quadrature weights, value semantics |
| `02-scalar-field` | a field from a function, and `Integrate` |
| `03-spin-weight` | the upper index (spin weight), why `conj` reverses it, and what will not compile |
| `04-lazy-evaluation` | expressions, `Materialise`, evaluation in place, `Map` |
| `05-views` | a field over storage the library did not allocate, contiguous or strided |
| `06-transform` | forward and inverse, and why transforming an arbitrary field is a *projection* |
| `07-batches-and-threads` | many fields at once, the batch descriptor, chunking, the explicit threading policy |

### Tensor fields

| | |
|---|---|
| `08-tensor-fields` | components by multi-index, upper indices, symmetry as storage |
| `09-tensor-algebra` | trace against the metric, transpose, symmetrise, materialise |
| `10-real-tensors` | the reality reduction, and the storage counts it achieves |
| `11-elasticity` | an elastic tensor applied to a strain, the rank-4 case end to end |

### The spectral side and derivatives

| | |
|---|---|
| `12-expansions` | coefficients indexed by degree and order, for spin fields and tensors |
| `13-raising-and-lowering` | `ð` and `ð̄`, the surface Laplacian, and the commutator identity |
| `14-wigner-functions` | the `d`-functions themselves: access, the stored normalisation, the relations that pin it, and a plot |
| `15-surface-gradient` | the contravariant derivative of Phinney & Burridge and Dahlen & Tromp: what it is, and why it is not `ð` applied to each component |

### Three-dimensional fields

| | |
|---|---|
| `16-layered-fields` | fields on a ball or shell: slices, the radial axis as the batch axis, radial operators, and the full gradient |
| `17-tangential-tensors` | the bundle with no radial slot: what it saves, `Embed` and `Tangential`, and the intrinsic derivative, which is closed on it |
| `18-radial-operators` | the ready-made radial derivatives, the rules for writing your own, and working on contiguous radial lines |
| `19-layered-models` | elements and interfaces: a two-valued derivative at a discontinuity, and resampling that does not cross one |

### Numerical tools

| | |
|---|---|
| `20-transform-kernels` | the two Legendre kernels, why both are kept, what the matrix kernel refuses, and BLAS threading |
| `21-interpolation` | a field as a callable of the two angles: the exact scheme, the local ones, the poles and the longitude wrap, and remeshing in one line |
| `22-wigner-3j` | coupling coefficients: the table, the stack, two closed forms, stretched triangles, and the independent check that says the values are right |

Example 22 stands apart from the rest: the 3-j symbols touch no grid,
transform or field, and it can be read at any point.
