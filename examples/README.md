# Examples

A series, meant to be read in order: each one introduces a single aspect and
assumes the ones before it. All are short, and all print something that shows
the point rather than only asserting it.

| | |
|---|---|
| `01-grid` | the grid, `ForBand` and oversampling, quadrature weights, value semantics |
| `02-scalar-field` | a field from a function, and `Integrate` |
| `03-spin-weight` | why `conj` reverses the upper index, and what will not compile |
| `04-lazy-evaluation` | expressions, `Materialise`, evaluation in place, `Map` |
| `05-views` | a field over storage the library did not allocate, contiguous or strided |
| `06-transform` | forward and inverse, and why transforming an arbitrary field is a *projection* |
| `07-batches-and-threads` | many fields at once, the batch descriptor, chunking, the explicit threading policy |
| `08-tensor-fields` | components by multi-index, upper indices, symmetry as storage |
| `09-tensor-algebra` | trace against the metric, transpose, symmetrise, materialise |
| `10-real-tensors` | the reality reduction, and the storage counts it achieves |
| `11-elasticity` | an elastic tensor applied to a strain, the rank-4 case end to end |
| `12-expansions` | the spectral side, indexed by degree and order |
| `13-raising-and-lowering` | `ð` and `ð̄`, and the surface Laplacian |
| `14-wigner-functions` | the `d`-functions themselves: access, the value convention, the relations that pin it, and a plot |
| `15-surface-gradient` | the contravariant derivative of P&B: what it is, and why it is not `ð` |

`wig` and `wigner3j_tests` are standalone Wigner 3-j experiments and are not
part of the series.

The mathematics and the numerical methods are documented separately, in
`docs/`.
