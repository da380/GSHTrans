# This is a rough list of things to look into:

- tangential tensor fields. Should be easy to implement and so specialise operations. 
- 3D tensor fields. Built from a fixed angular grid and then a specified set of radii. These are just "dumb" storage objects in the library, but allow for things like transformations, and pointwise and tangential operations. Thought needs to be given here a bit as the aim would be to interact with, say, a finite-difference or finite-element discretisation in the radial direction that provides a means for doing radial differentiation. 
- Specialisations of the tensor classes for common objects. VectorFields, SymmetricSecondOrderTensorFields, things like that. 