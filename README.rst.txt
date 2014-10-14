=====
adjac
=====

Automatic Differentiation (forward-mode) for generating sparse Jacobians.

Provides three AD data types:

- adjac_double: double precision AD variable
- adjac_complex: double complex AD variable
- adjac_complexan: double complex analytic AD variable

and the support routines:

- adjac_set_independent: initialize independent variable (dy_i/dx_j = delta_ij)
- adjac_get_value: get values from a dependent variables
- adjac_get_nnz: get number of nonzero entries in the Jacobian
- adjac_get_csr_jacobian: get Jacobian in Compressed Sparse Row format
- adjac_get_coo_jacobian: get Jacobian in sparse coordinate format
- adjac_get_dense_jacobian: get Jacobian as a full matrix

The complex analytic adjac_complexan generates complex-valued
Jacobians corresponding to the complex derivative, whereas
adjac_complex can be used to generate real-valued Jacobians
corresponding to separate derivatives vs. real and imaginary parts
of the variables. In complex-analytic cases, the results will be
equivalent, but adjac_complexan is more efficient computationally.

The data types support operations =,*,+,-,matmul,exp,sin,cos,log,dble,aimag,conjg.
adjac_complexan does not support operations that break complex analyticity.

See examples/*.f90 for a usage examples.
