=====
adjac
=====

Automatic Differentiation (forward-mode) for generating sparse
Jacobians, using Fortran 95 and operator overloading.

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
However, adjac_complexan does not support operations that break complex analyticity.

For more information about automatic differentiation, and other AD
software (there are many and adjac does not do anything unusual), see
http://autodiff.org/


Example
-------

Adjac enables computation of the Jacobian of a multivariate function,
requiring only slightly modified code computing the *value* of the
function.

For example, consider the following::

    subroutine my_func(x, y)
        implicit none
        double complex, dimension(3), intent(in) :: x
        double complex, dimension(2), intent(out) :: y

        integer :: j
        do j = 1, 2
            y(j) = log(x(j) / (1 + cos(x(j+1))**2))
        end do
    end subroutine my_func

The following function calculates the same as the above, and in
addition the partial derivatives with respect to `x`::

    subroutine my_func_jac(x_value, y_value, dy_dx)
        use adjac
        implicit none
        double complex, dimension(3), intent(in) :: x_value
        double complex, dimension(2), intent(out) :: y_value
        double complex, dimension(2,3), intent(out) :: dy_dx

	type(adjac_complexan), dimension(3) :: x
	type(adjac_complexan), dimension(2) :: y
        integer :: j

	call adjac_set_independent(x, x_value)

        do j = 1, 2
            y(j) = log(x(j) / (1 + cos(x(j+1))**2))
        end do

	call adjac_get_value(y, y_value)
	call adjac_get_dense_jacobian(y, dy_dx)
    end subroutine my_func_jac

Note that the computational part of the code is unchanged. In general,
only data type replacements of the form ``double precision ->
adjac_double`` are usually necessary to make things work.

See ``examples/*.f95`` for mode a usage examples.

