# NUMERL

Numerical Erlang is a small API for matrix operations in Erlang, done in NIFs.

# INSTALLATION

Currently, this project is done trought a makefile. Install the following packages:
'''
    sudo apt-get install gcc erlang erlang-eunit liblapacke-dev libgslcblas0 
'''

From there, run 'make' to build the project.


# Matrix creation

Matrix creation is done as such:

    > L0 = numerl:list_to_matrix([[1, 2.0], [3, 4.0]]).
    > Eye = numerl:eye(2).
    > Zeros = numerl:zeros(2).

The eye and zero functions take as argument a positive number; list\_to\_matrix takes a list of list of numbers (floats or ints).
    - get, row, col
    - operator ==, +, -, *

# Operators

The following operators are implemented: comparison, addition, and multiplication.
    > NotTrue = numerl:'=='(L0, Eye).
    > Eye2 = numerl:'+'(Eye, Zeros).
    > Neye = numerl:'-'(Zeros, Eye).
    > Mult = numerl:'*'(L0, L0).

# BLAS

Support for the following BLAS function was added:
 -ddot
 -daxpy
 -dgemv
 -dgemm
 -dgesv
