# NUMERL

NumErl is a small API for matrix operations in Erlang.

# Installation

This project is built trough a makefile. Install the following packages:

    sudo apt-get install gcc erlang erlang-eunit liblapacke-dev libgslcblas0 


From there, run 'make' to build the project.


# Matrix creation

Matrix are created the following ways:

    L0 = numerl:list_to_matrix([[1, 2.0], [3, 4.0]]).
    Eye = numerl:eye(2).
    Zeros = numerl:zeros(2, 2).

The eye and zero functions take as argument positive numbers; list\_to\_matrix takes a list of list of numbers (floats or ints).

# Operators

The following operators are implemented: comparison, addition, and multiplication.

    NotTrue = numerl:'=='(L0, Eye).
    Eye2 = numerl:'+'(Eye, Zeros).
    Neye = numerl:'-'(Zeros, Eye).
    Mult = numerl:'*'(L0, L0).

# Accessors

Access to elements / columns / rows of matrices is done as such:

    One = numerl:get(Eye, 1,1).
    OneO = numerl:row(Eye, 1).
    OOne = numerl:col(Eye, 2).
    
        
The function print can be used to return an atom representation of a matrix.
    
    numerl:print(numerl:eye(2)).
        '[[1.00000 0.00000][0.00000 1.00000]]'

# BLAS

BLAS are stardard, highly optimized functions used to compute operations between matrices and vectors. In numerl, vectors are simply matrices with at least one dimension of size 1:

    X = numerl:matrix([[1.0, 2.0]]).
    Y = numerl:matrix([[3], [4]]).
    
The following variables are used to give examples.

    N = 2.
    Alpha = 1.0.
    Beta = 2.0.
    A = numerl:eye(2).
    B = numerl:eye(2).
    C = numerl:zeros(2,2).


ddot returns the dot produtc of the first N values contained in X and Y, vectors.

    DotRes = numerl:ddot(N,X,Y).

daxpy returns the result of Alpha\*X + Y taking into account the first N coordinates of vectors X and Y, alpha being a number.

    DaxRes = numerl:daxpy(N, Alpha, X, Y).

dgemv returns the result of Alpha\*A\*X + Beta\*Y. Alpha and Beta are numbers; M is a matrix; X and Y are vectors.

    GemvRes = numerl:dgemv(Alpha, A, X, Beta, Y).

dgemm returns the result of Alpha\*A\*B + Beta\*C. Alpha and Beta are numbers; A,B and C are matrices.

    GemmRes = numerl:dgemv(Alpha, A, B, Beta, C).
