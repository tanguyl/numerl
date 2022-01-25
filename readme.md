# NUMERL

NumErl is a small API for matrix operations in Erlang.

# Usage

This project should be used as a rebar3 dependency:

```erlang
    {deps, [{numerl, {git, "https://github.com/tanguyl/numerl.git", "master"}}]}.
```

# Installation
Assuming you have a working Erlang + rebar3 installation, you will need openblas and lapacke:

Ubuntu-like linux:
```sh
    sudo apt-get install libopenblas-dev liblapacke-dev
```

Mac: lapacke should be usable via the accelerate framework; only openblas is required.
```sh
    brew install openblas
```

Windows:
Isn't supported 'out of the box'; but can use Numerl trough WSL.

# API

## Matrix creation

Matrices are created the following ways:

```erlang
L0 = numerl:matrix([[1, 2.0], [3, 4.0]]).
Eye = numerl:eye(2).
Zeros = numerl:zeros(2, 2).
```

The eye and zero functions take as argument positive numbers; list\_to\_matrix takes a list of list of numbers (floats or ints).

## Operators

The following operators are implemented: comparison, addition, and multiplication.

```erlang
NotTrue = numerl:equals(L0, Eye).
Eye2    = numerl:add(Eye, Zeros).
Neye    = numerl:sub(Zeros, Eye).
Mult    = numerl:mult(L0, L0). %Arg2 either matrix (elemt-wise op) or number
D       = numerl:div(L0, 2).   %Arg2 number.
Inv     = numerl:inv(L0).       %Arg2 
Tr      = numerl:transpose(L0).
```

## Accessors

Access to elements / columns / rows of matrices is done as such:

```erlang
One     = numerl:get(Eye, 1,1).
OneO    = numerl:row(Eye, 1).
OOne    = numerl:col(Eye, 2).
[1]     = numerl:mtfli(numerl:matrix([[1.0]])). %Matrix to flattened list of int.
1.0     = numerl:at(Eye,1). %Returns N'th value in flattened array.
```
        
The function print can be used to return an atom representation of a matrix.

```erlang
1> numerl:print(numerl:eye(2)).
'[[1.00000 0.00000][0.00000 1.00000]]'
```
## BLAS

BLAS are stardard, highly optimized functions used to compute operations between matrices and vectors. In numerl, vectors are simply matrices with at least one dimension of size 1:

```erlang
X = numerl:matrix([[1.0, 2.0]]).
Y = numerl:matrix([[3], [4]]).
```
    
The following variables are used to give examples.

```erlang
N = 2.
Alpha = 1.0.
Beta = 2.0.
A = numerl:eye(2).
B = numerl:eye(2).
C = numerl:zeros(2,2).
```

ddot returns the dot produtc of the first N values contained in X and Y, vectors.

```erlang
DotRes = numerl:ddot(N,X,Y).
```

daxpy returns the result of Alpha\*X + Y taking into account the first N coordinates of vectors X and Y, alpha being a number.

```erlang
DaxRes = numerl:daxpy(N, Alpha, X, Y).
```

dgemv returns the result of Alpha\*A\*X + Beta\*Y. Alpha and Beta are numbers; M is a matrix; X and Y are vectors.

```erlang
GemvRes = numerl:dgemv(Alpha, A, X, Beta, Y).
```

dgemm returns the dot product of A*B.

```erlang
GemmRes = numerl:dot(A, B).
```

Support for LAPACKE (dgesv) was temporarily dropped for easier OSX integration.
