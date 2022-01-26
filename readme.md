# NUMERL

NumErl is a small API for matrix operations in Erlang.

# Usage

This project should be used as a rebar3 dependency:

```erlang
    {deps, [{numerl, {git, "https://github.com/tanguyl/numerl.git", "master"}}]}.
```

# Installation
Assuming a working Erlang + rebar3 installation, openblas and lapacke are required.  

Ubuntu-like os'es:
```sh
    sudo apt-get install libopenblas-dev liblapacke-dev
```

macOS: lapacke is included trough the accelerate framework; only openblas is required.
```sh
    brew install openblas
```

Windows:
Isn't tested natively, but can be installed via WSL.

# MAN
## API
The following functions are available:

```erlang
Mat_builds = [matrix/1, rnd_matrix/1, eye/1, zeros/2].
Comparator = [equals/2].
Accesors   = [mtfli/1, mtfl/1, get/3, at/2, row/2, col/2].
Misc_ops   = [transpose/1].
BLAPACK    = [inv/1, nrm2/1, vec_dot/2, dot/2].
```

## Matrices
Matrices are defined as follows:

```erlang
-record(matrix,{n_rows, n_cols, bin}).
```
The bin field is a Binary, containing the represented matrices values stored as doubles.

## Matrix creation

The following functions can be used to create matrices:

```erlang
Mat_builds = [matrix/1, rnd_matrix/1, eye/1, zeros/2].
```

They can be used as follows:

```erlang                           
% FCT                              INPUT(S)
numerl:matrix([[1, 2.0],[3, 4.0]]).% L,lists of lists of numberss.
numerl:rnd_matrix(2).              % N,  create random   matrix of size NxN.
numerl:eye(2).                     % N,  create identity matrix of size NxN.
numerl:zeros(2, 2).                % N,M create empty    matrix of size NxM.
```

## Comparison

Matrices are stored as arrays of doubles, the granularity of which makes the =/2 operator unadviced. Instead, numerl provides it's own comparison operator:

```erlang
Comparator = [equals/2].
```

It can be used as follows:

```erlang                                    
E = numerl:eye(2),                 %
Z = nmerl:zeros(2,2),              %
%FCT                               INPUTS
Boolean = numerl:equals(E,Z).      % M1,M2: compared matrices.
```

## Accessors
Matrices content can be extracted with the following functions:
```erlang
Accesors = [mtfli/1, mtfl/1, row/2, col/2, get/3, at/2].
```

```erlang                             
M = numerl:rnd_matrix(5),          %
N = 1,                             %
O = 1,                             %
%FCT                               OUTPUT
numerl:mtfli(M),                   % M as a flattened list of ints
numerl:mtfl(M),                    % M as a flattened list of doubles
numerl:row(M,N),                   % The N'th row of M as a matrix.
numerl:col(M,N),                   % The N'th col of M as a matrix.
numerl:get(N,O,M).                 % The element at position N,O of M.
numerl:at(M,N).                    % The N'th element of the matrix.
```

## Element-wise operations
The following operations can be done element-wise on matrices:

```erlang
Element_wise_ops = [add/2, sub/2 ,mult/2, divide/2].
```

They have the following structure:

```erlang
%FCT                               INPUTS
op(Lval, Rval).                    % Lval: a matrix, Rval: a matrix || a number
```

In case of ```erlang divide/2 ``` operator, for an ```erlang Rval ``` either null or containing a null value, a badarg is thrown.   

These ops can be combined with the ```erlang eval/1 ``` function:

```erlang
M = numerl:rnd_matrix(N),
%FCT                               INPUTS
eval([M, add, 1, div, 2, mult, M]).% L: a list of V1,Op1,V2,OP2...
```

## MISC
```erlang
Misc_ops = [transpose/1].
```

## BLAS and LAPACKE
```erlang
BLAPACK = [inv/1, nrm2/1, vec_dot/2, dot/2].
```