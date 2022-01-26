# NUMERL

NUMeric ERLang is an API for matrix operations in Erlang. All functions (except ```[eval/1, rnd_matrix/1]```) are implemented in NIFs; some make use of BLAS and LAPACKE ```[inv/1, nrm2/1, vec_dot/2, dot/2]``` for even better performances.  

## 1ms rule
This library does not yet take into account the 1ms rule of NIFs; DO NOT USE BIG MATRICES yet.   
Work is in progress to make ```numerl```'s functions execute on chunks of matrices in parallel, solving this issue and improving performances.   
This will be solved by the end of 2022's Q1.

# Usage

This library is a rebar3 dependency (https://github.com/erlang/rebar3):

```erlang
    {deps, [{numerl, {git, "https://github.com/tanguyl/numerl.git", {branch, "master"}}}]}.
```

A ```numerl``` example can be found at https://github.com/tanguyl/raytracer.

# Installation
Openblas and lapacke are mandatory to use this library. 

Ubuntu-like os'es:
```sh
    sudo apt-get install libopenblas-dev liblapacke-dev
```

macOS:   
lapacke is included trough the accelerate framework; only openblas is required.
```sh
    brew install openblas
```

Windows:   
Tested on WSL :)

# MAN
## API
The following functions are available:

```erlang
<<<<<<< HEAD
Mat_builds = [matrix/1, rnd_matrix/1, eye/1, zeros/2].
Comparator = [equals/2].
Accesors   = [mtfli/1, mtfl/1, get/3, at/2, row/2, col/2].
BLAPACK    = [inv/1, nrm2/1, vec_dot/2, dot/2, transpose/1].
=======
Constructors     = [matrix/1, rnd_matrix/1, eye/1, zeros/2].
Comparator       = [equals/2].
Accesors         = [mtfli/1, mtfl/1, get/3, at/2, row/2, col/2].
Element_wise_ops = [add/2, sub/2 ,mult/2, divide/2].
Ops              = [inv/1, nrm2/1, vec_dot/2, dot/2, transpose/1].
>>>>>>> 199ce279a4b4e7b6d98add6ddfd096c6706475cd
```

## Matrices
Matrices are defined as follows:

```erlang
-record(matrix,{n_rows, n_cols, bin}).
```
The bin field is a Binary, containing the represented matrices values stored as doubles.

## Constructors

The following functions can be used to create matrices:

```erlang
Constructors = [matrix/1, rnd_matrix/1, eye/1, zeros/2].
```

They can be used as follows:

```erlang                           
% FCT                              INPUT(S)   OUTPUT    
<<<<<<< HEAD
numerl:matrix([[1, 2.0],[3, 4.0]]).% L
=======
numerl:matrix([[1, 2.0],[3, 4.0]]).% L    
>>>>>>> 199ce279a4b4e7b6d98add6ddfd096c6706475cd
numerl:rnd_matrix(2).              % N,     random   matrix of size NxN.
numerl:eye(2).                     % N,     identity matrix of size NxN.
numerl:zeros(2, 2).                % N,M    empty    matrix of size NxM.
```

The ```matrix/1``` takes as input a list of rows.

<<<<<<< HEAD
## Comparison

Matrices are stored as arrays of doubles, the granularity of which makes the ```=/2``` operator unadviced. Instead, numerl provides its own comparison operator:
=======
## Comparator

Matrices are stored as arrays of doubles, the granularity of which makes the ```==/2``` operator unadviced. Instead, numerl provides its own comparison operator:
>>>>>>> 199ce279a4b4e7b6d98add6ddfd096c6706475cd

```erlang
Comparator = [equals/2].
```

It should be used as follows:

```erlang                                    
E = numerl:eye(2),                 %
<<<<<<< HEAD
Z = numerl:zeros(2,2),              %
=======
Z = numerl:zeros(2,2),             %
>>>>>>> 199ce279a4b4e7b6d98add6ddfd096c6706475cd
%FCT                               INPUTS
Boolean = numerl:equals(E,Z).      % M1,M2: compared matrices.
```

## Accessors
Matrices content can be extracted with the following functions:
```erlang
Accesors = [mtfli/1, mtfl/1, row/2, col/2, get/3, at/2].
```

They are used as such:
```erlang                             
M = numerl:matrix([[1,2,3]]),
P = numerl:matrix([[1]]),   
N = 1,                             
O = 1,                             
%Output         Fct                 Input
[1,2,3]       = numerl:mtfli(M),   % matrix M
[1.0,2.0,3.0] = numerl:mtfl(M),    % Matrix M
M             = numerl:row(M,N),   % N  in [1, M.n_rows], M
P             = numerl:col(M,N),   % N  in [1, M.n_cols], M
1.0           = numerl:get(M,N,O), % N in
1.0           = numerl:at(M,N).    % The N'th element of the matrix.
```

## Element-wise operations
The following operations can be done element-wise on matrices:

```erlang
Element_wise_ops = [add/2, sub/2 ,mult/2, divide/2].
```

They have the following usage:

```erlang
%FCT                               INPUTS
op(Lval, Rval).                    % Lval: a matrix, Rval: a matrix || a number
```

For ``` divide/2 ``` op used with ``` Rval ``` either null or containing a null value, a badarg is thrown.   

These ops can be combined with the ``` eval/1 ``` function:

```erlang
M = numerl:rnd_matrix(2),
%FCT                               INPUTS
numerl:eval([M, add, 1,            % L: a list of V1,Op1,V2,OP2...
             divide, 2,
              mult, M]).
<<<<<<< HEAD
```

## More functions

```erlang
Ops = [inv/1, nrm2/1, vec_dot/2, dot/2, transpose/1].
```

### Inv

Returns the inverse of input function.
```erlang
M = numerl:rnd_matrix(2),
P = numerl:inv(M).
```
A badarg is thrown for input non-square matrix; and error "nif_inv: could not invert singular matrix." is thrown in case of input singular matrix.   

Implementation is based upon a LAPACKE LU's decomposition/inversion.

### nrm2
 Returns the root square of the sum of the squared elements of input matrix.
 ```erlang
P   = numerl:matrix([[-3]]),
3.0 = numerl:nrm2(P).
```

### vec_dot
 Returns the sum of element-wise multiplication of input matrices.
 ```erlang
Q    = numerl:matrix([[1,2]]),
R    = numerl:matrix([[3,4]]),
11.0 = numerl:vec_dot(Q,R).
```
Input matrices need to contain the same number of elements; but their dimensions do not need to match.

### dot
Returns the product of input matrices. 
 ```erlang
I    = numerl:eye(2),
I2   = numerl:mult(I,2),
M    = numerl:matrix([[1,2], [3,4]]),
R    = numerl:dot(M,I2).
```

### transpose
Returns the transpose of input matrix.
```erlang
I = numerl:eye(2),
It = numerl:transpose(I),
true = numerl:equals(I,It).
```
=======
```

## More functions

```erlang
Ops = [inv/1, nrm2/1, vec_dot/2, dot/2, transpose/1].
```

### Inv

Returns the inverse of input function.
```erlang
M = numerl:rnd_matrix(2),
P = numerl:inv(M).
```
A badarg is thrown for input non-square matrix; and error "nif_inv: could not invert singular matrix." is thrown in case of input singular matrix.   

Implementation is based upon a LAPACKE LU's decomposition/inversion.

### nrm2
 Returns the root square of the sum of the squared elements of input matrix.
 ```erlang
P   = numerl:matrix([[-3]]),
3.0 = numerl:nrm2(P).
```

### vec_dot
 Returns the sum of element-wise multiplication of input matrices.
 ```erlang
Q    = numerl:matrix([[1,2]]),
R    = numerl:matrix([[3,4]]),
11.0 = numerl:vec_dot(Q,R).
```
Input matrices need to contain the same number of elements; but their dimensions do not need to match.

### dot
Returns the product of input matrices. 
 ```erlang
I    = numerl:eye(2),
I2   = numerl:mult(I,2),
M    = numerl:matrix([[1,2], [3,4]]),
R    = numerl:dot(M,I2).
```

### transpose
Returns the transpose of input matrix.
```erlang
I = numerl:eye(2),
It = numerl:transpose(I),
true = numerl:equals(I,It).
```
>>>>>>> 199ce279a4b4e7b6d98add6ddfd096c6706475cd
