-module(numerl).
-on_load(init/0).
-export([ eye/1, zeros/2, 'equals'/2, add/2, sub/2,dot/2, matrix/1, get/3, row/2, col/2, transpose/1, inv/1, print/1, ddot/3, daxpy/4, dgemv/5, dgemm/5, dgesv/2]).

%Matrices are represented as such:
%-record(matrix, {n_rows, n_cols, bin}).

init()->
  Dir = case code:priv_dir(numerl) of
              {error, bad_name} ->
                  filename:join(
                    filename:dirname(
                      filename:dirname(
                        code:which(?MODULE))), "priv");
              D -> D
          end,
    SoName = filename:join(Dir, atom_to_list(?MODULE)),
    erlang:load_nif(SoName, 0).

%%Creates a matrix.
%List: List of doubles, of length N.
%Return: a matrix of dimension MxN, containing the data.
matrix(_) ->
    nif_not_loaded.

%%Returns a value from a matrix.
get(_,_,_) ->
    nif_not_loaded.

%%Returns requested row. 
row(_,_) ->
    nif_not_loaded.


%%Returns requested col.
col(_,_) ->
    nif_not_loaded.


%%Equality test between matrixes.
'equals'(_, _) ->
    nif_not_loaded.


%%Addition of matrix.

add(_, _) ->
    nif_not_loaded.


%%Substraction of matrix.
sub(_, _) ->
    nif_not_loaded.


%% Matrix multiplication.
dot(A,B) when is_number(A) -> '*_num'(A,B);
dot(A,B) -> '*_matrix'(A,B).

'*_num'(_,_)->
    nif_not_loaded.

'*_matrix'(_, _)->
    nif_not_loaded.


%% build a null matrix of size NxM
zeros(_, _) ->
    nif_not_loaded.

%%Returns an Identity matrix NxN.
eye(_)->
    nif_not_loaded.

%%Print in stdout the matrix's content.
print(_)->
    nif_not_loaded.

%Returns the transpose of the given square matrix.
transpose(_)->
    nif_not_loaded.

%Returns the inverse of asked square matrix.
inv(_)->
    nif_not_loaded.


%------CBLAS--------

% ddot: dot product of two vectors
% Arguments: int n, vector x, vector y.
%   n is the number of coordinates we should take into account (n < min(len(x), len(y)).
%   x and y are matrices, with one of their dimension equal to 1.
% Returns the result of the dot product of the first N coordinates of x, y.
% The returned vector is in the same dimension of y (column or row).
ddot(_,_, _)->
    nif_not_loaded.

% daxpy: alpha*x + y
% Arguments: int n, number alpha, vector x, vector y.
%   n is the number of coordinates we should take into account (n < min(len(x), len(y)).
%   alpha is a number, converted to a double.
%   x and y are matrices, with one of their dimension equal to 1.
% Returns the result of the operation alpha*x + y.
% The returned vector is in the same dimension of y (column or row).
daxpy(_,_,_,_)->
    nif_not_loaded.

% dgemv: alpha*A*x + beta*y
% Arguments: number alpha, Matrix A, vector x, number beta, vector y.
%   alpha and beta re a numbers, converted to doubles.
%   A is a Matrix.
%   x and y are matrices, with one of their dimension equal to 1.
% Returns the result of the operation alpha*A*x + beta*y.
% The returned vector is in the same dimension of y (column or row).
dgemv(_,_,_,_,_)->
    nif_not_loaded.

% dgemm: alpha * A * B + beta * C.
% Arguments: number alpha, Matrix A, Matrix B, number beta, matrix C.
%   alpha, beta: numbers (float or ints) used as doubles.
%   A,B,C: matrices.
% Returns the matrice resulting of the operations alpha * A * B + beta * C.
dgemm(_,_,_,_,_)->
    nif_not_loaded.

%
%dgesv: solve x for A*x = b.
dgesv(_,_)->
    nif_not_loaded.
