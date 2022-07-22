-module(numerl).
-on_load(init/0).
-export([ eval/1, eye/1, zeros/2, equals/2, add/2, sub/2,mult/2, divide/2, sqrt/1, matrix/1, rnd_matrix/1, get/3, at/2, mtfli/1, mtfl/1, row/2, col/2, transpose/1, inv/1, nrm2/1, vec_dot/2, dot/2]).

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

%Creates a random matrix.
rnd_matrix(N)->
    L = [[rand:uniform(20) || _ <- lists:seq(1,N) ] || _ <- lists:seq(1,N)],
    matrix(L).

%Combine multiple functions.
eval([L,O,R|T])->
    F = fun numerl:O/2,
    eval([F(L,R) |T]);
eval([Res])->
    Res.

%%Creates a matrix.
%List: List of doubles, of length N.
%Return: a matrix of dimension MxN, containing the data.
matrix(_) ->
    nif_not_loaded.

%%Returns the Nth value contained within Matrix.
at(_Matrix,_Nth)->
    nif_not_loaded.

%%Returns the matrix as a flattened list of ints.
mtfli(_mtrix)->
    nif_not_loaded.

%%Returns the matrix as a flattened list of doubles.
mtfl(_mtrix)->
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
equals(_, _) ->
    nif_not_loaded.


%%Addition of matrix.
add(_, _) ->
    nif_not_loaded.


%%Substraction of matrix.
sub(_, _) ->
    nif_not_loaded.


%% Matrix multiplication.
mult(A,B) when is_number(B) -> '*_num'(A,B);
mult(A,B) -> '*_matrix'(A,B).

'*_num'(_,_)->
    nif_not_loaded.

'*_matrix'(_, _)->
    nif_not_loaded.

%Matrix division by a number
divide(_,_)->
    nif_not_loaded.

%Returns the square root of each element of a matrix.
sqrt(_)->
    nif_not_loaded.


%% build a null matrix of size NxM
zeros(_, _) ->
    nif_not_loaded.

%%Returns an Identity matrix NxN.
eye(_)->
    nif_not_loaded.

%Returns the transpose of the given square matrix.
transpose(_)->
    nif_not_loaded.

%Returns the inverse of asked square matrix.
inv(_)->
    nif_not_loaded.


%------CBLAS--------

%nrm2
%Calculates the squared root of the sum of the squared contents.
nrm2(_)->
    nif_not_loaded.

% : dot product of two vectors
% Arguments: vector x, vector y.
%   x and y are matrices
% Returns the dot product of all the coordinates of X,Y.
vec_dot(_, _)->
    nif_not_loaded.

% dgemm: A dot B 
% Arguments: Matrix A, Matrix B.
%   alpha, beta: numbers (float or ints) used as doubles.
%   A,B,C: matrices.
% Returns the matrice resulting of the operations alpha * A * B + beta * C.
dot(_,_)->
    nif_not_loaded.