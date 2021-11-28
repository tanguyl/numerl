-module(numerl).
-on_load(init/0).
-export([ eye/1, zeros/2, '=='/2, '+'/2, '-'/2,'*'/2, matrix/1, get/3, row/2, col/2, tr/1, inv/1, print/1, ddot/3, daxpy/4, dgemv/5, dtrsv/4, dgemm/5]).

-record(matrix, {n_rows, n_cols, info}).

%%  Load nif.
init()->
    erlang:load_nif("./numerl_nif", 0).


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
-spec '=='(M0, M1) -> boolean() when
    M0 :: #matrix{},
    M1 :: #matrix{}.

'=='(_, _) ->
    nif_not_loaded.


%%Addition of matrix.
-spec '+'(M1, M2) -> #matrix{} when
    M1 :: #matrix{},
    M2 :: #matrix{}.

'+'(#matrix{n_rows=N_ROWS, n_cols=N_COLS}, #matrix{n_rows=N_ROWS, n_cols=N_COLS}) ->
    nif_not_loaded.


%%Substraction of matrix.
-spec '-'(M1, M2) -> #matrix{} when
    M1 :: #matrix{},
    M2 :: #matrix{}.

'-'(#matrix{n_rows=N_ROWS, n_cols=N_COLS}, #matrix{n_rows=N_ROWS, n_cols=N_COLS}) ->
    nif_not_loaded.



%% Matrix multiplication.
'*'(A,B) when is_number(A) -> '*_num'(A,B);
'*'(A,B) -> '*_matrix'(A,B).

'*_num'(_,_)->
    nif_not_loaded.

'*_matrix'(#matrix{n_cols=N}, #matrix{n_rows=N})->
    nif_not_loaded.

%% build a null matrix of size NxM
-spec zeros(pos_integer(), pos_integer()) -> #matrix{}.

zeros(_, _) ->
    nif_not_loaded.

%%Returns an Identity matrix NxN.
eye(_)->
    nif_not_loaded.

%%Print in stdout the matrix's content.
print(_)->
    nif_not_loaded.

%Returns the transpose of the given square matrix.
tr(_)->
    nif_not_loaded.

%Returns the inverse of asked square matrix.
inv(_)->
    nif_not_loaded.

%------CBLAS--------

%ddot: dot product of two vectors
%Vectors are matrices with one of their dimension set to 1.
ddot(_,_, _)->
    nif_not_loaded.

%daxpy
daxpy(_,_,_,_)->
    nif_not_loaded.

%dgemv
dgemv(_,_,_,_,_)->
    nif_not_loaded.

%dtrsv
dtrsv(_,_,_,_)->
    nif_not_loaded.

%dgemm
dgemm(_,_,_,_,_)->
    nif_not_loaded.
