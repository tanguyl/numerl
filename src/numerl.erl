-module(numerl).
-on_load(init/0).
-export([ eye/1, zeros/2, '=='/2, '+'/2, '-'/2,'*'/2, matrix/1, get/3, row/2, col/2, tr/1, inv/1, print/1]).

-record(array, {content, info}).

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
    M0 :: #array{},
    M1 :: #array{}.

'=='(_, _) ->
    nif_not_loaded.


%%Addition of matrix.
-spec '+'(M1, M2) -> _ when
    M1 :: #array{},
    M2 :: #array{}.

'+'(_, _) ->
    nif_not_loaded.

%%Substraction of matrix.
-spec '-'(M1, M2) -> _ when
    M1 :: #array{},
    M2 :: #array{}.

'-'(_, _) ->
    nif_not_loaded.

%% build a null matrix of size NxM
-spec zeros(N, M) -> Zeros when
    N :: pos_integer(),
    M :: pos_integer(),
    Zeros :: #array{}.

%% Matrix multiplication.
-spec '*'(M1, M2) -> _ when
    M1 :: #array{},
    M2 :: #array{}.

'*'(_, _)->
    nif_not_loaded.

%Returns an empty matrix.
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