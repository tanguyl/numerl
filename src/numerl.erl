-module(numerl).
-on_load(init/0).
-export([ eye/1, zeros/2, '=='/2, '+'/2, '-'/2, list_to_matrix/1, get/3, row/2, col/2]).

-record(array, {content, info}).

%%  Load nif.
init()->
    erlang:load_nif("./numerl_nif", 0).


%%Creates a matrix.
%List: List of doubles, of length N.
%Return: a matrix of dimension MxN, containing the data.
list_to_matrix(_) ->
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

%Returns an empty matrix.
zeros(N, M) ->
    nif_not_loaded.

%%Returns an Identity matrix NxN.
eye(N)->
    nif_not_loaded.