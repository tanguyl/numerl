-module(numerl).
-on_load(init/0).
-export([ array_eq/2, list_to_matrix/1, get/3]).

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

%%Equality test between matrixes.
array_eq(_, _) ->
    nif_not_loaded.