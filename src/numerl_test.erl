-module(numerl_test).
-include_lib("eunit/include/eunit.hrl").
    

list_to_matrix_test() ->
    M0 = [[1.0, 0.0], [0.0, 1.0]],
    _ = numerl:list_to_matrix(M0).

get_test() ->
    %Testing access on square matrix
    M0 = [[1.0, 2.0], [3.0, 4.0]],
    CM0 = numerl:list_to_matrix(M0),
    V00 = mat:get(1,1,M0),
    V00 = numerl:get(0,0,CM0),
    V11 = mat:get(2,2,M0),
    V11 = numerl:get(1,1,CM0),
    V01 = mat:get(1,2,M0),
    V01 = numerl:get(0,1,CM0),
    V10 = mat:get(2,1,M0),
    V10 = numerl:get(1,0, CM0).

equal_test() ->
    M0 = [[1.0, 2.0], [3.0, 4.0]],
    M1 = [[1.0, 2.0]],
    M2 = [[1.0], [2.0]],
    CM0 = numerl:list_to_matrix(M0),
    CM1 = numerl:list_to_matrix(M1),
    CM2 = numerl:list_to_matrix(M2),
    true = numerl:array_eq(CM1, CM1),
    true = numerl:array_eq(CM0, CM0),
    false = numerl:array_eq(CM1, CM2),
    false = numerl:array_eq(CM0, CM2),
    false = numerl:array_eq(CM0, CM1).



