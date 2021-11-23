-module(numerl_test).
-include_lib("eunit/include/eunit.hrl").

matrix_test() ->
    M0 = [[1.0, 0.0], [0.0, 1.0]],
    _ = numerl:matrix(M0).

print_test() ->
    M0 = numerl:matrix([[1.0/3.0, 0.0], [0.0, 1.0/3.0]]),
    numerl:print(M0).


get_test() ->
    %Testing access on square matrix
    M0 = [[1.0, 2.0], [3.0, 4.0]],
    CM0 = numerl:matrix(M0),
    V00 = mat:get(1,1,M0),
    V00 = numerl:get(1,1,CM0),

    V11 = mat:get(2,2,M0),
    V11 = numerl:get(2,2,CM0),

    V01 = mat:get(1,2,M0),
    V01 = numerl:get(1,2,CM0),

    V10 = mat:get(2,1,M0),
    V10 = numerl:get(2,1, CM0),

    M1 = [[1.0], [2.0]],
    CM1 = numerl:matrix(M1),
    W00 = mat:get(1, 1, M1),
    W00 = numerl:get(1, 1, CM1),
    W01 = mat:get(2, 1, M1),
    W01 = numerl:get(2, 1, CM1),

    M2 = [[1.0, 2.0]],
    CM2 = numerl:matrix(M2),
    X00 = mat:get(1,2,M2),
    X00 = numerl:get(1,2,CM2).


equal_test() ->
    M0 = [[1.0, 2.0], [3.0, 4.0]],
    M1 = [[1.0, 2.0]],
    M2 = [[1.0], [2.0]],
    CM0 = numerl:matrix(M0),
    CM1 = numerl:matrix(M1),
    CM2 = numerl:matrix(M2),
    true = numerl:'=='(CM1, CM1),
    true = numerl:'=='(CM0, CM0),
    false = numerl:'=='(CM1, CM2),
    false = numerl:'=='(CM0, CM2),
    false = numerl:'=='(CM0, CM1).

row_test() ->
    M0 = [[1.0, 2.0], [3.0, 4.0]],
    R0 = mat:row(2, M0),
    CM0 = numerl:matrix(M0),
    CR0 = numerl:matrix(R0),
    true = numerl:'=='(CR0, numerl:row(2, CM0)).


col_test() ->
    M0 = [[1.0, 2.0], [3.0, 4.0]],
    C0 = mat:col(2, M0),
    CM0 = numerl:matrix(M0),
    CC0 = numerl:matrix(C0),
    true = numerl:'=='(CC0, numerl:col(2, CM0)).

plus_test()->
     M0 = [[1.0, 2.0], [3.0, 4.0]],
     M1 = [[2.0, 4.0], [6.0, 8.0]],
     CM0 = numerl:matrix(M0),
     CM1 = numerl:matrix(M1),
     CM0p = numerl:'+'(CM0, CM0),
     true = numerl:'=='(CM1, CM0p).

minus_test()->
     M0 = [[1, 2], [3, 4]],
     M1 = [[0, 0], [0, 0]],
     CM0 = numerl:matrix(M0),
     CM1 = numerl:matrix(M1),
     CM0p = numerl:'-'(CM0, CM0),
     true = numerl:'=='(CM1, CM0p).

zero_test() ->
    M0 = mat:zeros(1,5),
    CM0 = numerl:zeros(1,5),
    G0 = float(mat:get(1, 5, M0)),
    G0 = numerl:get(1, 5, CM0).

eye_test() ->
    CM0 = numerl:eye(5),
    M0 = numerl:get(5,5,CM0),
    M0 = 1.0,
    M0 = numerl:get(1,1,CM0),
    M0 = numerl:get(2,2,CM0),
    M0 = numerl:get(4,4,CM0).

mult_num_test()->
    M0 = numerl:matrix([[1.0, 2.0]]),
    true = numerl:'=='(numerl:matrix([[2,4]]), numerl:'*'(2, M0)),
    true = numerl:'=='(numerl:matrix([[0,0]]), numerl:'*'(0, M0)),
    true = numerl:'=='(numerl:matrix([[-1, -2]]), numerl:'*'(-1, M0)).

mult_matrix_test() ->
    CM0 = numerl:eye(2),
    CM1 = numerl:matrix([[1.0, 2.0], [3.0, 4.0]]),
    CM3 = numerl:matrix([[1.0], [2.0]]),
    CM5 = numerl:matrix([[5.0], [11.0]]),
    CM4 = numerl:matrix([[1.0, 2.0]]),
    CM6 = numerl:matrix([[7.0, 10.0]]),
    true = numerl:'=='(CM1, numerl:'*'(CM1, CM0)),
    true = numerl:'=='(CM5, numerl:'*'(CM1, CM3)),
    true = numerl:'=='(CM6 ,numerl:'*'(CM4, CM1)).

mult_matrix_num_test() ->
    M0 = [[1.0, 2.0, 3.0]],
    CR0 = numerl:matrix(mat:'*'(2.0, M0)),
    CR = numerl:'*'(2.0, numerl:matrix(M0)),
    true = numerl:'=='(CR0, CR).



tr_test() ->
    CM0 = numerl:eye(2),
    CM0 = numerl:tr(CM0),
    CM1 = numerl:matrix([[1.0, 2.0],[3.0, 4.0]]),
    CM2 = numerl:matrix([[1.0, 3.0], [2.0, 4.0]]),
    CM3 = numerl:matrix([[1.0, 2.0]]),
    true = numerl:'=='(CM3, numerl:tr(numerl:matrix([[1.0], [2.0]]))),
    CM1 = numerl:tr(CM2).

inv_test() ->
    M = numerl:matrix([[2.0, -1.0, 0.0], [-1.0, 2.0, -1.0], [0.0, -1.0, 2.0]]),
    M_inv = numerl:inv(M),
    io:fwrite(numerl:print(M_inv),[]),
    true = numerl:'=='(numerl:'*'(M, M_inv), numerl:eye(3)).


