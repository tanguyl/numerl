-module(numerl_tests).
-compile(export_all).
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
    V00 = 1.0,
    V00 = numerl:get(1,1,CM0),

    V11 = 4.0,
    V11 = numerl:get(2,2,CM0),

    V01 = 2.0,
    V01 = numerl:get(1,2,CM0),

    V10 = 3.0,
    V10 = numerl:get(2,1, CM0).


equal_test() ->
    M0 = [[1.0, 2.0], [3.0, 4.0]],
    M1 = [[1.0, 2.0]],
    M2 = [[1.0], [2.0]],
    CM0 = numerl:matrix(M0),
    CM1 = numerl:matrix(M1),
    CM2 = numerl:matrix(M2),
    true = numerl:equals(CM1, CM1),
    true = numerl:equals(CM0, CM0),
    false = numerl:equals(CM1, CM2),
    false = numerl:equals(CM0, CM2),
    false = numerl:equals(CM0, CM1).

row_test() ->
    M0 = [[1.0, 2.0], [3.0, 4.0]],
    R0 = [[3.0, 4.0]],
    CM0 = numerl:matrix(M0),
    CR0 = numerl:matrix(R0),
    true = numerl:equals(CR0, numerl:row(2, CM0)).


col_test() ->
    M0 = [[1.0, 2.0], [3.0, 4.0]],
    C0 = [[2.0], [4.0]],
    CM0 = numerl:matrix(M0),
    CC0 = numerl:matrix(C0),
    true = numerl:equals(CC0, numerl:col(2, CM0)).

plus_test()->
     M0 = [[1.0, 2.0], [3.0, 4.0]],
     M1 = [[2.0, 4.0], [6.0, 8.0]],
     CM0 = numerl:matrix(M0),
     CM1 = numerl:matrix(M1),
     CM0p = numerl:add(CM0, CM0),
     true = numerl:equals(CM1, CM0p).

minus_test()->
     M0 = [[1, 2], [3, 4]],
     M1 = [[0, 0], [0, 0]],
     CM0 = numerl:matrix(M0),
     CM1 = numerl:matrix(M1),
     CM0p = numerl:sub(CM0, CM0),
     true = numerl:equals(CM1, CM0p).

zero_test() ->
    CM0 = numerl:zeros(1,5),
    0.0 = numerl:get(1, 5, CM0).

eye_test() ->
    CM0 = numerl:eye(5),
    M0 = numerl:get(5,5,CM0),
    M0 = 1.0,
    M0 = numerl:get(1,1,CM0),
    M0 = numerl:get(2,2,CM0),
    M0 = numerl:get(4,4,CM0).

mult_num_test()->
    M0 = numerl:matrix([[1.0, 2.0]]),
    true = numerl:equals(numerl:matrix([[2,4]]), numerl:dot(2, M0)),
    true = numerl:equals(numerl:matrix([[0,0]]), numerl:dot(0, M0)),
    true = numerl:equals(numerl:matrix([[-1, -2]]), numerl:dot(-1, M0)).

mult_matrix_test() ->
    CM0 = numerl:eye(2),
    CM1 = numerl:matrix([[1.0, 2.0], [3.0, 4.0]]),
    CM3 = numerl:matrix([[1.0], [2.0]]),
    CM5 = numerl:matrix([[5.0], [11.0]]),
    CM4 = numerl:matrix([[1.0, 2.0]]),
    CM6 = numerl:matrix([[7.0, 10.0]]),
    true = numerl:equals(CM1, numerl:dot(CM1, CM0)),
    true = numerl:equals(CM5, numerl:dot(CM1, CM3)),
    true = numerl:equals(CM6 ,numerl:dot(CM4, CM1)).

mult_matrix_num_test() ->
    M0 = [[1.0, 2.0, 3.0]],
    CR0 = numerl:matrix([[2, 4, 6]]),
    CR = numerl:dot(2.0, numerl:matrix(M0)),
    true = numerl:equals(CR0, CR).



tr_test() ->
    CM0 = numerl:eye(2),
    CM0 = numerl:transpose(CM0),
    CM1 = numerl:matrix([[1.0, 2.0],[3.0, 4.0]]),
    CM2 = numerl:matrix([[1.0, 3.0], [2.0, 4.0]]),
    CM3 = numerl:matrix([[1.0, 2.0]]),
    true = numerl:equals(CM3, numerl:transpose(numerl:matrix([[1.0], [2.0]]))),
    CM1 = numerl:transpose(CM2).

inv_test() ->
    M = numerl:matrix([[2.0, -1.0, 0.0], [-1.0, 2.0, -1.0], [0.0, -1.0, 2.0]]),
    M_inv = numerl:inv(M),
    io:fwrite(numerl:print(M_inv),[]),
    true = numerl:equals(numerl:dot(M, M_inv), numerl:eye(3)).

ddot_test() ->
    Incs = numerl:matrix([[1, 2, 3, 4]]),
    Ones = numerl:matrix([[1], [1], [1], [1]]),
    10.0 = numerl:ddot(4, Incs, Ones),
    30.0 = numerl:ddot(4, Incs, Incs),
    4.0 = numerl:ddot(4, Ones, Ones),
    1.0 = numerl:ddot(1, Incs, Ones).

daxpy_test()->
    Ones = numerl:matrix([[1, 1, 1, 1]]),
    Incs = numerl:matrix([[1, 2, 3, 4]]),
    true = numerl:equals(numerl:matrix([[3, 4, 5, 6]]), numerl:daxpy(4, 2, Ones, Incs)).

dgemv_test()->
    V10 = numerl:matrix([[1,2]]),
    V01 = numerl:matrix([[0], [1]]),
    M = numerl:matrix([[1,2],[3,4]]),
    true = numerl:equals(numerl:matrix([[10],[26]]), numerl:dgemv(2,M,V10, 4, V01)).

dgemm_test()->
    A = numerl:matrix([[1,2]]),
    B = numerl:matrix([[3,4], [5,6]]),
    C = numerl:matrix([[10, 12]]),
    true = numerl:equals(numerl:matrix([[31, 38]]), numerl:dgemm(2,A,B,0.5,C)).

dgesv_test()->
    A = numerl:matrix([[1,2],[0,3]]),
    B = numerl:matrix([[2,4],[0,6]]),
    X = numerl:matrix([[2,0], [0,2]]),
    true = numerl:equals(X, numerl:dgesv(A, B)).
