-module(numerl_tests).
-include_lib("eunit/include/eunit.hrl").

matrix_test() ->
    M0 = [[1.0, 0.0], [0.0, 1.0]],
    _ = numerl:matrix(M0).

get_test() ->
    %Testing access on square matrix
    M0 = [[1.0, 2.0], [3.0, 4.0]],
    CM0 = numerl:matrix(M0),
    V00 = 1.0,
    V00 = numerl:get(CM0, 1,1),

    V11 = 4.0,
    V11 = numerl:get(CM0,2,2),

    V01 = 2.0,
    V01 = numerl:get(CM0,1,2),

    V10 = 3.0,
    V10 = numerl:get(CM0,2,1).

at_test()->
    0.0 = numerl:at(numerl:matrix([[1,0]]), 2),
    2.0 = numerl:at(numerl:matrix([[2,3,4]]), 1).

mtfli_test()->
    [1,2,3] = numerl:mtfli(numerl:matrix([[1.1, 2.9, 3]])).

mtfl_test()->
    [1.1,2.1,3.4] = numerl:mtfl(numerl:matrix([[1.1, 2.1, 3.4]])).

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
    false = numerl:equals(CM0, CM1),
    false = numerl:equals(CM0,1).

row_test() ->
    M0 = [[1.0, 2.0], [3.0, 4.0]],
    R0 = [[3.0, 4.0]],
    CM0 = numerl:matrix(M0),
    CR0 = numerl:matrix(R0),
    true = numerl:equals(CR0, numerl:row(CM0,2)).


col_test() ->
    M0 = [[1.0, 2.0], [3.0, 4.0]],
    C0 = [[2.0], [4.0]],
    CM0 = numerl:matrix(M0),
    CC0 = numerl:matrix(C0),
    true = numerl:equals(CC0, numerl:col(CM0,2)).


zero_test() ->
    CM0 = numerl:zeros(1,5),
    0.0 = numerl:get(CM0, 1, 5).

eye_test() ->
    CM0 = numerl:eye(5),
    M0 = numerl:get(CM0, 5,5),
    M0 = 1.0,
    M0 = numerl:get(CM0, 1,1),
    M0 = numerl:get(CM0, 2,2),
    M0 = numerl:get(CM0, 4,4).


add_test()->
     CM0 = numerl:matrix([[1, 2], [3, 4]]),
     CM1 = numerl:matrix([[2, 4], [6, 8]]),
     CM3 = numerl:matrix([[2, 3], [4, 5]]),
     true = numerl:equals(CM1, numerl:add(CM0, CM0)),
     true = numerl:equals(CM3, numerl:add(CM0,1)).


sub_test()->
     CM0 = numerl:matrix([[1, 2], [3, 4]]),
     CM1 = numerl:matrix([[2, 4], [6, 8]]),
     CM3 = numerl:matrix([[-1, -2], [-3, -4]]),
     CM4 = numerl:matrix([[0, 1], [2, 3]]),
     true = numerl:equals(CM3, numerl:sub(CM0, CM1)),
     true = numerl:equals(CM4, numerl:sub(CM0,1)).

mult_test()->
    CM0 = numerl:matrix([[1, 2], [3, 4]]),
    CM1 = numerl:matrix([[2, 4], [6, 8]]),
    CM3 = numerl:matrix([[2, 8], [18, 32]]),
    CM4 = numerl:matrix([[2, 4], [6, 8]]),
    true = numerl:equals(CM3, numerl:mult(CM0, CM1)),
    true = numerl:equals(CM4, numerl:mult(CM0,2)).

divide_test()->
    CM0 = numerl:matrix([[1, 2], [3, 4]]),
    CM1 = numerl:matrix([[2, 4], [6, 8]]),
    CM3 = numerl:matrix([[0.5, 0.5], [0.5, 0.5]]),
    CM4 = numerl:matrix([[0.5, 1], [1.5, 2]]),
    true = numerl:equals(CM3, numerl:divide(CM0, CM1)),
    true = numerl:equals(CM4, numerl:divide(CM0,2)),
    badarg_error = 
    try
        numerl:divide(numerl:matrix([[1]]), 0),
        no_error
    catch 
        error:badarg -> badarg_error
    end,
    badarg_error = 
    try
        numerl:divide(numerl:matrix([[1],[2]]), numerl:matrix([[1],[0]])),
        no_error
    catch 
        error:badarg -> badarg_error
    end.


tr_test() ->
    CM0 = numerl:eye(2),
    CM0 = numerl:transpose(CM0),
    CM1 = numerl:matrix([[1.0, 2.0],[3.0, 4.0]]),
    CM2 = numerl:matrix([[1.0, 3.0], [2.0, 4.0]]),
    CM3 = numerl:matrix([[1.0, 2.0]]),
    true = numerl:equals(CM3, numerl:transpose(numerl:matrix([[1.0], [2.0]]))),
    CM1 = numerl:transpose(CM2).

inv_test() ->
    %Might cause an error if randomly generated matrix is singular.
    F = fun()->
        N = rand:uniform(20),
        M = numerl:rnd_matrix(N),
        M_inv = numerl:inv(M),
        numerl:equals(numerl:dot(M, M_inv), numerl:eye(N))
    end,
    List = [F() || _ <- lists:seq(1,50)],
    lists:all(fun(_)-> true end, List).



vec_dot_test() ->
    Incs = numerl:matrix([[1, 2, 3, 4]]),
    Ones = numerl:matrix([[1], [1], [1], [1]]),
    10.0 = numerl:vec_dot(Incs, Ones),
    30.0 = numerl:vec_dot(Incs, Incs),
    4.0 = numerl:vec_dot(Ones, Ones).

dot_test()->
    A = numerl:matrix([[1,2]]),
    B = numerl:matrix([[3,4], [5,6]]),
    true = numerl:equals(numerl:matrix([[13, 16]]), numerl:dot(A,B)).

memleak_test()->
    %For input matrices of size 10: run all function once, check memory, run a couple more times, check if memory increase.
    N = 10,
    AllFcts = fun()->
        M = numerl:rnd_matrix(N),
        _ = numerl:at(M,1),
        _ = numerl:mtfli(M),
        _ = numerl:mtfl(M),
        _ = numerl:equals(M,M),
        _ = numerl:add(M,M),
        _ = numerl:add(M,2),
        _ = numerl:sub(M,M),
        _ = numerl:sub(M,2),
        _ = numerl:mult(M,M),
        _ = numerl:mult(M,2),
        _ = numerl:divide(M,M),
        _ = numerl:divide(M,2),
        _ = numerl:transpose(M),
        _ = numerl:inv(M),
        _ = numerl:nrm2(M),
        _ = numerl:vec_dot(M,M),
        _ = numerl:dot(M,M)
        end,
    
    AllFcts(),
    erlang:garbage_collect(),
    {memory, M_first_run} = erlang:process_info(self(), memory),

    AllFcts(), AllFcts(), AllFcts(), AllFcts(),

    erlang:garbage_collect(),
    {memory, M_second_run} = erlang:process_info(self(), memory),

    M_first_run >= M_second_run.


    
    