-module(numerl_tests).
-include_lib("eunit/include/eunit.hrl").


add_test()->
    X = numerl:array([3], [1,2,3]),
    Y = numerl:array([1], [1]),
    
    Z = numerl:add(Y,X),
    %W = numerl:add(X,Y),
    A = numerl:add(X,X),
    [2.0,3.0,4.0]   = numerl:to_list(Z),
    %[2.0,3.0,4.0]   = numerl:to_list(W),
    [2.0, 4.0, 6.0] = numerl:to_list(A).


    
    