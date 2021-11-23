-module(benchmark).
-compile(export_all).
-import(math, [log10/1]).
-import(timer, [sleep/1, tc/3]).


%Used to measure the average run time required for a function to run.
%Measure in micro seconds.
bench(_,_,_,0,R)->
    R;

bench(F, Args,N,I,R)->
    {Time, _} = timer:tc(F, Args),
    bench(F, Args, N,I-1, R+(Time/float(N))).

%Run the function N time to get it's average performance.
bench(F ,Args ,N) ->
    bench(F, Args, N, N, 0.0).

%Run the function F a number of times.
bench(F, Args)->
    bench(F, Args, 5000).


%Creates a random matrix of size NxN
rnd_matrix(N)->
    [ [rand:uniform(100) || _ <- lists:seq(1, N)] || _ <- lists:seq(1,N)].

%Creates a random int
rnd() ->
    rand:uniform(100).

%Prints the given results.
show_results(Name, T_e, T_n)->
    io:format("~nTesting ~w\nErlang native:  ~f\nNif: ~f\nFactor:~f~n", [Name,T_e, T_n, T_e/T_n]).

%Runs test for Mat_fct, Num_fct functions having as argument I, an int.
%The test is printed alongside Name.
cmp_fct_1i(Size, Name, Mat_fct, Num_fct)->
    T_e = bench(Mat_fct, [Size]),
    T_n = bench(Num_fct, [Size]),
    show_results(Name, T_e, T_n).

%Runs test for Mat_fct, Num_fct functions having as argument a square matrix of size Size x Size.
%The test is printed alongside Name.
cmp_fct_1(Size, Name, Mat_fct, Num_fct)->
    M_e = rnd_matrix(Size),
    M_n = numerl:matrix(M_e),
    T_e = bench(Mat_fct, [M_e]),
    T_n = bench(Num_fct, [M_n]),
    show_results(Name, T_e, T_n).

%Runs test for Mat_fct, Num_fct having as argument a number + a square matrix of Size x Size.
cmp_fct_2im(Size, Name, Mat_fct, Num_fct)->
    Num = rnd() / 2.0,
    M_e = rnd_matrix(Size),
    M_n = numerl:matrix(M_e),
    T_e = bench(Mat_fct, [Num,  M_e]),
    T_n = bench(Num_fct, [Num, M_n]),
    show_results(Name, T_e, T_n).


%Runs test for Mat_fct, Num_fct functions having as argument two square matrix of size Size x Size.
%The test is printed alongside Name.
cmp_fct_2(Size, Name, Mat_fct, Num_fct)->
    M_e = rnd_matrix(Size),
    M_e2 = rnd_matrix(Size),
    M_n = numerl:matrix(M_e),
    M_n2 = numerl:matrix(M_e2),
    T_e = bench(Mat_fct, [M_e, M_e2]),
    T_n = bench(Num_fct, [M_n, M_n2]),
    show_results(Name, T_e, T_n).

%Test how long does it takes to create a matrix
bench_matrix(Size)->
    M_e = rnd_matrix(Size),
    T_n = bench(fun numerl:matrix/1, [M_e]),
    show_results(matrix, 0.0, T_n).

%Compare performances for all benchmarks.
run(Size)->
    io:format("Input size set to ~B~n", [Size]),
    cmp_fct_1i(Size, eye, fun mat:eye/1, fun numerl:eye/1),
    cmp_fct_1(Size, inversion, fun mat:inv/1, fun numerl:inv/1),
    cmp_fct_1(Size, transpose, fun mat:tr/1, fun numerl:tr/1),
    cmp_fct_2(Size, mult, fun mat:'*'/2, fun numerl:'*'/2),
    cmp_fct_2im(Size, mult, fun mat:'*'/2, fun numerl:'*'/2).

