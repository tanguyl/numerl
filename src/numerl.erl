-module(numerl).
-export([to_list/1, array/2, array/3, add/2, add/3]).

-record(array, {encoding, shape, data}).


array(Shape, Content) when is_list(Shape)->
    array(float64, Shape, Content);
array(Encoding, Shape) when is_atom(Encoding), is_list(Shape)->
    N_elem = lists:foldr(fun(Dim,Acc)->Dim*Acc end, 1, Shape),
    #array{encoding=Encoding, shape=Shape, data=blas:new(N_elem*elem_size(Encoding))}.

array(Encoding, Shape, Content)->
    Data = blas:new(Encoding, Content),
    #array{encoding=Encoding, shape=Shape, data=Data}.


to_list(#array{encoding=Encoding, data=Data})->
    Binary = blas:to_bin(Data),
    case Encoding of 
        int64       -> [ V || <<V:64/native-integer>> <= Binary ];
        float32     -> [ V || <<V:32/native-float>>   <= Binary ];
        float64     -> [ V || <<V:64/native-float>>   <= Binary ];
        complex64   -> [ V || <<V:32/native-float>>   <= Binary ];
        complex128  -> [ V || <<V:64/native-float>>   <= Binary ]
    end.


elem_size(Encoding)->
    case Encoding of
        int64       -> 8;
        float32     -> 4;
        float64     -> 8;
        complex64   -> 8;
        complex128  -> 16
    end.

%Copy from lhs to rhs the N elements.
copy(#array{encoding=Encoding, data=Ld}, #array{encoding=Encoding, data=Rd}, N)->
    case Encoding of 
        float32    -> blas:run({scopy, N, Ld, 1, Rd, 1});
        float64    -> blas:run({dcopy, N, Ld, 1, Rd, 1});
        complex64  -> blas:run({ccopy, N, Ld, 1, Rd, 1});
        complex128 -> blas:run({zcopy, N, Ld, 1, Rd, 1})
    end.

add(Lhs=#array{encoding=Encoding,shape=LShape}, Rhs=#array{encoding=Encoding,shape=RShape})->
    % Unify the shapes: unified left,right, out shapes
    Shapes = {ULS,URS,UOS} = lists:unzip3(lists:foldr(
        fun
            ({V,V}, Acc) -> [{V,V,V}|Acc];
            ({V,1}, Acc) -> [{V,1,V}|Acc];
            ({1,V}, Acc) -> [{1,V,V}|Acc]
        end,
        [],
        lists:reverse(lists:zip(lists:reverse(LShape), lists:reverse(RShape), {pad, {1,1}}))
    )),
    %io:format("Shapes are ~w~n", [Shapes]),
    Out   = array(Encoding, UOS),
    broadcast(Lhs#array{shape=ULS}, Rhs#array{shape=URS}, Out, fun add/3).

           
add({array, Encoding, I, Ld}, Rhs, D={array, Encoding, N, Dd})->
    {LI,LN} = {lists:last(I), lists:last(N)},
    K = if LI == 1 -> 0; true -> 1 end,
    L = if LN == 1 -> 0; true -> 1 end,
    copy(Rhs,D, LN),
    case Encoding of
        float32    -> blas:run({saxpy, LN, 1.0, Ld, K, Dd, L});
        float64    -> blas:run({daxpy, LN, 1.0, Ld, K, Dd, L});
        complex64  -> blas:run({caxpy, LN, [1.0, 0.0], Ld, K, Dd, L});
        complex128 -> blas:run({zaxpy, LN, [1.0, 0.0], Ld, K, Dd, L})
    end,
    D.



shift(Arr={array, Encoding, Shape, Data}, Steps)->
    % In shape, for each step, shift accordingly the array.
    {DimSize,_} = lists:mapfoldr(fun(V,Acc)-> {Acc, V*Acc} end, elem_size(Encoding), Shape),
    Shift = lists:foldr(
        fun
            ({1, _, _},          Acc)                     -> Acc;
            ({IShape, IStep, _}, Acc) when IShape < IStep -> Acc;
            ({_,IStep, ISize},   Acc)                     -> Acc + ISize*IStep 
        end,
        0,
        apply(lists, zip3, lists:map(fun lists:droplast/1, [Shape, Steps, DimSize]))
    ),
    io:format("Applying shift ~w to ~w ~n", [Shift, Arr]),
    Arr#array{data=blas:shift(Shift, Data)}.

broadcast(Lhs=#array{encoding=Encoding},
          Rhs=#array{encoding=Encoding},
          Out=#array{encoding=Encoding, shape=OShape},
          Fct
)->    
    CurIt  = [0 || _ <- OShape],

    Iteration = fun LocalIt(ILhs,IRhs,IOut, ItCounter)->
        io:format("~nNew iteration: it is ~w.~n", [ItCounter]),
        Fct(ILhs, IRhs, IOut),
        if 
            length(OShape) == 1 ->
                Out;
            true ->
                {NextIt, Stop} = lists:mapfoldr(
                    fun({Step,Max},Acc) -> 
                        Val = Step+Acc,
                        if 
                            Val >= Max -> {0,1};
                            true       -> {Val,0}
                        end
                    end,
                    lists:last(OShape) + 1,
                    lists:zip(ItCounter, OShape)
                ),
                io:format("Finished? ~w~n", [Stop]),
                if 
                    Stop == 1 -> Out;
                    true           -> LocalIt(shift(Lhs,NextIt), shift(Rhs,NextIt), shift(Out,NextIt), NextIt)
                end
            end
        end,
        Iteration(Lhs, Rhs, Out, CurIt).


%matmul({array, Encoding, [N,M], Ld}, {array, Encoding, [M,K], Rd})->
%    Ed = blas:new(N*K*elem_size(Encoding)),
%%    case Encoding of 
%        float32    -> blas:run({sgemm, blasRowMajor, n,n, N,M,K, 1.0,  Ld,N, Rd,M, 0.0,  Ed,N});
%        float64    -> blas:run({dgemm, blasRowMajor, n,n, N,M,K, 1.0,  Ld,N, Rd,M, 0.0,  Ed,N});
%%        complex64  -> blas:run({cgemm, blasRowMajor, n,n, N,M,K, [1,0],Ld,N, Rd,M, [0,0],Ed,N});
%        complex128 -> blas:run({zgemm, blasRowMajor, n,n, N,M,K, [1,0],Ld,N, Rd,M, [0.0],Ed,N})
%    end,
%    {array, Encoding, [N,K], Ed}.

%dot({array, Encoding, [N], Ld}, {array, Encoding, [N], Rd})->
%    case Encoding of 
%        float32    -> blas:run({sdot, N, Ld, 1, Rd, 1});
%        float64    -> blas:run({ddot, N, Ld, 1, Rd, 1});
%        complex64  -> blas:run({cdot, N, Ld, 1, Rd, 1});
%        complex128 -> blas:run({zdot, N, Ld, 1, Rd, 1})
%    end.