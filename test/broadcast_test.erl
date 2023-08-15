-module(broadcast_test).
-include_lib("eunit/include/eunit.hrl").

str_to_list(S)->
    NumsStr = string:lexemes(S, " "),
    lists:map(fun(Str) -> list_to_integer(Str) end, NumsStr).

get_3_arrays(File)->
    Loop = 
        fun
        IterateFile(0, R) -> lists:reverse(R);
        IterateFile(N, R) ->
            case io:get_line(File, "") of
                eof  -> eof;
                Line ->
                    if length(Line) < 4->
                        IterateFile(N, R);
                    true ->
                        [Shape, Content] =string:lexemes(string:trim(Line, trailing, "\n"), "content"),
                        Array = numerl:array(str_to_list(Shape), str_to_list(Content)),
                        IterateFile(N-1, [Array|R])
                    end
            end
        end,
    Loop(3, []).

broadcast_cmp_numpy_test() ->
    PrivDir = case code:priv_dir(my_application) of
                 {error, bad_name} -> "priv";
                 P -> P end,
    {ok, File} = file:open(filename:join(PrivDir, "test_cases.txt"), [read]),

    Loop = 
        fun L(I)->
            case get_3_arrays(File) of
                eof -> 
                    ok;
                [Rhs, Lhs, Result] ->
                    io:format("Testing line : ~b~n", [I]),
                    Obtained = numerl:to_list(numerl:add(Lhs, Rhs)),
                    io:format("Obtained ~w",[Obtained]),
                    Obtained = numerl:to_list(Result),
                    L(I+4);
                _ -> 
                    L(I+1)
            end
        end,            
    Loop(1),
    file:close(File).

