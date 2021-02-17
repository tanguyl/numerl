-module(numerl).

-on_load(init).


-init_nif()->
    nif_not_loaded.

-init()->
    erlang:load_nif("./numerl.c", 0),
    init_nif().