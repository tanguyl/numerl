#include <erl_nif.h>
#include <stdio.h>
#include "gsl_cblas.h"



/*
----------------------------------------------------------------------------------------------------|
                        ------------------------------------------                                  |
                        |               INIT FC                  |                                  |
                        ------------------------------------------                                  |
----------------------------------------------------------------------------------------------------|
*/


ERL_NIF_TERM atom_ok;
ERL_NIF_TERM atom_nok;
ERL_NIF_TERM atom_vector;
ERL_NIF_TERM atom_int;
ERL_NIF_TERM atom_double;
ERL_NIF_TERM atom_true;
ERL_NIF_TERM atom_false;

ErlNifResourceType *MULT_YIELDING_ARG = NULL;

int load(ErlNifEnv* env, void** priv_data, ERL_NIF_TERM load_info){
    atom_ok = enif_make_atom(env, "ok\0");
    atom_nok = enif_make_atom(env, "nok");
    atom_vector = enif_make_atom(env, "vector\0");
    atom_int = enif_make_atom(env, "i\0");
    atom_double = enif_make_atom(env, "d\0");
    atom_true = enif_make_atom(env, "true\0");
    atom_false = enif_make_atom(env, "false\0");

    return 0;
}


//First prototype:
//A binary for data,
//A binary for "info".
//Info layout: dimensions, strides.

typedef struct{
    ErlNifBinary content_b;
    ErlNifBinary info_b;

    double* content;
    int* info;
} Array;

int array_get_content_size(Array a){ return a.content_b.size / sizeof(double);}
int array_get_info_size(Array a){return a.info_b.size / sizeof(int);}
int array_get_stride(Array a, int stride){ return a.info[stride + (array_get_info_size(a)/2)]; }
int array_get_dim(Array a, int dimension) {return a.info[dimension];}

int array_eq(Array a, Array b){
    if((array_get_content_size(a) != array_get_content_size(b)) ||
        (array_get_info_size(a) != array_get_info_size(b)))
        return 0;

    for(int i = 0; i<array_get_content_size(a); i++)
        if(a.content[i] != b.content[i]) return 0;

    for(int i = 0; i<array_get_info_size(a); i++)
        if(a.info[i] != b.info[i]) return 0;
        
    return 1;
}

void array_print(Array a){
    printf("Content ");
    for(int i = 0; i<array_get_content_size(a); i++)
        printf(" %f ", a.content[i]);
    
    printf("; Info ");
    for(int i = 0; i<array_get_info_size(a); i++)
        printf(" %d ", a.info[i]);

    printf("\n"); 
}

ERL_NIF_TERM array_to_erl(ErlNifEnv* env, Array a){
    ERL_NIF_TERM content_t = enif_make_binary(env, &a.content_b);
    ERL_NIF_TERM info_t = enif_make_binary(env, &a.info_b);

    return enif_make_tuple2(env, content_t, info_t);
}

Array erl_to_array(ErlNifEnv* env, ERL_NIF_TERM term){
    Array array;
    int arity;
    const ERL_NIF_TERM* content;

    enif_get_tuple(env, term, &arity, &content);
    enif_inspect_binary(env, content[0], &array.content_b);
    enif_inspect_binary(env, content[1], &array.info_b);
    array.content = (double*) array.content_b.data;
    array.info = (int*) array.info_b.data;

    return array;
}


/*
----------------------------------------------------------------------------------------------------|
                        ------------------------------------------                                  |
                        |       Array reading/creation           |                                  |
                        ------------------------------------------                                  |
----------------------------------------------------------------------------------------------------|
*/

/*
@arg 0: Array.
@arg 1: Array.
@return: true if both array are true.
*/
ERL_NIF_TERM nif_array_eq(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
    Array a = erl_to_array(env, argv[0]),
            b = erl_to_array(env, argv[1]);
    
    int equal = array_eq(a,b);

    return equal? atom_true:atom_false;
}

/*
@arg 0: List of Lists of doubles.
@return: Matrix of dimension
*/
ERL_NIF_TERM nif_list_to_matrix(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
    unsigned n, cur_m, m = -1;
    ERL_NIF_TERM list = argv[0], row, elem;
    Array matrix;

    //Reading incoming matrix.
    if(!enif_get_list_length(env, list, &n) && n > 0)
        return atom_nok;

    for(int i = 0; enif_get_list_cell(env, list, &row, &list); i++){
        if(!enif_get_list_length(env, row, &cur_m)) 
            return atom_nok;
        
        if(m == -1){
            //Allocate binary, make matrix accessor.
            m = cur_m;
            enif_alloc_binary(sizeof(double)*m*n, &matrix.content_b);
            enif_alloc_binary(sizeof(int)*4, &matrix.info_b);

            matrix.content = (double*) matrix.content_b.data;
            matrix.info = (int*) matrix.info_b.data;

            matrix.info[0] = m;
            matrix.info[1] = n;
            matrix.info[2] = 1;
            matrix.info[3] = m;
        }

        if(m != cur_m)
            return atom_nok;

        for(int j = 0; enif_get_list_cell(env, row, &elem, &row); j++){
            enif_get_double(env, elem, matrix.content++);
        }
    }

    return array_to_erl(env, matrix);
}


/*
@arg 0: matrix.
@arg 1: int, coord m.
@arg 2: int, coord n.
@return: double, at given coords.
*/
ERL_NIF_TERM nif_get(ErlNifEnv * env, int argc, const ERL_NIF_TERM argv[]){
    ErlNifBinary bin;
    int m,n;
    enif_get_int(env, argv[0], &n);
    enif_get_int(env, argv[1], &m);
    Array matrix = erl_to_array(env, argv[2]);
    

    int index = m*array_get_stride(matrix,0) + n*array_get_stride(matrix, 1);
    return enif_make_double(env, matrix.content[index]);
}


//------------------------------------------------------------------------


ErlNifFunc nif_funcs[] = {
    {"array_eq", 2, nif_array_eq},
    {"list_to_matrix", 1, nif_list_to_matrix},
    {"get", 3, nif_get}
};


ERL_NIF_INIT(numerl, nif_funcs, load, NULL, NULL, NULL)