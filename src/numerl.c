#include <erl_nif.h>
#include <stdio.h>
#include <string.h>
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
ERL_NIF_TERM atom_array;

ErlNifResourceType *MULT_YIELDING_ARG = NULL;

int load(ErlNifEnv* env, void** priv_data, ERL_NIF_TERM load_info){
    atom_ok = enif_make_atom(env, "ok\0");
    atom_nok = enif_make_atom(env, "nok");
    atom_vector = enif_make_atom(env, "vector\0");
    atom_int = enif_make_atom(env, "i\0");
    atom_double = enif_make_atom(env, "d\0");
    atom_true = enif_make_atom(env, "true\0");
    atom_false = enif_make_atom(env, "false\0");
    atom_array = enif_make_atom(env, "array\0");

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
int array_get_content_index(Array a, int i, int j){ return array_get_stride(a,0)*i + array_get_stride(a,1)*j;}
double* array_get_content_at(Array a, int i, int j){ return &a.content[array_get_content_index(a,i,j)];}

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

int array_same_dim(Array a, Array b){ 
    if((array_get_content_size(a) != array_get_content_size(b)) ||
        (array_get_info_size(a) != array_get_info_size(b)))
        return 0;
    
     for(int i = 0; i<array_get_content_size(a); i++)
        if(a.content[i] != b.content[i]) return 0;

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


//Allocates memory space of for matrix of requested size.
//Still requires it to get 'sent' to erlang.
Array array_alloc_matrix(int m, int n){
    Array matrix;

    enif_alloc_binary(sizeof(double)*m*n, &matrix.content_b);
    enif_alloc_binary(sizeof(int)*4, &matrix.info_b);

    matrix.content = (double*) matrix.content_b.data;
    matrix.info = (int*) matrix.info_b.data;

    matrix.info[0] = m;
    matrix.info[1] = n;
    matrix.info[2] = 1;
    matrix.info[3] = m;

    return matrix;
}

ERL_NIF_TERM array_to_erl(ErlNifEnv* env, Array a){
    ERL_NIF_TERM content_t = enif_make_binary(env, &a.content_b);
    ERL_NIF_TERM info_t = enif_make_binary(env, &a.info_b);

    return enif_make_tuple3(env, atom_array, content_t, info_t);
}


Array erl_to_array(ErlNifEnv* env, ERL_NIF_TERM term){
    Array array;
    int arity;
    const ERL_NIF_TERM* content;

    enif_get_tuple(env, term, &arity, &content);
    enif_inspect_binary(env, content[1], &array.content_b);
    enif_inspect_binary(env, content[2], &array.info_b);
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
@arg 0: int.
@arg 1: int.
@return: empty array of dimension [0,0]..
*/
ERL_NIF_TERM nif_array_zero(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
    int m,n;
    enif_get_int(env, argv[0], &m);
    enif_get_int(env, argv[1], &n);

    Array a = array_alloc_matrix(m,n);
    memset(a.content, 0, sizeof(double)*m*n);
    return array_to_erl(env, a);
}

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
@arg 0: Array.
@arg 1: Array.
@return Array resulting of element wise + operation.
*/
ERL_NIF_TERM nif_array_plus(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
   
    Array a = erl_to_array(env, argv[0]),
            b = erl_to_array(env, argv[1]);
    
    if (!array_same_dim(a,b)){
        return atom_nok;
    }

    Array result = array_alloc_matrix(array_get_dim(a,0), array_get_dim(a,1));
    
    for(int i = 0; i < array_get_dim(a,0); i++){
        for(int j = 0; j < array_get_dim(a,1); j++){
            (*array_get_content_at(result, i, j)) = (*array_get_content_at(a, i, j)) + (*array_get_content_at(b,i,j));
        }
    }

    return array_to_erl(env, result);
}

/*
@arg 0: Array.
@arg 1: Array.
@return Array resulting of element wise + operation.
*/
ERL_NIF_TERM nif_array_minus(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
   
    Array a = erl_to_array(env, argv[0]),
            b = erl_to_array(env, argv[1]);
    
    if (!array_same_dim(a,b)){
        return atom_nok;
    }

    Array result = array_alloc_matrix(array_get_dim(a,0), array_get_dim(a,1));
    
    for(int i = 0; i < array_get_dim(a,0); i++){
        for(int j = 0; j < array_get_dim(a,1); j++){
            (*array_get_content_at(result, i, j)) = (*array_get_content_at(a, i, j)) - (*array_get_content_at(b,i,j));
        }
    }

    return array_to_erl(env, result);
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
            matrix = array_alloc_matrix(m,n);
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
    enif_get_int(env, argv[0], &m);
    enif_get_int(env, argv[1], &n);
    Array matrix = erl_to_array(env, argv[2]);

    if(m < 0 || m >= array_get_dim(matrix, 0) || n < 0 || n > array_get_dim(matrix, 1))
        return atom_nok;
    

    int index = array_get_content_index(matrix, m, n);
    return enif_make_double(env, matrix.content[index]);
}

/*
@arg 0: int.
@arg 1: Array.
@return: returns an array, containing requested row.
*/
ERL_NIF_TERM nif_row(ErlNifEnv * env, int argc, const ERL_NIF_TERM argv[]){
    int row_req;
    enif_get_int(env, argv[0], &row_req);
    Array matrix = erl_to_array(env, argv[1]);

    int row_size = array_get_dim(matrix, 0);
    Array row = array_alloc_matrix(row_size, 1);


    memcpy(row.content, matrix.content + (row_req * row_size), row_size*sizeof(double));

    return array_to_erl(env, row);
}


/*
@arg 0: int.
@arg 1: Array.
@return: returns an array, containing requested col.
*/
ERL_NIF_TERM nif_col(ErlNifEnv * env, int argc, const ERL_NIF_TERM argv[]){
    int col_req;
    enif_get_int(env, argv[0], &col_req);
    Array matrix = erl_to_array(env, argv[1]);

    int row_size = array_get_dim(matrix, 0);
    int col_size = array_get_dim(matrix, 1);
    Array col = array_alloc_matrix(1, col_size);

    for(int i = 0; i < col_size; i++){
        col.content[i] = matrix.content[i*row_size + col_req];
    }

    return array_to_erl(env, col);
}

//------------------------------------------------------------------------

ErlNifFunc nif_funcs[] = {
    {"zeros", 2, nif_array_zero},
    {"==", 2, nif_array_eq},
    {"+", 2, nif_array_plus},
    {"-", 2, nif_array_minus},
    {"list_to_matrix", 1, nif_list_to_matrix},
    {"get", 3, nif_get},
    {"row", 2, nif_row},
    {"col", 2, nif_col}
};


ERL_NIF_INIT(numerl, nif_funcs, load, NULL, NULL, NULL)