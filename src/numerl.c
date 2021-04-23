#include <erl_nif.h>
#include <stdio.h>
#include <string.h>
#include "gsl_cblas.h"
#include <math.h>



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


typedef struct{
    //Size 0: size of dimensions, size 1: size of content.
    int* sizes;
    //The size of each dimensions.
    int* dimensions;
    //Single array content of the dimensions, in row major format.
    double* content;

    //Erlang binary's data is pointed at by sizes, etc
    ErlNifBinary bin;
} Array;


//Allocates memory space of for matrix of requested size.
//Caution: does not initialise memory!
/*
Array.sizes
Array.dimension[0]: n_rows
Array.dimension[1]: n_cols
*/
Array alloc_matrix(int n_rows, int n_cols){
    Array matrix;

    int content_size = n_rows*n_cols;
    enif_alloc_binary(sizeof(double)*content_size + sizeof(int)*4, &matrix.bin);

    matrix.sizes = (int*)matrix.bin.data;
    matrix.sizes[0] = 2;
    matrix.sizes[1] = content_size;

    matrix.dimensions = matrix.sizes + 2;
    matrix.dimensions[0] = n_rows;
    matrix.dimensions[1] = n_cols;

    matrix.content = (double*) (matrix.bin.data + 4 * sizeof(int)); 

    return matrix;
}

/*
'Uploads' an array to it's erlang representation.
*/
ERL_NIF_TERM array_to_erl(ErlNifEnv* env, Array a){
    ERL_NIF_TERM term = enif_make_binary(env, &a.bin);

    return enif_make_tuple2(env, atom_array, term);
}


/*
'Downloads' an Erlang represented array to an Array.
*/
Array erl_to_array(ErlNifEnv* env, ERL_NIF_TERM term){
    Array array;
    int arity;
    const ERL_NIF_TERM* content;

    enif_get_tuple(env, term, &arity, &content);
    enif_inspect_binary(env, content[1], &array.bin);

    //Reading
    array.sizes = (int*) array.bin.data;
    array.dimensions = array.sizes + 2;
    array.content = (double*) array.bin.data + sizeof(int) * (array.sizes[0] + 2);
    
    return array;
}

int equal_ai(int* i, int* j, int cur){
    while(--cur >= 0 && i[cur] == j[cur]);
    return cur < 0;
}

int equal_d(double i, double j){ return fabs(i-j) <= 0.00000001;}
int equal_ad(double* i, double* j, int cur){
    while(--cur >= 0 && equal_d(i[cur], j[cur]));
    return cur < 0;
}


//----------------------------------------------------------------------------------------------------|
//                        ------------------------------------------                                  |
//                        |       Array reading/creation           |                                  |
//                        ------------------------------------------                                  |
//----------------------------------------------------------------------------------------------------|


//@arg 0: List of Lists of numbers.
//@return: Matrix of dimension
ERL_NIF_TERM nif_matrix(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
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

//@arg 0: Matrix.
//@return Nothing.
ERL_NIF_TERM nif_matrix_print(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
   
    Array a = erl_to_array(env, argv[0]);
    printf("Matrix [%d, %d], (%d) [", a.dimensions[0], a.dimensions[1], a.sizes[1]);

    for(int i = 0; i < a.sizes[1]; i++){
        printf("%lf ,", a.content[i]);
    } 

    printf("]\n");
    return atom_ok;
}

/*
//@arg 0: int.
//@arg 1: int.
//@return: empty array of dimension [0,0]..
ERL_NIF_TERM nif_array_zero(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
    int m,n;
    enif_get_int(env, argv[0], &m);
    enif_get_int(env, argv[1], &n);

    Array a = alloc_matrix(m,n);
    memset(a.content, 0, sizeof(double)*m*n);
    return array_to_erl(env, a);
}


//@arg 0: int.
//@arg 1: int.
//@return: empty array of dimension [0,0]..
ERL_NIF_TERM nif_array_eye(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
    int m;
    enif_get_int(env, argv[0], &m);

    Array a = array_alloc_matrix(m,m);
    memset(a.content, 0, sizeof(double)*m*m);
    for(int i = 0; i<m; i++){
        a.content[i*m+i] = 1.0;
    }
    return array_to_erl(env, a);
}



//@arg 0: Array.
//@arg 1: Array.
//@return: true if arrays share content.
ERL_NIF_TERM nif_array_eq(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
    Array a = erl_to_array(env, argv[0]),
            b = erl_to_array(env, argv[1]);
    
    if(a.sizes[0] != b.sizes[0] || )

    return equal? atom_true:atom_false;
}


//@arg 0: Array.
//@arg 1: Array.
//@return Array resulting of element wise + operation.
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


//@arg 0: Array.
//@arg 1: Array.
//@return Array resulting of element wise + operation.
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


//@arg 0: matrix.
//@arg 1: int, coord m.
//@arg 2: int, coord n.
//@return: double, at given coords.
ERL_NIF_TERM nif_get(ErlNifEnv * env, int argc, const ERL_NIF_TERM argv[]){
    ErlNifBinary bin;
    int m,n;
    enif_get_int(env, argv[0], &n);
    enif_get_int(env, argv[1], &m);
    Array matrix = erl_to_array(env, argv[2]);

    if(m < 0 || m >= array_get_dim(matrix, 1) || n < 0 || n > array_get_dim(matrix, 0))
        return atom_nok;
    

    int index = array_get_content_index(matrix, m, n);
    return enif_make_double(env, matrix.content[index]);
}


//@arg 0: int.
//@arg 1: Array.
//@return: returns an array, containing requested row.
ERL_NIF_TERM nif_row(ErlNifEnv * env, int argc, const ERL_NIF_TERM argv[]){
    int row_req;
    enif_get_int(env, argv[0], &row_req);
    Array matrix = erl_to_array(env, argv[1]);

    int row_size = array_get_dim(matrix, 0);
    Array row = array_alloc_matrix(row_size, 1);


    memcpy(row.content, matrix.content + (row_req * row_size), row_size*sizeof(double));

    return array_to_erl(env, row);
}



//@arg 0: int.
//@arg 1: Array.
//@return: returns an array, containing requested col.
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



//@arg 0: Array.
//@arg 1: Array.
//@return Array resulting of multiplication.
ERL_NIF_TERM nif_array_mult(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
   
    Array a = erl_to_array(env, argv[0]),
            b = erl_to_array(env, argv[1]);
        
    int m = array_get_dim(a,1);
    int n = array_get_dim(b,0);

    if(array_get_dim(a,0) != array_get_dim(b,1))
        return atom_nok;

    Array result = array_alloc_matrix(m, n);

    printf("About to work\n");
    
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            printf("(%d, %d); ", i, j);
            (*array_get_content_at(result, j, i)) = 0.0;
            for(int k = 0; k < array_get_dim(a,1); k++){
                printf("(%lf += %lf * %lf);  ", (*array_get_content_at(result, i, j)), (*array_get_content_at(a, k, i)), (*array_get_content_at(b,j,k)));
                (*array_get_content_at(result, j, i)) += (*array_get_content_at(a, k, i)) * (*array_get_content_at(b,j,k));   
            }
        }
    }

    return array_to_erl(env, result);
}


*/
//------------------------------------------------------------------------

ErlNifFunc nif_funcs[] = {
    
    {"list_to_matrix", 1, nif_list_to_matrix},
    {"print", 1, nif_array_print},
    //{"get", 3, nif_get},
    /*
    {"eye", 1, nif_array_eye},
    {"zeros", 2, nif_array_zero},
    {"==", 2, nif_array_eq},
    {"+", 2, nif_array_plus},
    {"-", 2, nif_array_minus},
    {"*", 2, nif_array_mult},
    {"row", 2, nif_row},
    {"col", 2, nif_col},
    */
};


ERL_NIF_INIT(numerl, nif_funcs, load, NULL, NULL, NULL)