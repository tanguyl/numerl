#include <erl_nif.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
----------------------------------------------------------------------------------------------------|
                        ------------------------------------------                                  |
                        |               INIT FC                  |                                  |
                        ------------------------------------------                                  |
----------------------------------------------------------------------------------------------------|
*/


ERL_NIF_TERM atom_nok;
ERL_NIF_TERM atom_true;
ERL_NIF_TERM atom_false;
ERL_NIF_TERM atom_matrix;

ErlNifResourceType *MULT_YIELDING_ARG = NULL;

int load(ErlNifEnv* env, void** priv_data, ERL_NIF_TERM load_info){
    atom_nok = enif_make_atom(env, "nok\0");
    atom_true = enif_make_atom(env, "true\0");
    atom_false = enif_make_atom(env, "false\0");
    atom_matrix = enif_make_atom(env, "matrix\0");
    return 0;
}

int upgrade(ErlNifEnv* caller_env, void** priv_data, void** old_priv_data, ERL_NIF_TERM load_info){

}

//Gives easier access to an ErlangBinary containing a matrix.
typedef struct{
    int n_rows;
    int n_cols;
    
    //Content of the matrix, in row major format.
    double* content;

    //Erlang binary containing the matrix.
    ErlNifBinary bin;
} Matrix;


//Access asked coordinates of matrix
double* matrix_at(int col, int row, Matrix m){
    return &m.content[row*m.n_cols + col];
}

//Allocates memory space of for matrix of dimensions n_rows, n_cols.
//The matrix_content can be modified, until a call to array_to_erl.
//Matrix content is stored in row major format.
Matrix alloc_matrix(int n_rows, int n_cols){
    ErlNifBinary bin;

    enif_alloc_binary(sizeof(double)*n_rows*n_cols, &bin);

    Matrix matrix;
    matrix.n_cols = n_cols;
    matrix.n_rows = n_rows;
    matrix.bin = bin;
    matrix.content = (double*) bin.data;

    return matrix;
}

//Free an allocated matrix that was not sent back to Erlang.
void free_matrix(Matrix m){
    enif_release_binary(&m.bin);
}

//Constructs a matrix erlang term.
//No modifications can be made afterwards to the matrix.
ERL_NIF_TERM matrix_to_erl(ErlNifEnv* env, Matrix m){
    ERL_NIF_TERM term = enif_make_binary(env, &m.bin);
    return enif_make_tuple4(env, atom_matrix, enif_make_int(env,m.n_rows), enif_make_int(env,m.n_cols), term);
}

//Reads an erlang term as a matrix.
//As such, no modifications can be made to the red matrix.
//Returns true if it was possible to read a matrix, false otherwise
int enif_get_matrix(ErlNifEnv* env, ERL_NIF_TERM term, Matrix *dest){
    
    int arity;
    const ERL_NIF_TERM* content;

    enif_get_tuple(env, term, &arity, &content);
    if(arity != 4)
        return 0;
    
    if(content[0] != atom_matrix
        || !enif_get_int(env, content[1], &dest->n_rows)
        || !enif_get_int(env, content[2], &dest->n_cols)
        || !enif_inspect_binary(env, content[3], &dest->bin))
    {
        return 0;
    }

    dest->content = (double*) (dest->bin.data);    
    return 1;
}

//Equal all doubles
//Compares wether all doubles are approximatively the same.
int equal_ad(double* a, double* b, int size){
    for(int i = 0; i<size; i++){
        if(fabs(a[i] - b[i])> 1e-6)
            return 0;
    }
    return 1;
}

//Used to translate at once a number of ERL_NIF_TERM.
//Data types are infered via content of format string:
//  n: number (int or double) translated to double.
//  m: matrix
//  i: int
int enif_get(ErlNifEnv* env, const ERL_NIF_TERM* erl_terms, const char* format, ...){
    va_list valist;
    va_start(valist, format);
    int valid = 1;

    while(valid && *format != '\0'){
        switch(*format++){
            case 'n':
                //Read a number as a double.
                ;
                double *d = va_arg(valist, double*);
                int i;
                if(!enif_get_double(env, *erl_terms, d))
                    if(enif_get_int(env, *erl_terms, &i))
                        *d = (double) i;
                    else valid = 0;
                break;

            case 'm':
                //Reads a matrix.
                valid = enif_get_matrix(env, *erl_terms, va_arg(valist, Matrix*));
                break;
            
            case 'i':
                //Reads an int.
                valid = enif_get_int(env, *erl_terms, va_arg(valist, int*));
                break;
            
            default:
                //Unknown type... give an error.
                valid = 0;
                break;
        }
        erl_terms ++;
    }

    va_end(valist);
    return valid;
}


//----------------------------------------------------------------------------------------------------|
//                        ------------------------------------------                                  |
//                        |                   NIFS                 |                                  |
//                        ------------------------------------------                                  |
//----------------------------------------------------------------------------------------------------|


//@arg 0: List of Lists of numbers.
//@return: Matrix of dimension
ERL_NIF_TERM nif_matrix(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
    unsigned n_rows, line_length, dest = 0, n_cols = -1;
    ERL_NIF_TERM list = argv[0], row, elem;
    Matrix m;

    //Reading incoming matrix.
    if(!enif_get_list_length(env, list, &n_rows) && n_rows > 0)
        return enif_make_badarg(env);

    for(int i = 0; enif_get_list_cell(env, list, &row, &list); i++){
        if(!enif_get_list_length(env, row, &line_length)) 
            return enif_make_badarg(env);
        
        if(n_cols == -1){
            //Allocate binary, make matrix accessor.
            n_cols = line_length;
            m = alloc_matrix(n_rows,n_cols);
        }

        if(n_cols != line_length)
            return enif_make_badarg(env);

        for(int j = 0; enif_get_list_cell(env, row, &elem, &row); j++){
            if(!enif_get_double(env, elem, &m.content[dest])){
                int i;
                if(enif_get_int(env, elem, &i)){
                    m.content[dest] = (double) i;
                }
                else{
                    return enif_make_badarg(env);
                }
            }
            dest++;
        }
    }

    return matrix_to_erl(env, m);
}

#define PRECISION 10
//Arg0: a valid matrix.
//Returns: an atom representation of the input matrix.
ERL_NIF_TERM matrix_to_atom(ErlNifEnv *env, Matrix m){
    char *content = enif_alloc(sizeof(char)*((2*m.n_cols-1)*m.n_rows*PRECISION + m.n_rows*2 + 3));
    content[0] = '[';
    content[1] = '\0';
    char converted[PRECISION];

    for(int i=0; i<m.n_rows; i++){
        strcat(content, "[");
        for(int j = 0; j<m.n_cols; j++){
            snprintf(converted, PRECISION-1, "%.5lf", m.content[i*m.n_cols+j]);
            strcat(content, converted);
            if(j != m.n_cols-1)
                strcat(content, " ");
        }
        strcat(content, "]");
    }
    strcat(content, "]");
    
    ERL_NIF_TERM result = enif_make_atom(env, content);
    enif_free(content);
    return result;
}

//@arg 0: Matrix.
//@return Nothing.
ERL_NIF_TERM nif_matrix_print(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
   
    Matrix m;
    if(!enif_get(env, argv, "m", &m))
        return enif_make_badarg(env);
    return matrix_to_atom(env, m);
}


//@arg 0: matrix.
//@arg 1: int, coord m: row
//@arg 2: int, coord n: col
//@return: double, at cord Matrix(m,n).
ERL_NIF_TERM nif_get(ErlNifEnv * env, int argc, const ERL_NIF_TERM argv[]){
    ErlNifBinary bin;
    int m,n;
    Matrix matrix;

    if(!enif_get(env, argv, "iim", &m, &n, &matrix))
        return enif_make_badarg(env);
    m--, n--;

    if(m < 0 || m >= matrix.n_rows || n < 0 || n >= matrix.n_cols)
        return enif_make_badarg(env);

    int index = m*matrix.n_cols+n;
    return enif_make_double(env, matrix.content[index]);
}


//@arg 0: Array.
//@arg 1: Array.
//@return: true if arrays share content, false if they have different content || size..
ERL_NIF_TERM nif_eq(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
    Matrix a,b;
    if(!enif_get(env, argv, "mm", &a, &b))
        return enif_make_badarg(env);

    //Compare number of columns and rows
    if((a.n_cols != b.n_cols || a.n_rows != b.n_rows))
        return atom_false;

    //Compare content of arrays
    if(!equal_ad(a.content, b.content, a.n_cols*a.n_rows))
        return atom_false;
    
    return atom_true;
}


//@arg 0: int.
//@arg 1: Array.
//@return: returns an array, containing requested row.
ERL_NIF_TERM nif_row(ErlNifEnv * env, int argc, const ERL_NIF_TERM argv[]){
    int row_req;
    Matrix matrix;
    if(!enif_get(env, argv, "im", &row_req, &matrix))
        return enif_make_badarg(env);

    row_req--;
    if(row_req<0 || row_req >= matrix.n_rows)
        return enif_make_badarg(env);

    Matrix row = alloc_matrix(1, matrix.n_cols);
    memcpy(row.content, matrix.content + (row_req * matrix.n_cols), matrix.n_cols*sizeof(double));

    return matrix_to_erl(env, row);
}


//@arg 0: int.
//@arg 1: Array.
//@return: returns an array, containing requested col.
ERL_NIF_TERM nif_col(ErlNifEnv * env, int argc, const ERL_NIF_TERM argv[]){
    int col_req;
    Matrix matrix;
    if(!enif_get(env, argv, "im", &col_req, &matrix))
        return enif_make_badarg(env);

    col_req--;    
    if(col_req<0 || col_req >= matrix.n_cols)
        return enif_make_badarg(env);


    Matrix col = alloc_matrix(matrix.n_rows, 1);

    for(int i = 0; i < matrix.n_rows; i++){
        col.content[i] = matrix.content[i * matrix.n_rows + col_req];
    }

    return matrix_to_erl(env, col);
}


//@arg 0: Matrix.
//@arg 1: Matrix.
//@return Matrix resulting of element wise + operation.
ERL_NIF_TERM nif_plus(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
   
    Matrix a,b;
    if(!enif_get(env, argv, "mm", &a, &b))
        return enif_make_badarg(env);
    
    if ((a.n_cols != b.n_cols || a.n_rows != b.n_rows)){
        return enif_make_badarg(env);
    }

    Matrix result = alloc_matrix(a.n_rows, a.n_cols);
    
    for(int i = 0; i < a.n_cols; i++){
        for(int j = 0; j < a.n_rows; j++){
            *matrix_at(i,j, result) = *matrix_at(i,j,a) + *matrix_at(i,j,b);
        }
    }

    return matrix_to_erl(env, result);
}

//@arg 0: Matrix.
//@arg 1: Matrix.
//@return Matrix resulting of element wise - operation.
ERL_NIF_TERM nif_minus(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
   
    Matrix a,b;
    if(!enif_get(env, argv, "mm", &a, &b))
        return enif_make_badarg(env);

    if (a.n_cols != b.n_cols || a.n_rows != b.n_rows){
        return enif_make_badarg(env);
    }

    Matrix result = alloc_matrix(a.n_rows, a.n_cols);
    
    for(int i = 0; i < a.n_cols; i++){
        for(int j = 0; j < a.n_rows; j++){
            *matrix_at(i,j, result) = *matrix_at(i,j,a) - *matrix_at(i,j,b);
        }
    }

    return matrix_to_erl(env, result);
}

//@arg 0: int.
//@arg 1: int.
//@return: empty Matrix of requested dimension..
ERL_NIF_TERM nif_zero(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
    int m,n;
    if(!enif_get(env, argv, "ii", &m, &n))
        return enif_make_badarg(env);

    Matrix a = alloc_matrix(m,n);
    memset(a.content, 0, sizeof(double)*m*n);
    return matrix_to_erl(env, a);
}

//@arg 0: int.
//@arg 1: int.
//@return: empty matrix of dimension [arg 0, arg 1]..
ERL_NIF_TERM nif_eye(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
    int m;
    if(!enif_get_int(env, argv[0], &m))
        return enif_make_badarg(env);
    
    if(m <= 0)
        return enif_make_badarg(env);

    Matrix a = alloc_matrix(m,m);
    memset(a.content, 0, sizeof(double)*m*m);
    for(int i = 0; i<m; i++){
        a.content[i*m+i] = 1.0;
    }
    return matrix_to_erl(env, a);
}

//Transpose of a matrix.
Matrix tr(Matrix a){
    Matrix result = alloc_matrix(a.n_cols ,a.n_rows);

    for(int j = 0; j < a.n_rows; j++){
        for(int i = 0; i < a.n_cols; i++){
            result.content[i*result.n_cols+j] = a.content[j*a.n_cols+i];
        }
    }
    return result;
}

//arg 0: double or int
//arg 1: Matrix
//@return the result of multiplying each matrix element by arg 0.
ERL_NIF_TERM nif_mult_num(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){

    double a;
    Matrix b;
    if(!enif_get(env, argv, "nm", &a, &b))
        return enif_make_badarg(env);

    Matrix c = alloc_matrix(b.n_rows, b.n_cols);

    for(int i = 0; i<b.n_cols*b.n_rows; i++){
        c.content[i] = a * b.content[i];
    }

    return matrix_to_erl(env, c);
}

//@arg 0: Matrix.
//@arg 1: Matrix.
//@return Matrix resulting of multiplication.
ERL_NIF_TERM _nif_mult_matrix(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
   
    Matrix a,b;
    if(!enif_get(env, argv, "mm", &a, &b))
        return enif_make_badarg(env);
        
    int n_rows = a.n_rows;
    int n_cols = b.n_cols;

    if(a.n_cols != b.n_rows)
        return atom_nok;

    Matrix result = alloc_matrix(n_rows, n_cols);
    memset(result.content, 0.0, n_rows*n_cols * sizeof(double));

    //Will this create a mem leak?
    Matrix b_tr = tr(b); 
    
    for(int i = 0; i < n_rows; i++){
        for(int j = 0; j < n_cols; j++){
           for(int k = 0; k<a.n_cols; k++){
               result.content[j+i*result.n_cols] += a.content[k+i*a.n_cols] * b_tr.content[k+j*b_tr.n_cols];
           }
        }
    }

    free_matrix(b_tr);

    return matrix_to_erl(env, result);
}


//@arg0: Matrix.
ERL_NIF_TERM nif_tr(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){

    Matrix a;
    if(!enif_get_matrix(env, argv[0], &a))
        return enif_make_badarg(env);

    return matrix_to_erl(env, tr(a));
}


//arg0: Matrix.
ERL_NIF_TERM nif_inv(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
    
    Matrix a;
    if(!enif_get_matrix(env, argv[0], &a))
        return enif_make_badarg(env);

    if(a.n_cols != a.n_rows){
        return atom_nok;
    }
    int n_cols = 2*a.n_cols;

    double* gj = (double*) enif_alloc(n_cols*a.n_rows*sizeof(double));
    for(int i=0; i<a.n_rows; i++){
        memcpy(gj+i*n_cols, a.content+i*a.n_cols, sizeof(double)*a.n_cols);
        memset(gj+i*n_cols + a.n_cols, 0, sizeof(double)*a.n_cols);
        gj[i*n_cols + a.n_cols + i] = 1.0;
    }

    //Elimination de Gauss Jordan:
    //https://fr.wikipedia.org/wiki/%C3%89limination_de_Gauss-Jordan
   
    //Row of last found pivot
    int r = -1;
    //j for all indexes of column
    for(int j=0; j<a.n_cols; j++){

        //Find the row of the maximum in column j
        int pivot_row = -1;
        for(int cur_row=r; cur_row<a.n_rows; cur_row++){
            if(pivot_row<0 || fabs(gj[cur_row*n_cols + j]) > fabs(gj[pivot_row*n_cols+j])){
                pivot_row = cur_row;
            }
        }
        double pivot_value = gj[pivot_row*n_cols+j];

        if(pivot_value != 0){
            r++;
            for(int cur_col=0; cur_col<n_cols; cur_col++){
                gj[cur_col+pivot_row*n_cols] /= pivot_value;
            }
            gj[pivot_row*n_cols + j] = 1.0; //make up for rounding errors

            //Do we need to swap?
            if(pivot_row != r){
                for(int i = 0; i<n_cols; i++){
                    double cpy = gj[pivot_row*n_cols+i];
                    gj[pivot_row*n_cols+i]= gj[r*n_cols+i];
                    gj[r*n_cols+i] = cpy;
                }
            }

            //We can simplify all the rows
            for(int i=0; i<a.n_rows; i++){
                if(i!=r){
                    double factor = gj[i*n_cols+j];
                    for(int col=0; col<n_cols; col++){
                        gj[col+i*n_cols] -= gj[col+r*n_cols]*factor;
                    }
                    gj[i*n_cols+j] = 0.0;    //make up for rounding errors
                }
            }
        }

    }
    
    Matrix inv = alloc_matrix(a.n_rows, a.n_cols);
    for(int l=0; l<inv.n_rows; l++){
        int line_start = l*n_cols + a.n_cols;
        memcpy(inv.content + inv.n_cols*l, gj + line_start, sizeof(double)*inv.n_cols);
    }

    enif_free(gj);
    return matrix_to_erl(env, inv);
}


//------------------------------------------------------------------------


ErlNifFunc nif_funcs[] = {
    
    {"matrix", 1, nif_matrix},
    {"print", 1, nif_matrix_print},
    {"get", 3, nif_get},
    {"==", 2, nif_eq},
    {"row", 2, nif_row},
    {"col", 2, nif_col},
    {"+", 2, nif_plus},
    {"-", 2, nif_minus},
    {"zeros", 2, nif_zero},
    {"eye", 1, nif_eye},
    {"*_matrix", 2, _nif_mult_matrix},
    {"*_num", 2, nif_mult_num},
    {"tr", 1, nif_tr},
    {"inv", 1, nif_inv}
};

ERL_NIF_INIT(numerl, nif_funcs, load, NULL, NULL, NULL)