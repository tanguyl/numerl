#include <erl_nif.h>
#include <stdio.h>
#include <string.h>
#include <math.h> 
#include <cblas.h>


/*
--------------------------------------------------------------------|
            ------------------------------                          |
            |          LAPACKE           |                          |
            ------------------------------                          |
--------------------------------------------------------------------|
*/

void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);


/*
--------------------------------------------------------------------|
            ------------------------------                          |
            |           INIT FC          |                          |
            ------------------------------                          |
--------------------------------------------------------------------|
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
    return 0;
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
Matrix matrix_alloc(int n_rows, int n_cols){
    ErlNifBinary bin;

    enif_alloc_binary(sizeof(double)*n_rows*n_cols, &bin);

    Matrix matrix;
    matrix.n_cols = n_cols;
    matrix.n_rows = n_rows;
    matrix.bin = bin;
    matrix.content = (double*) bin.data;

    return matrix;
}

//Creates a duplicate of a matrix.
//This duplicate can be modified untill uploaded.
Matrix matrix_dup(Matrix m){
    Matrix d = matrix_alloc(m.n_rows, m.n_cols);
    memcpy(d.content, m.content, d.bin.size);
    return d;
}

//Free an allocated matrix that was not sent back to Erlang.
void matrix_free(Matrix m){
    enif_release_binary(&m.bin);
}

//Constructs a matrix erlang term.
//No modifications can be made afterwards to the matrix.
ERL_NIF_TERM matrix_to_erl(ErlNifEnv* env, Matrix m){
    ERL_NIF_TERM term = enif_make_binary(env, &m.bin);
    return enif_make_tuple4(env, atom_matrix, enif_make_int(env,m.n_rows), enif_make_int(env,m.n_cols), term);
}

int enif_is_matrix(ErlNifEnv* env, ERL_NIF_TERM term){
    int arity;
    const ERL_NIF_TERM* content;

    if(!enif_is_tuple(env, term))
        return 0;

    enif_get_tuple(env, term, &arity, &content);
    if(arity != 4)
        return 0;
    
    if(content[0] != atom_matrix
        || !enif_is_number(env, content[1])
        || !enif_is_number(env, content[2])
        || !enif_is_binary(env, content[3]))
            
            return 0;
    return 1;
}

//Reads an erlang term as a matrix.
//As such, no modifications can be made to the red matrix.
//Returns true if it was possible to read a matrix, false otherwise
int enif_get_matrix(ErlNifEnv* env, ERL_NIF_TERM term, Matrix *dest){
    
    int arity;
    const ERL_NIF_TERM* content;

    if(!enif_is_tuple(env, term))
        return 0;

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
            m = matrix_alloc(n_rows,n_cols);
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

//Used for debug purpose.
void debug_write(char info[]){
    FILE* fp = fopen("debug.txt", "a");
    fprintf(fp, info);
    fprintf(fp, "\n");
    fclose(fp);
}

void debug_write_matrix(Matrix m){
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
    debug_write(content);
    enif_free(content);

}


//@arg 0: matrix.
//@arg 1: int, coord m: row
//@arg 2: int, coord n: col
//@return: double, at cord Matrix(m,n).
ERL_NIF_TERM nif_get(ErlNifEnv * env, int argc, const ERL_NIF_TERM argv[]){
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

ERL_NIF_TERM nif_at(ErlNifEnv * env, int argc, const ERL_NIF_TERM argv[]){
    int n;
    Matrix matrix;

    if(!enif_get(env, argv, "mi", &matrix, &n))
        return enif_make_badarg(env);
    n--;

    if( n < 0 || n >= matrix.n_cols * matrix.n_rows)
        return enif_make_badarg(env);

    return enif_make_double(env, matrix.content[n]);
}

//Matrix to flattened list of ints
ERL_NIF_TERM nif_mtfli(ErlNifEnv * env, int argc, const ERL_NIF_TERM argv[]){
    Matrix M;
    if(!enif_get(env, argv, "m", &M)){
        return enif_make_badarg(env);
    }

    int n_elems = M.n_cols * M.n_rows;
    ERL_NIF_TERM *arr = enif_alloc(sizeof(ERL_NIF_TERM)*n_elems);
    for(int i = 0; i<n_elems; i++){
        arr[i] = enif_make_int(env, (int)M.content[i]);
    }
    
    ERL_NIF_TERM result = enif_make_list_from_array(env, arr, n_elems);
    enif_free(arr);
    return result;
}

//Matrix to flattened list of ints
ERL_NIF_TERM nif_mtfl(ErlNifEnv * env, int argc, const ERL_NIF_TERM argv[]){
    Matrix M;
    if(!enif_get(env, argv, "m", &M)){
        return enif_make_badarg(env);
    }

    int n_elems = M.n_cols * M.n_rows;
    ERL_NIF_TERM *arr = enif_alloc(sizeof(ERL_NIF_TERM)*n_elems);
    for(int i = 0; i<n_elems; i++){
        arr[i] = enif_make_double(env, M.content[i]);
    }
    
    ERL_NIF_TERM result = enif_make_list_from_array(env, arr, n_elems);
    enif_free(arr);
    return result;
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

//@arg 0: Array.
//@arg 1: Array.
//@return: true if arrays share content, false if they have different content || size..
ERL_NIF_TERM nif_equals(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
    Matrix a,b;

    if(!enif_get(env, argv, "mm", &a, &b))
        return atom_false;

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

    Matrix row = matrix_alloc(1, matrix.n_cols);
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


    Matrix col = matrix_alloc(matrix.n_rows, 1);

    for(int i = 0; i < matrix.n_rows; i++){
        col.content[i] = matrix.content[i * matrix.n_rows + col_req];
    }

    return matrix_to_erl(env, col);
}


//@arg 0: int.
//@arg 1: int.
//@return: empty Matrix of requested dimension..
ERL_NIF_TERM nif_zeros(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
    int m,n;
    if(!enif_get(env, argv, "ii", &m, &n))
        return enif_make_badarg(env);

    Matrix a = matrix_alloc(m,n);
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

    Matrix a = matrix_alloc(m,m);
    memset(a.content, 0, sizeof(double)*m*m);
    for(int i = 0; i<m; i++){
        a.content[i*m+i] = 1.0;
    }
    return matrix_to_erl(env, a);
}

//Either element-wise multiplication, or multiplication of a matrix by a number.
ERL_NIF_TERM nif_mult(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){

    Matrix a,b;
    double c;
    
    if(enif_get(env, argv, "mn", &a, &c)){
        Matrix d = matrix_alloc(a.n_rows, a.n_cols);

        for(int i = 0; i<a.n_cols*a.n_rows; i++){
            d.content[i] = a.content[i] * c;
        }
        return matrix_to_erl(env, d);
    }
    else 
    if(enif_get(env, argv, "mm", &a, &b)
            && a.n_rows*a.n_cols == b.n_rows*b.n_cols){
        
        Matrix d = matrix_alloc(a.n_rows, a.n_cols);

        for(int i = 0; i < a.n_rows*a.n_cols; i++){
            d.content[i] = a.content[i] * b.content[i];
        }
        return matrix_to_erl(env, d);
    }

    return enif_make_badarg(env);
}


//Either element-wise addition, or addition of a matrix by a number.
ERL_NIF_TERM nif_add(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){

    Matrix a,b;
    double c;
    
    if(enif_get(env, argv, "mn", &a, &c)){
        Matrix d = matrix_alloc(a.n_rows, a.n_cols);

        for(int i = 0; i<a.n_cols*a.n_rows; i++){
            d.content[i] = a.content[i] + c;
        }
        return matrix_to_erl(env, d);
    }
    else 
    if(enif_get(env, argv, "mm", &a, &b)
            && a.n_rows*a.n_cols == b.n_rows*b.n_cols){
        
        Matrix d = matrix_alloc(a.n_rows, a.n_cols);

        for(int i = 0; i < a.n_rows*a.n_cols; i++){
            d.content[i] = a.content[i] + b.content[i];
        }
        return matrix_to_erl(env, d);
    }

    return enif_make_badarg(env);
}


//Either element-wise substraction, or substraction of a matrix by a number.
ERL_NIF_TERM nif_sub(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){

    Matrix a,b;
    double c;
    
    if(enif_get(env, argv, "mn", &a, &c)){
        Matrix d = matrix_alloc(a.n_rows, a.n_cols);

        for(int i = 0; i<a.n_cols*a.n_rows; i++){
            d.content[i] = a.content[i] - c;
        }
        return matrix_to_erl(env, d);
    }
    else 
    if(enif_get(env, argv, "mm", &a, &b)
            && a.n_rows*a.n_cols == b.n_rows*b.n_cols){
        
        Matrix d = matrix_alloc(a.n_rows, a.n_cols);

        for(int i = 0; i < a.n_rows*a.n_cols; i++){
            d.content[i] = a.content[i] - b.content[i];
        }
        return matrix_to_erl(env, d);
    }

    return enif_make_badarg(env);
}


//Either element-wise division, or division of a matrix by a number.
ERL_NIF_TERM nif_divide(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){

    Matrix a,b;
    double c;
    
    if(enif_get(env, argv, "mn", &a, &c)){
        
        if(c!=0.0){
            Matrix d = matrix_alloc(a.n_rows, a.n_cols);

            for(int i = 0; i<a.n_cols*a.n_rows; i++){
                d.content[i] = a.content[i] / c;
            }
            return matrix_to_erl(env, d);
        }
    }
    else if(enif_get(env, argv, "mm", &a, &b)
            && a.n_rows*a.n_cols == b.n_rows*b.n_cols){
        
        Matrix d = matrix_alloc(a.n_rows, a.n_cols);

        for(int i = 0; i < a.n_rows*a.n_cols; i++){
            if(b.content[i] == 0.0){
                return enif_make_badarg(env);
            }
            else{
                d.content[i] = a.content[i] / b.content[i];
            }
        }
        return matrix_to_erl(env, d);
    }

    return enif_make_badarg(env);
}


//Transpose of a matrix.
Matrix tr(Matrix a){
    Matrix result = matrix_alloc(a.n_cols ,a.n_rows);

    for(int j = 0; j < a.n_rows; j++){
        for(int i = 0; i < a.n_cols; i++){
            result.content[i*result.n_cols+j] = a.content[j*a.n_cols+i];
        }
    }
    return result;
}
ERL_NIF_TERM nif_transpose(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){

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

    Matrix inv = matrix_dup(a);

    int N = a.n_rows;
    int *IPIV = enif_alloc(sizeof(int)*N);
    int LWORK = N*N;
    double* WORK = enif_alloc(sizeof(int)*N*N);
    int INFO1, INFO2;

    dgetrf_(&N,&N,inv.content,&N,IPIV,&INFO1);
    dgetri_(&N,inv.content,&N,IPIV,WORK,&LWORK,&INFO2);

    
    enif_free(IPIV);
    enif_free(WORK);
    ERL_NIF_TERM result;

    if(INFO1 > 0 || INFO2 > 0){
        result = enif_raise_exception(env, enif_make_atom(env, "nif_inv: could not invert singular matrix."));
        matrix_free(inv);
    }
    else if(INFO1 < 0 || INFO2 < 0){
        result = enif_raise_exception(env, enif_make_atom(env, "nif_inv: LAPACK error."));
        matrix_free(inv);
    }
    else result = matrix_to_erl(env, inv);


    return result;
}

//----------------------------------------------------------------------------------------------------|
//                        ------------------------------------------                                  |
//                        |                   CBLAS                |                                  |
//                        ------------------------------------------                                  |
//----------------------------------------------------------------------------------------------------|
//Some CBLAS wrappers.

//Calculates the norm of input vector/matrix, aka the square root of the sum of its composants.
ERL_NIF_TERM nif_dnrm2(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
    Matrix x;

    if(!enif_get(env, argv, "m", &x)){
        return enif_make_badarg(env);
    }

    double result = cblas_dnrm2(x.n_cols*x.n_rows, x.content, 1);

    return enif_make_double(env, result);
}


//Performs blas_ddot
//Input: two vectors / matrices
ERL_NIF_TERM nif_ddot(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
    Matrix x,y;

    if(!enif_get(env, argv, "mm", &x, &y)){
        return enif_make_badarg(env);
    }

    int n = fmin(x.n_rows*x.n_cols, y.n_rows*y.n_cols);
    if(n <= 0){
        return enif_make_badarg(env);
    }

    double result = cblas_ddot(n, x.content, 1, y.content, 1);

    return enif_make_double(env, result);
}


//Performs blas_daxpy
//Input: a number, vectors X and Y
//Output: a vector of same dimension then Y, containing alpha X + Y
//------------------------------------------------------------------------
ERL_NIF_TERM nif_daxpy(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
    Matrix x,y;
    int n;
    double alpha;

    if(!enif_get(env, argv, "inmm", &n, &alpha, &x, &y)){
        return enif_make_badarg(env);
    }
    
    if(fmin(x.n_rows, x.n_cols) * fmin(y.n_rows, y.n_cols) != 1){
        //We are not using vectors...
        return enif_make_badarg(env);
    }

    Matrix ny = matrix_dup(y);

    cblas_daxpy(n, alpha, x.content, 1, ny.content, 1);

    return matrix_to_erl(env, ny);
}

// Arguments: alpha, A, x, beta, y
// Performs alpha*A*x + beta*y.
// alpha, beta are numbers
// A, x, y are matrices (x and y being vectors)
// x and y are expected to have a length of A.n_cols
ERL_NIF_TERM nif_dgemv(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
    Matrix A,x,y;
    double alpha, beta;

    if(!enif_get(env, argv, "nmmnm", &alpha, &A, &x, &beta, &y)){
        enif_make_badarg(env);
    }

    //Check dimensions compatibility
    int vec_length = fmin(fmax(x.n_cols, x.n_rows), fmax(y.n_cols, y.n_rows));
    if(vec_length < A.n_cols || fmin(x.n_cols, x.n_rows) != 1 || fmin(y.n_cols, y.n_rows) != 1){
        enif_make_badarg(env);
    }

    Matrix ny = matrix_dup(y);

    cblas_dgemv(CblasRowMajor, CblasNoTrans, A.n_rows, A.n_cols, alpha, A.content, A.n_rows, x.content, 1,  beta, ny.content, 1);

    return matrix_to_erl(env, ny);
}



//Arguments: double alpha, matrix A, matrix B, double beta, matrix C
ERL_NIF_TERM nif_dgemm(ErlNifEnv *env, int argc, const ERL_NIF_TERM argv[]){
    Matrix A, B;
    double alpha = 1.0;
    double beta = 0.0;

    if(!enif_get(env, argv, "mm", &A, &B)
            || A.n_cols != B.n_rows){

        return enif_make_badarg(env);
    }
    
    Matrix C = matrix_alloc(A.n_rows, B.n_cols);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    A.n_rows, B.n_cols, A.n_cols,
    alpha, A.content, A.n_cols, B.content, B.n_rows, 
    1.0, C.content, C.n_cols);

    return matrix_to_erl(env, C);
}

ErlNifFunc nif_funcs[] = {
    {"matrix", 1, nif_matrix},
    {"get", 3, nif_get},
    {"at", 2, nif_at},
    {"mtfli", 1, nif_mtfli},
    {"mtfl", 1, nif_mtfl},
    {"equals", 2, nif_equals},
    {"row", 2, nif_row},
    {"col", 2, nif_col},
    {"zeros", 2, nif_zeros},
    {"eye", 1, nif_eye},
    {"mult", 2, nif_mult},
    {"add",  2, nif_add},
    {"sub",  2, nif_sub},
    {"divide", 2, nif_divide},
    {"transpose", 1, nif_transpose},
    {"inv", 1, nif_inv},
    
    //--- BLAS----------
    {"nrm2", 1, nif_dnrm2},
    {"vec_dot", 2, nif_ddot},
    {"dot", 2, nif_dgemm}
};

ERL_NIF_INIT(numerl, nif_funcs, load, NULL, upgrade, NULL)