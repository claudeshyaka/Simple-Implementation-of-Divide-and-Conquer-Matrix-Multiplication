/**
 * Copyright (c) 2018 I-Ting Angelina Lee  
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 **/

// This switch allows for the quick and easy switching between the serial eliason and the
// parallel implementation with cilk runtime.
#if 1
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#else
#define cilk_spawn 
#define cilk_sync
#endif

#include <numa.h>

// This controls the debug infrastructure added to the code base. The debugPrintf and
// print_mm macros are utilized to support the easy enable/disable of debug
// infrastructure.
#if 0
#define debugPrintf printf
#define print_mm print_matrix
#else
#define debugPrintf(...)
#define print_mm(...)
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "getoptions.h"
#include "ktiming.h"
#include "papi.h"

#ifndef RAND_MAX
#define RAND_MAX 32767
#endif

#define REAL double 

#define EPSILON (1.0E-6)
#define BASE_BLOCK_SIZE 8
#define Z_WIDTH 16

unsigned long rand_nxt = 0;

int cilk_rand(void) {
    int result;
    rand_nxt = rand_nxt * 1103515245 + 12345;
    result = (rand_nxt >> 16) % ((unsigned int) RAND_MAX + 1);
    return result;
}

/*
 * Naive sequential algorithm, for comparison purposes
 */
void matrixmul(REAL *C, REAL *A, REAL *B, int n) {

  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      REAL s = (REAL)0;
      for(int k = 0; k < n; ++k) {
        s += A[i*n+k] * B[k*n+j];
      }
      C[i*n+j] = s;
    }
  }
}

void print_matrix( REAL *A, int n, int orig_n )
{
    for (int i = 0; i < n; ++i )
    {
        for( int j = 0; j < n; ++j )
        {
            printf("%16.2f ", A[i * orig_n + j]);
        }
        printf("\n");
    }
}

void mm_base( REAL *C, REAL *A, REAL *B, int block_size, int original_size )
{
    REAL s = (REAL)0;
    int a_index = 0;
    int b_index = 0;

    
    debugPrintf("\n\nmm_base; a_addr=0x%08x, b_addr=0x%08x\n", A, B);

    for( int i = 0; i < block_size; ++i )
    {
        // a_index = i * original_size;
        for( int j = 0; j < block_size; ++j )
        {
            s = (REAL)0;
            b_index = 0;
            for( int k = 0; k < block_size; ++k )
            {
                debugPrintf("%f * %f\n", A[ a_index + k ],B[ b_index + j ]);
                s += A[ a_index + k ] * B[ b_index + j ];
                b_index += original_size;
            }
            C[a_index + j] += s;
        }
        a_index += original_size;
    }
}

void mm_morton_base( REAL *C, REAL *A, REAL *B, int block_size, int original_size )
{
    int a_index = 0;
    int b_index = 0;
    REAL s = 0;

    // for every row in matrix A
    for( int i = 0; i < block_size; ++i )
    {
        // for each column in matrix B
        b_index = 0;
        for( int j = 0; j < block_size; ++j )
        {
            s = (REAL)0;
            for( int k = 0; k < block_size; ++k )
            {
                debugPrintf("%f * %f\n", A[ a_index + k ], B[ b_index + k ]);
                s += A[ a_index + k ] * B[ b_index + k ];
            }
            C[a_index + j] += s;
            b_index += block_size;
        }
        a_index += block_size;
    }

    debugPrintf("\n\nmm_morton_base: c=0x%08x, block_size=%d\n", C, block_size);
    print_mm( C, block_size, block_size);
}

/*
 * Compare two matrices.  Print an error message if they differ by
 * more than EPSILON.
 * return 0 if no error, and return non-zero if error.
 */
static int compare_matrix(REAL *A, REAL *B, int n) {

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      // compute the relative error c 
      REAL c = A[n*i + j] - B[n*i + j];
      if (c < 0.0) c = -c;

      c = c / A[n*i + j];
      if (c > EPSILON) { 
        return -1; 
      }
    }
  }
  return 0;
}

//init the matrix to all zeros 
void zero(REAL *M, int n){
    int i;
    for(i = 0; i < n * n; i++) {
        M[i] = 0.0;
    }
}

//init the matrix to random numbers
void init(REAL *M, int n){
    int i,j;
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++){
            M[i*n+j] = (REAL) cilk_rand();
        }
    }
}

//recursive parallel solution to matrix multiplication - row major order
void mat_mul_par(REAL *A, REAL *B, REAL *C, int n, int orig_n) {

    if(n <= BASE_BLOCK_SIZE) {
        mm_morton_base( C, A, B, n, orig_n );
        // mm_base( C, A, B, n, orig_n );
        return;
    }

    int sub_block_length = n * n;
    int half = sub_block_length >> 1;
    int fourth = sub_block_length >> 2;

    //partition each matrix into 4 sub matrices
    //each sub-matrix points to the start of the z pattern
    REAL *A1 = &A[0];
    REAL *A2 = &A[fourth];
    REAL *A3 = &A[half];
    REAL *A4 = &A3[fourth];

    REAL *B1 = &B[0];
    REAL *B2 = &B[half];
    REAL *B3 = &B[fourth];
    REAL *B4 = &B2[fourth];

    REAL *C1 = &C[0];
    REAL *C2 = &C[fourth];
    REAL *C3 = &C[half];
    REAL *C4 = &C3[fourth];

    //recrusively call the sub-matrices for evaluation in parallel
    cilk_spawn mat_mul_par(A1, B1, C1, n >> 1, orig_n);
    cilk_spawn mat_mul_par(A1, B2, C2, n >> 1, orig_n);
    cilk_spawn mat_mul_par(A3, B1, C3, n >> 1, orig_n);
    mat_mul_par(A3, B2, C4, n >> 1, orig_n);
    cilk_sync; //wait here for first round to finish

    cilk_spawn mat_mul_par(A2, B3, C1, n >> 1, orig_n);
    cilk_spawn mat_mul_par(A2, B4, C2, n >> 1, orig_n);
    cilk_spawn mat_mul_par(A4, B3, C3, n >> 1, orig_n);
    mat_mul_par(A4, B4, C4, n >> 1, orig_n);
    cilk_sync; //wait here for all second round to finish
}

void transformMatrixA( REAL *src, REAL *z_dest, int z_dest_size, int row_index, int col_index, int matrix_width, int src_original_n )
{
    if( matrix_width == BASE_BLOCK_SIZE )
    {
        // reached base case
        int row_idx = 0;
        int col_idx = 0;
        int z_idx = 0;
        int calc_col_idx = 0;
        int calc_row_idx = 0;

        for( row_idx = 0; row_idx < BASE_BLOCK_SIZE; row_idx++ )
        {
            calc_row_idx = row_index + row_idx;
            for( col_idx = 0; col_idx < BASE_BLOCK_SIZE; col_idx++ )
            {
                
                calc_col_idx = col_index + col_idx;
                z_dest[z_idx++] = src[(calc_row_idx * src_original_n) + calc_col_idx];
            }
        }
        return;
    }

    // compute new location within the morton z buffer
    int fourth_size = z_dest_size >> 2;
    REAL *top_left     = z_dest;
    REAL *top_right    = z_dest + (fourth_size);
    REAL *bottom_left  = z_dest + (z_dest_size >> 1);
    REAL *bottom_right = z_dest + 3 * (fourth_size);

    // recursively sub-partition the matrix until the base case is reached
    transformMatrixA( src, top_left,     fourth_size, row_index, col_index,                                             matrix_width >> 1, src_original_n );
    transformMatrixA( src, top_right,    fourth_size, row_index, col_index + (matrix_width >> 1),                       matrix_width >> 1, src_original_n );
    transformMatrixA( src, bottom_left,  fourth_size, row_index + (matrix_width >> 1), col_index,                       matrix_width >> 1, src_original_n );
    transformMatrixA( src, bottom_right, fourth_size, row_index + (matrix_width >> 1), col_index + (matrix_width >> 1), matrix_width >> 1, src_original_n );
}

void transformMatrixB( REAL *src, REAL *z_dest, int z_dest_size, int row_index, int col_index, int matrix_width, int src_original_n )
{
    if( matrix_width == BASE_BLOCK_SIZE )
    {

        // reached base case
        int row_idx = 0;
        int col_idx = 0;
        int z_idx = 0;
        int calc_col_idx = 0;
        int calc_row_idx = 0;

        for( col_idx = 0; col_idx < BASE_BLOCK_SIZE; col_idx++ )
        {
            calc_col_idx = col_index + col_idx;
            for( row_idx = 0; row_idx < BASE_BLOCK_SIZE; row_idx++ )
            {
                calc_row_idx = row_index + row_idx;
                z_dest[z_idx++] = src[(calc_row_idx * src_original_n) + calc_col_idx];
            }
        }
        return;
    }

    // compute new location within the morton z buffer
    int fourth_size = z_dest_size >> 2;
    REAL *top_left     = z_dest;
    REAL *bottom_left  = z_dest + fourth_size;
    REAL *top_right    = z_dest + (z_dest_size >> 1);
    REAL *bottom_right = z_dest + 3 * fourth_size;

    // recursively sub-partition the matrix until the base case is reached
    transformMatrixB( src, top_left,     fourth_size, row_index, col_index,                                             matrix_width >> 1, src_original_n );
    transformMatrixB( src, bottom_left,  fourth_size, row_index + (matrix_width >> 1), col_index,                       matrix_width >> 1, src_original_n );
    transformMatrixB( src, top_right,    fourth_size, row_index, col_index + (matrix_width >> 1),                       matrix_width >> 1, src_original_n );
    transformMatrixB( src, bottom_right, fourth_size, row_index + (matrix_width >> 1), col_index + (matrix_width >> 1), matrix_width >> 1, src_original_n );
}

void extractResults( REAL *results_dest, REAL *morton_results, int morton_results_size, int row_index, int col_index, int matrix_width, int original_n  )
{

    if( matrix_width == BASE_BLOCK_SIZE )
    {
        // reached base case
        int row_idx = 0;
        int col_idx = 0;
        int z_idx = 0;
        int calc_col_idx = 0;
        int calc_row_idx = 0;

        debugPrintf("\n\nextractResults: morton_results=0x%08x, size=%d, width=%d\n", morton_results, morton_results_size, matrix_width);
        print_mm( morton_results, matrix_width, matrix_width);

        for( row_idx = 0; row_idx < BASE_BLOCK_SIZE; row_idx++ )
        {
            calc_row_idx = row_index + row_idx;
            for( col_idx = 0; col_idx < BASE_BLOCK_SIZE; col_idx++ )
            {
                calc_col_idx = col_index + col_idx;
                results_dest[(calc_row_idx * original_n) + calc_col_idx] = morton_results[z_idx++];
            }
        }
        return;
    }

    // compute new location within the morton z buffer
    int fourth_size    = morton_results_size >> 2;
    REAL *top_left     = morton_results;
    REAL *top_right    = morton_results + fourth_size;
    REAL *bottom_left  = morton_results + (morton_results_size >> 1);
    REAL *bottom_right = morton_results + 3 * (fourth_size);

    // recursively sub-partition the matrix until the base case is reached
    extractResults( results_dest, top_left,     fourth_size, row_index, col_index,                                             matrix_width >> 1, original_n );
    extractResults( results_dest, top_right,    fourth_size, row_index, col_index + (matrix_width >> 1),                       matrix_width >> 1, original_n );
    extractResults( results_dest, bottom_left,  fourth_size, row_index + (matrix_width >> 1), col_index,                       matrix_width >> 1, original_n );
    extractResults( results_dest, bottom_right, fourth_size, row_index + (matrix_width >> 1), col_index + (matrix_width >> 1), matrix_width >> 1, original_n );
}

const char *specifiers[] = {"-n", "-c", "-h", 0};
int opt_types[] = {INTARG, BOOLARG, BOOLARG, 0};

int usage(void) {
  fprintf(stderr, 
      "\nUsage: mm_dac [-n #] [-c]\n\n"
      "Multiplies two randomly generated n x n matrices. To check for\n"
      "correctness use -c\n");
  return 1;
}

int main(int argc, char *argv[]) {

    int n = 2048;  
    int verify = 0;  
    int help = 0;

    get_options(argc, argv, specifiers, opt_types, &n, &verify, &help);
    if (help || argc == 1) return usage();

    REAL *A, *B, *C;
    REAL *A_MORTON, *B_MORTON, *C_MORTON;

    // determine if we need to apply the interleaving memory policy
    int numWorkers = __cilkrts_get_nworkers();
    printf("numWorkers=%d\n", numWorkers);
    if(numWorkers > 8)
    {
        printf("Enabling Interleave Page Policy\n");
        numa_set_interleave_mask( numa_all_cpus_ptr );
    }

    A = (REAL *) malloc(n * n * sizeof(REAL)); //source matrix 
    B = (REAL *) malloc(n * n * sizeof(REAL)); //source matrix
    C = (REAL *) malloc(n * n * sizeof(REAL)); //result matrix
    A_MORTON = (REAL *) malloc(n * n * sizeof(REAL)); //source matrix 
    B_MORTON = (REAL *) malloc(n * n * sizeof(REAL)); //source matrix
    C_MORTON = (REAL *) malloc(n * n * sizeof(REAL)); //result matrix
    
    init(A, n); 
    init(B, n);
    zero(C, n);
    zero(C_MORTON, n);

    transformMatrixA( A, A_MORTON, n * n, 0, 0, n, n );
    transformMatrixB( B, B_MORTON, n * n, 0, 0, n, n );
    
    clockmark_t begin_rm = ktiming_getmark(); 
    mat_mul_par(A_MORTON, B_MORTON, C_MORTON, n, n);
    clockmark_t end_rm = ktiming_getmark();

    printf("Elapsed time in seconds: %f\n", ktiming_diff_sec(&begin_rm, &end_rm));

    extractResults( C, C_MORTON, n * n, 0, 0, n, n );

    debugPrintf("\n\nComputed Results:\n");
    print_mm( C, n, n);

    if(verify) {
        printf("Checking results ... \n");
        REAL *C2 = (REAL *) malloc(n * n * sizeof(REAL));
        matrixmul(C2, A, B, n);
        verify = compare_matrix(C, C2, n);

        debugPrintf("\n\nCorrect Results:\n");
        print_mm( C2, n, n);
        
        free(C2);
    }

    if(verify) {
        printf("WRONG RESULT!\n");
    } else {
        printf("\nCilk Example: matrix multiplication\n");
        printf("Options: n = %d\n\n", n);
    }

    //clean up memory
    delete [] A;
    delete [] B;
    delete [] C;
    delete [] A_MORTON;
    delete [] B_MORTON;
    delete [] C_MORTON;
	
    return 0;
}
