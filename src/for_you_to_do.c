#include "../include/for_you_to_do.h"
#include <math.h>

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 126;
  
}

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf(double *A, int *ipiv, int n) 
{
    // add your code here
    int i, j, k, max_index, tmp2;
    double max, tmp1;
    // used for swapping rows
    double * tmp_row = (double *)malloc(sizeof(double) * n);
    
    for (i = 0; i < n; i++)
    {
        max_index = i;
        max = fabs(A[i * n + i]);
        for (j = i + 1; j < n; j++)
        {
            tmp1 = fabs(A[j * n + i]);
            if (max < tmp1)
            {
                max_index = j;
                max = tmp1;
            }
        }
        
        // if the matrix is singular
        if (max == 0)
        {
	    perror("LU factorization failed: coefficient matrix is singular.\n");
            return -1;
        }
           
        if (max_index != i)
        {	
            tmp2 = ipiv[i];
            ipiv[i] = ipiv[max_index];
            ipiv[max_index] = tmp2;
		
            // swap rows of A
            memcpy(tmp_row, A + i * n, sizeof(double) * n);
            memcpy(A + i * n, A + max_index * n, sizeof(double) * n);
            memcpy(A + max_index * n, tmp_row, sizeof(double) * n);
        }
    
        for (j = i + 1; j < n; j++)
        {
            A[j * n + i] = A[j * n + i] / A[i * n + i];
            for (k = i + 1; k < n; k++)
            {
                A[j * n + k] -= A[j * n + i] * A[i * n + k];
            }
        }
    }
    free(tmp_row);
    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    // add your code here
    int i, j;
    double sum = 0;
    double * tmp_B = (double *)malloc(sizeof(double) * n);
    
    if (UPLO == 'L')
    {
	// Ly = b, forward subtitution    
        for (i = 0; i < n; i++)
        {
            tmp_B[i] = B[ipiv[i]];
        }
        
        for (i = 0; i < n; i++)
        {
            sum = tmp_B[i];
            for (j = 0; j < i; j++)
            {
                sum -= B[j] * A[i * n + j];
            }
            B[i] = sum;
        }
    }
    else if (UPLO == 'U')
    {
	// Ux = y, backward subtitution     
        for (i = n - 1; i >= 0; i--)
        {
            sum = 0;
            for (j = i + 1; j < n; j++)
            {
                sum += B[j] * A[i * n + j];
            }
            B[i] = (B[i] - sum) / A[i * n + i];
        }
    }
    free(tmp_B);
    return;
}



/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{    
    register int i1, j1, k1; 
    // used for optimizing	  
    register int in1, in2, in3, j2, j3, kn;	
    register int n1 = i + b > n ? n : i + b;
    register int n2 = j + b > n ? n : j + b;
    register int n3 = k + b > n ? n : k + b;

    for (i1 = i; i1 < n1; i1 += 3)
    {
        for (j1 = j; j1 < n2; j1 += 3)
        {
	    in1 = i1 * n;
	    in2 = (i1 + 1) * n;
	    in3 = (i1 + 2) * n;
	    j2 = j1 + 1;
	    j3 = j1 + 2;	
            register double C_0_0 = C[in1 + j1];
            register double C_0_1 = C[in1 + j2];
            register double C_0_2 = C[in1 + j3];
            register double C_1_0 = C[in2 + j1];
            register double C_1_1 = C[in2 + j2];
            register double C_1_2 = C[in2 + j3];
            register double C_2_0 = C[in3 + j1];
            register double C_2_1 = C[in3 + j2];
            register double C_2_2 = C[in3 + j3];
		
	    for (k1 = k; k1 < n3; k1++)
            {
		kn = k1 * n + j1;    
                register double A_0 = A[in1 + k1];
                register double A_1 = A[in2 + k1];
                register double A_2 = A[in3 + k1];
                register double B_0 = B[kn];
                register double B_1 = B[kn + 1];
                register double B_2 = B[kn + 2];

                C_0_0 -= A_0 * B_0;
                C_0_1 -= A_0 * B_1;
                C_0_2 -= A_0 * B_2;
                C_1_0 -= A_1 * B_0;
                C_1_1 -= A_1 * B_1;
                C_1_2 -= A_1 * B_2;
                C_2_0 -= A_2 * B_0;
                C_2_1 -= A_2 * B_1;
                C_2_2 -= A_2 * B_2;
            }
		
            C[in1 + j1] = C_0_0;
            C[in1 + j2] = C_0_1;
            C[in1 + j3] = C_0_2;
            C[in2 + j1] = C_1_0;
            C[in2 + j2] = C_1_1;
            C[in2 + j3] = C_1_2;    
            C[in3 + j1] = C_2_0;
            C[in3 + j2] = C_2_1;
            C[in3 + j3] = C_2_2;
    	}	
    }
    return;
}


/**
 * 
 * this function computes LU decomposition
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *    
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *     
 *      b                , block size   
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
    register int ib, i, j, k, max_index, tmp2;
    register double max, tmp1, sum;	
    // used for optimizing
    register int ib2, n1, in, maxn, jn; 	
    register int size_n = sizeof(double) * n;
    // used for swapping rows	
    register double * tmp_row = (double *)malloc(size_n);

    // LU factorization		
    for (ib = 0; ib < n; ib += b)
    {
	ib2 = ib + b; 
	n1 = ib2 > n ? n : ib2;   
        for (i = ib; i < n1; i++)
        {
            max_index = i;
	    in = i * n;
            max = fabs(A[in + i]);

     	    // Pivoting       
            for (j = i + 1; j < n; j++)
            {
		tmp1 = fabs(A[j * n + i]);    
                if (tmp1 > max)
                {
                    max_index = j;
                    max = tmp1;
                }
            }
            
            // singularity check
            if (max == 0)
            {
		perror("LU factorization failed: Singular Matrix.\n");
                return -1;
            }
            else
            {    
                if (max_index != i)
                {
                    tmp2 = ipiv[i];
                    ipiv[i] = ipiv[max_index];
                    ipiv[max_index] = tmp2;    	
                    // swapping rows
		    maxn = max_index * n;	
                    memcpy(tmp_row, A + in, size_n);
                    memcpy(A + in, A + maxn, size_n);
                    memcpy(A + maxn, tmp_row, size_n);
                }
            }
		
	    // Gauss Elimination 
            for (j = i + 1; j < n; j++)
            {
		jn = j * n;    
                A[jn + i] = A[jn + i] / A[in + i];
                  
                for (k = i + 1; k < n1; k++)
                {
                    A[jn + k] -= A[jn + i] * A[in + k];
                }
            }
        }

        // A(ib:end, end+1:n) = LL-1 * A(ib:end, end+1:n)
        for (i = ib; i < n1; i++)
        {
	    in = i * n;	
            for (j = ib2; j < n; j++)
            {
                sum = 0;
                for (k = ib; k < i; k++)
                {
                    sum += A[in + k] * A[k * n + j];
                }
                A[in + j] -= sum;
            }
        }

        // A(end+1:n, end+1:n) -= A(end+1:n, ib:end) * A(ib:end, end+1:n)
        for (i = ib2; i < n; i += b)
        {
            for (j = ib2; j < n; j += b)
            {
                mydgemm(A, A, A, n, i, j, ib, b);
            }
        }
    }
    free(tmp_row);
    return 0;
}
