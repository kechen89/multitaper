#include "DPSS.h"

void DPSS(int NMAX, int KMAX, int N, float W, double **V)
{
/*
 Calculate Discrete Prolate Spheroidal Sequences (DPSS)
 
 Bell, B., Percival, D.B. and Walden, A.T., Calculating Thomson's Spectral
 Multitapers by Inverse Iteration, J. Comput. and Graph. Stat., 1993.
 
 INPUT:
     NMAX      maximum length of taper
     KMAX      maximum order of dpss
     N         length of time series
     W         half-bandwidth W < 1/2
 
 OUTPUT:
     V         columns contain tapers
 */
    int IT;
    int MAXIT;
    int J, JJ, K, M, ISIG, ILOW, IHIG;
    double *SINES, *VOLD, *U, *SCR1, *SIG;
    
    SINES = alloc1double(N);
    VOLD = alloc1double(N);
    U = alloc1double(N);
    SCR1 = alloc1double(N);
    SIG = alloc1double(KMAX + 1);
    
    
    //Check input parameters
    if (W > 0.5)
    {
        fprintf(stderr, "The bandwidth is larger than 0.5\n");
        exit(0);
    }
    
    if (N < 2)
    {
        fprintf(stderr, "Number of samples is too small\n");
        exit(0);
    }
    
    if (NMAX < N)
    {
        fprintf(stderr, "Matrix too small\n");
        exit(0);
    }
    
    if (KMAX < 0 || KMAX > N - 1)
    {
        fprintf(stderr, "KMAX < 0\n");
        exit(0);
    }
    
    //Set up S matrix, S(n,m) = SINES(n-m), for n not equal to m
    for (M = 1; M <= N-1; M++)
    {
        SINES[M] = sin(2.0*PI*W*M)/(PI*M);
    }
    
    //loop over DPSS orders 0 to KMAX
    
    //Set SINES[0]
    for (K = 0; K <= KMAX; K++)
    {
        if (K == 0)
        {
            SINES[0] = 2.0 * W - 1.0;
        }
        else
        {
            SINES[0] = 2.0 * W - (1.0 + SIG[K - 1]);
        }
        
        // Define starting vector for inverse iteration
        ISIG = 1;
        
        for (J = 0; J <= K; J++)
        {
            ILOW = J * N / (K + 1);
            IHIG = (J + 1) * N / (K + 1) - 1;
            
            for (JJ = ILOW; JJ <= IHIG; JJ++)
            {
                U[JJ] = ISIG * (1.0 / sqrt(N));
            }
            
            ISIG = -1 * ISIG;
        }
        
        if ((K % 2) > 0 && (N % 2) > 0) U[N/2 + 1] = 0.0;    //K is odd, N is odd
        
        /*Verified*/
        
        //Maximum number of iterations
        MAXIT = (K + 3) * sqrt(N);
        
        //Carry out inverse iteration
        for (IT = 1; IT <= MAXIT; IT++)
        {
            
            //copy U into old V
            for (J = 0; J < N; J++)
            {
                VOLD[J] = U[J];
            }
            
            //Solve symmetric Toeplitz matrix
            SYTOEP(N, SINES, VOLD, U, SCR1);
            
            //new vector must be orthogonal to previous eigenvectors
            if (K > 0)
            {
                for (K1 = 0; K1 < K; K1++)
                {
                    //projection of U onto V[*][K1]
                    PROJ = 0.0;
                    for (J = 0; J < N; J++)
                    {
                        PROJ = PROJ + U[J]*V[J][K1];
                    }
                    //subtract projection
                    for (J = 0; J < N; J++)
                    {
                        U[J] = U[J] - PROJ*V[J][K1];
                    }
                } //end for K1
            } //end if
            
            //normalize
            SNORM = 0.0;
            
            for (J = 0; J < N; J++)
            {
                SNORM += U[J]*U[J];
            }
            
            SSNORM = sqrt(SNORM);
            
            for (J = 0; J < N; J++)
            {
                U[J] = U[J]/SSNORM;
            }

            //Verified
            
            //check for convergence
            SUM = 0.0;
            DIFF = 0.0;
            
            for (J = 0; J < N; J++)
            {
                DIFF += pow(VOLD[J] - U[J], 2);
                SUM += pow(VOLD[J] + U[J], 2);
            }
            
            DELTA = sqrt(MIN(DIFF,SUM));
        
            if (DELTA <= EPS) break;
            
        }//end for inverse iteration
        
        if (SUM < DIFF)
        {
            if (K == 0)
            {
                SIG[0]= - 1.0/SSNORM;
                
            }
            else
            {
                SIG[K] = SIG[K-1] - 1.0/SSNORM;
                
            }
            
        }
        else
        {
            if(K == 0)
            {
                SIG[0]= 1.0/SSNORM;
                
            }
            else
            {
                SIG(K) = SIG(K-1) + 1.0/SSNORM;
                
            }
            
        }
        
        //ensure tapers satisfy Slepian convention
        SPOL(N, U, K);
        
        for (J = 0; J < N; J++)
        {
            V[J][K] = U[J];
        }
        
    } //end for K
}
