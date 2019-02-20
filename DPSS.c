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
    int IFAULT, TOTIT;
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
    
    //Set up S matrix
    for (M = 1; M <= N-1; M++)
    {
        SINES[M] = sin(2*PI*W*M)/(PI*M);
    }
    
    //Set total iteration counter
    TOTIT = 0;
    
    for (K = 0; K <= KMAX; K++)
    {
        if (K == 0)
        {
            SINES[0] = 2 * W - 1;
        }
        else
        {
            SINES[0] = 2 * W - (1 + SIG[K - 1]);
        }
        
        // Define starting vector for inverse iteration
        ISIG = 1;
        for (J = 0; J <= K; J++)
        {
            ILOW = J * N / (K + 1);
            IHIG = (J + 1) * N / (K + 1) - 1;
            for (JJ = ILOW; JJ <= IHIG; JJ++)
            {
                U[JJ] = ISIG * (1 / sqrt(N));
            }
            ISIG = -1 * ISIG;
        }
        if ((K % 2) > 0 && (N % 2) > 0) U[(N/2) + 1] = 0.0;
        
        //Maximum number of iterations
        MAXIT = (K + 3) * sqrt(N);
        
        //Carry out inverse iteration
        for (IT = )
        
    }
    
}
