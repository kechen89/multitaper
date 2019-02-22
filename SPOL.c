#include "DPSS.h"

void SPOL(int N, double *V, int K)
{
    /*
     Scale the discrete prolate spheroidal sequence and use Splepian's convention
     
     Bell, B., Percival, D.B. and Walden, A.T., Calculating Thomson's Spectral
     Multitapers by Inverse Iteration, J. Comput. and Graph. Stat., 1993.
     
     INPUT:
         N      length of dpss sequence
         V      eigenvector dpss with unit energy
         K      order of dpss 0=<K<=N-1
     
     OUTPUT:
         V      dpss conforming to Slepian's polarity convention
     */
    
    int L;
    double DSUM, DWSUM;
    
    DSUM = 0.0;
    DWSUM = 0.0;
    
    for (L = 1; L <= N; L++)
    {
        DSUM += V[L-1];
        DWSUM += V[L-1]*(N - 1.0 - 2.0*(L-1));
        
        if ( ((K%2) == 0 && (DSUM < 0.0)) || ((K%2) == 1 && (DWSUM < 0.0)))
        {
            for (L = 1; L <=N; L++)
            {
                V[L-1] = -V[L-1];
            }
        }
    }
}
