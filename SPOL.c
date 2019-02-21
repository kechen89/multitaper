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
    
    DSUM = 0.0;
    DWSUN = 0.0;
    
    for (L = 0; L < N; L++)
    {
        DSUM += V[L];
        DWSUM += V[L]*(N - 1.0 - 2.0*(L));
        
        if ( ((K%2) == 0 && (DSUM < 0.0)) || ((K%2) == 1) && (DSUM < 0.0))
        {
            for (L = 0; L < N; L++)
            {
                V[L] = -V[L];
            }
        }
    }
    
    return;
}
