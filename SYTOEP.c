#include "DPSS.h"

void SYTOEP(int N, double *R, double *G, double *F, double *W)
{
    /*
     Calculate filter corresponding to a symmetric Toeplitz matrix
     
     Bell, B., Percival, D.B. and Walden, A.T., Calculating Thomson's Spectral
     Multitapers by Inverse Iteration, J. Comput. and Graph. Stat., 1993.
     
     INPUT:
         N       dimension of Toeplitz matrix and cross-correlation vector
         R       autocovariance vector from lag 0 to N-1
         G       cross-corrleation vector
         W       work array
     
     OUTPUT:
         F       filter
     */
    
    int L, J;
    double V, D, Q;
    
    V = R[0];
    F[0] = G[0] / V;
    
    D = R[1];
    W[0] = 1.0；
    Q = F[0] * R[1];
    
    for (L = 1; L < N; L++)
    {
        W[L] = -D/V;
        
        if (L > 1)
        {
            L1 = (L - 1)/2;
            L2 = L1 + 1;
            
            if (L != 2)
            {
                for (J = 1; J < L2; J++)
                {
                    HOLD = W[J];
                    K = L - J;
                    W[J] = W[J] + W[L]*W[K];
                    W[K] = W[K] + W[L]*HOLD;
                }
            }//end if
            
            if ((2*L1 != L-1) || L == 2) W[L2] = W[L2] + W[L]*W[L2];
        
        }//end if
        
        V = V + W[L]*D;
        F[L] = (G[L] - Q)/V;
        
        for (J = 0; J < L; J++)
        {
            K = L - J;
            F[J] = F[J] + F[L]*W[K];
        }
        
        if (L == N - 1) exit(0);
        
        D = 0.0;
        Q = 0.0;
        
        for (I = 0; I < L; I++)
        {
            K = L - I + 1;
            D = D + W[I]*R[K];
            Q = Q + F[I]*R[K];
        }
    }
    return;
}
