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
    
    int L, J, K, I, L1, L2,L3;
    double V, D, Q, HOLD;
    
    V = R[0];
    F[0] = G[0] / V;
    
    D = R[1];
    W[0] = 1.0;
    Q = F[0] * R[1];
    
    for (L = 2; L <= N; L++)
    {
        W[L-1] = -D/V;
        
        if (L > 2)
        {
            L1 = (L - 2)/2;
            L2 = L1 + 1;
            
            if (L != 3)
            {
                for (J = 2; J <= L2; J++)
                {
                    HOLD = W[J-1];
                    K = L - J + 1;
                    W[J-1] = W[J-1] + W[L-1]*W[K-1];
                    W[K-1] = W[K-1] + W[L-1]*HOLD;
                }
            }//end if
            
            if ((2*L1 != L-2) || L == 3) W[L2] = W[L2] + W[L-1]*W[L2];
        
        }//end if
        
        V = V + W[L-1]*D;
        F[L-1] = (G[L-1] - Q)/V;
        L3 = L-1;
        for (J = 1; J <= L3; J++)
        {
            K = L - J + 1;
            F[J-1] = F[J-1] + F[L-1]*W[K-1];
        }
        
        if (L < N)
        {
        D = 0.0;
        Q = 0.0;
        
        for (I = 1; I <= L; I++)
        {
            K = L - I + 2;
            D = D + W[I-1]*R[K-1];
            Q = Q + F[I-1]*R[K-1];
        }
        
       }//end if
    }

}
