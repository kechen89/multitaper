#include "DPSS.h"

int main (int argc, char **argv)
{
    int NMAX;
    int KMAX;
    int N, I, J;
    double W;
    double **V;
    
    //initialize parameters
    N = 64;
    NMAX = 64;
    KMAX = 3;
    W = 1.0/16.0;
    
    V = alloc2double(KMAX + 1, NMAX);
    
    for (I = 0; I < NMAX; I++)
    {
        for (J = 0; J < KMAX+1; J++)
        {
            V[J][I] = 0.0;
        }
    }
    
    //compute DPSS
    DPSS(NMAX, KMAX, N, W, V);
    
    //save array
    FILE *file1;
    file1 = fopen("DPSS.bin","wb");
    
    for (I = 0; I < KMAX + 1; I++)
    {
        fwrite(V[I], sizeof(double), NMAX, file1);
    }
    fclose(file1);
}
