#include "DPSS.h"

int main (int argc, char **argv)
{
    int NMAX;
    int KMAX;
    int N;
    float W;
    double **V;
    
    //initialize parameters
    N = 64;
    NMAX = 64;
    KMAX = 3;
    W = 1.0/16.0;
    
    V = alloc2double(NMAX, KMAX + 1);
    
    //compute DPSS
    DPSS(NMAX, KMAX, N, W, V);
    
    //save array
    
}
