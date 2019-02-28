#include "DPSS.h"

/* allocate a 1-d double array*/
double *alloc1double(int n)
{
    double *p;
    p = malloc(n*sizeof(double));
    return p;
}

/* allocate a 2-d double array */
double **alloc2double (int n1, int n2)
{
    double **p;
    int i;
    p = (double**)malloc(n1*sizeof(double*));
    for (i = 0; i < n1; i++)
    {
        p[i] = (double*)malloc(n2*sizeof(double));
    }
	return p;
}

