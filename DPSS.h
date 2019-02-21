/*--------------------------------Macros-------------------------------------------------------*/
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define PI 3.141592653589793238462643

/*----------------------------include headers----------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double *alloc1double(int n);

double **alloc2double(int n1, int n2);

void DPSS(int NMAX, int KMAX, int N, float W, double **V);
