#include <iostream>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string>

#include "NRutil.h"
#include "read_write.h"

#define NR_END 1
#define FREE_ARG char*

double *Nvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;

  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  return v-nl+NR_END;
}

double **Nmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate pointers to rows */
  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

double interpolate (double xData[], double yData[], int size, double x, bool extrapolate) {
  if ((x > xData[size-1]) || (x < xData[0])) {
    std::cout << "xmin: " << xData[0] << " | xmax: " << xData[size-1] << " | x: " << x << std::endl;
    throw_error("ERROR: x out of range called in interpolate.");
  }
  int i = 0;
  if (x >= xData[size - 2]) {
    i = size - 2;
  } else {
    while (x > xData[i+1]) i++;
  }
  double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];
  if (!extrapolate) {
    if (x < xL) yR = yL;
    if (x > xR) yL = yR;
  }
  double dydx = ( yR - yL ) / ( xR - xL );
  return yL + dydx * ( x - xL );
}
