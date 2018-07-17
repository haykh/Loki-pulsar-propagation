#include <iostream>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string>
#include <math.h>


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

double trilinear(Point3 *p, double *d, int xsize, int ysize, int zsize, double def) {
#define DENS(X, Y, Z) d[(X)+xsize*((Y)+ysize*(Z))]

/* linear interpolation from l (when a=0) to h (when a=1)*/
/* (equal to (a*h)+((1-a)*l) */
#define LERP(a,l,h)	((l)+(((h)-(l))*(a)))

    int	       x0, y0, z0,
	       x1, y1, z1;
    double     *dp,
	       fx, fy, fz,
	       d000, d001, d010, d011,
	       d100, d101, d110, d111,
	       dx00, dx01, dx10, dx11,
	       dxy0, dxy1, dxyz;

    x0 = (int)floor(p->x); fx = p->x - x0;
    y0 = (int)floor(p->y); fy = p->y - y0;
    z0 = (int)floor(p->z); fz = p->z - z0;

    x1 = x0 + 1;
    y1 = y0 + 1;
    z1 = z0 + 1;

    if (x0 >= 0 && x1 < xsize &&
	y0 >= 0 && y1 < ysize &&
	z0 >= 0 && z1 < zsize)
    {
	dp = &DENS(x0, y0, z0);
	d000 = dp[0];
	d100 = dp[1];
	dp += xsize;
	d010 = dp[0];
	d110 = dp[1];
	dp += xsize*ysize;
	d011 = dp[0];
	d111 = dp[1];
	dp -= xsize;
	d001 = dp[0];
	d101 = dp[1];
    }
    else
    {
#	define INRANGE(X, Y, Z) \
		  ((X) >= 0 && (X) < xsize && \
		   (Y) >= 0 && (Y) < ysize && \
		   (Z) >= 0 && (Z) < zsize)

	d000 = INRANGE(x0, y0, z0) ? DENS(x0, y0, z0) : def;
	d001 = INRANGE(x0, y0, z1) ? DENS(x0, y0, z1) : def;
	d010 = INRANGE(x0, y1, z0) ? DENS(x0, y1, z0) : def;
	d011 = INRANGE(x0, y1, z1) ? DENS(x0, y1, z1) : def;

	d100 = INRANGE(x1, y0, z0) ? DENS(x1, y0, z0) : def;
	d101 = INRANGE(x1, y0, z1) ? DENS(x1, y0, z1) : def;
	d110 = INRANGE(x1, y1, z0) ? DENS(x1, y1, z0) : def;
	d111 = INRANGE(x1, y1, z1) ? DENS(x1, y1, z1) : def;
    }

    dx00 = LERP(fx, d000, d100);
    dx01 = LERP(fx, d001, d101);
    dx10 = LERP(fx, d010, d110);
    dx11 = LERP(fx, d011, d111);

    dxy0 = LERP(fy, dx00, dx10);
    dxy1 = LERP(fy, dx01, dx11);

    dxyz = LERP(fz, dxy0, dxy1);

    return dxyz;
}
