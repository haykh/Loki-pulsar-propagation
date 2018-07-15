#include <iostream>
#include <math.h>
#include "dsolve.h"
#include "NRutil.h"
#include "rk4.h"
#include "rkf45.h"

void dsolve (double ystart[], double x1, double x2, void (*derivs)(double, double *, double *)) {
  odeint(ystart, 2, x1, x2, 0.01, 1e-4, 1e-6, 0, 0, derivs);
}

void dsolve_r4_rkf45 (double ystart[], double x1, double x2, void (*derivs)(float, float *, float *)) {
  float abserr;
  int flag;
  int n_step;
  float relerr;
  float t;
  float t_out;
  float t_start;
  float t_stop;
  float *y;
  float *yp;

  y = new float[2];
  yp = new float[2];

  abserr = sqrt ( r4_epsilon ( ) );
  relerr = sqrt ( r4_epsilon ( ) );

  flag = 1;

  t_start = (float)(x1);
  t_stop = (float)(x2);

  n_step = 500;

  y[0] = (float)(ystart[0]);
  y[1] = (float)(ystart[1]);

  for (int i_step = 1; i_step <= n_step; i_step++ ) {
    t = ( ( float ) ( n_step - i_step + 1 ) * t_start
        + ( float ) (          i_step - 1 ) * t_stop )
        / ( float ) ( n_step              );

    t_out = ( ( float ) ( n_step - i_step ) * t_start
            + ( float ) (	   i_step ) * t_stop )
            / ( float ) ( n_step );

    flag = r4_rkf45 ( derivs, 2, y, yp, &t, t_out, &relerr, abserr, flag );
  }
  ystart[0] = (double)(y[0]);
  ystart[1] = (double)(y[1]);

  delete [] y;
  delete [] yp;
}

void dsolve_r8_rkf45 (double ystart[], double x1, double x2, void (*derivs)(double, double *, double *)) {
  double abserr;
  int flag;
  int i_step;
  int n_step;
  double relerr;
  double t;
  double t_out;
  double t_start;
  double t_stop;

  double *y;
  double *yp;

  y = new double[2];
  yp = new double[2];

  abserr = sqrt ( r8_epsilon ( ) );
  relerr = sqrt ( r8_epsilon ( ) );

  flag = 1;

  t_start = x1;
  t_stop = x2;

  n_step = 500;

  y[0] = ystart[0];
  y[1] = ystart[1];
  derivs(x1, y, yp);

  for (int i_step = 1; i_step <= n_step; i_step++) {
    t = ( ( double ) ( n_step - i_step + 1 ) * t_start
        + ( double ) (          i_step - 1 ) * t_stop )
        / ( double ) ( n_step              );

    t_out = ( ( double ) ( n_step - i_step ) * t_start
            + ( double ) (	   i_step ) * t_stop )
            / ( double ) ( n_step );
    // std::cout << i_step << " " << t << " " << t_out << " " << flag << "\n";
    flag = r8_rkf45 ( derivs, 2, y, yp, &t, t_out, &relerr, abserr, flag );
  }
  ystart[0] = y[0];
  ystart[1] = y[1];

  delete [] y;
  delete [] yp;
}
