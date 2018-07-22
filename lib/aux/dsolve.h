#pragma DSOLVE

void dsolve (double ystart[], double x1, double x2, void (*derivs)(double, double *, double *));
void dsolve_r4_rkf45 (double ystart[], double x1, double x2, void (*derivs)(float, float *, float *));

void dsolve_r8_rkf45 (double ystart[], double x1, double x2, void (*derivs)(double, double *, double *));
