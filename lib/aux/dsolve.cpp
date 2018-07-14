#include <iostream>
#include "dsolve.h"
#include "NRutil.h"
#include "rk4.h"

void dsolve (double ystart[], double x1, double x2, void (*derivs)(double, double *, double *)) {
  odeint(ystart, 2, x1, x2, 0.01, 1e-4, 1e-6, 0, 0, derivs);
}
