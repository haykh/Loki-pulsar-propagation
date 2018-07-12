#pragma NRUTIL

double *Nvector(long nl, long nh);
double **Nmatrix(long nrl, long nrh, long ncl, long nch);
double interpolate (double xData[], double yData[], int size, double x, bool extrapolate);
