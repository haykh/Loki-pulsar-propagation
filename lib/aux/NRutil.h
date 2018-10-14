#pragma NRUTIL

double *Nvector(long nl, long nh);
double **Nmatrix(long nrl, long nrh, long ncl, long nch);
double interpolate (double xData[], double yData[], int size, double x, bool extrapolate);

typedef struct Point3Struct {	/* 3d point */
	double x, y, z;
	} Point3;
typedef Point3 Vector3;

double trilinear(Point3 *p, double *d, int xsize, int ysize, int zsize, double def);
