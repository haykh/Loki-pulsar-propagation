#pragma PROCESS_FUNCTIONS
using namespace std;

vector <double> vB(double R);
vector <double> vBsplit(double R);

double delta (double R);
double BetaB (double R);
double Q (double R);
double Lambda (double R);

double dtau (double R);
double psi_m (double R);
vector <double> vR (double R);
double gFunc (double R);

vector <double> vMoment (double R);
vector <double> vUdr (double R);
double theta_kb (double R);

/* NEW STUFF HERE /> */

double r_perpFromR(double R1, double R2);
vector <double> vBdipoleXYZ(vector <double> r, vector <double> m);
vector <double> vBsplitXYZ (vector <double> r);
vector <double> Bxyz(vector <double> r, vector <double> m);
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
void CreateMas(int N);
void DelMas();
