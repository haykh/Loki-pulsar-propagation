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

double r_perpFromR(double R1, double R2);
vector <double> Bxyz(vector <double> r);
vector <double> M (vector <double> r);
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
void CreateMas(int N);
void DelMas();