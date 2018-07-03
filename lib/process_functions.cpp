#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "read_write.h"
#include "constants.h"
#include "functions.h"
#include "process_functions.h"
#include "initialize.h"
using namespace std;

// most of functions should be local

double* r_perp;
double* Rtr;


double sgn (double value) {
  if (value >= 0.0) {
    return 1.0;
  } else {
    return -1.0;
  }
}

vector <double> vMoment (double R) {
  vector <double> mvec(3);
  mvec[0] = sin(Globals::alpha) * cos(Globals::PHI0 + R / Globals::RLC);
  mvec[1] = sin(Globals::alpha) * sin(Globals::PHI0 + R / Globals::RLC);
  mvec[2] = cos(Globals::alpha);
  return mvec;
} // Moment vector
vector <double> vR (double R) {
  vector <double> n0(3);
  n0[0] = sin(Globals::theta_em) * cos(Globals::phi_em);
  n0[1] = sin(Globals::theta_em) * sin(Globals::phi_em);
  n0[2] = cos(Globals::theta_em);

  vector <double> o(3);
  o[0] = sin(Globals::dzeta);
  o[1] = 0.0;
  o[2] = cos(Globals::dzeta);
  return SUM(TIMES(Globals::R_em, n0), TIMES(R, o));
} // Propagation radius vector
double psi_m (double R) {
  return ANGLE(vR(R), vMoment(R));
}

vector <double> vBdipole (double R) {
  vector <double> m;
  vector <double> n;
  m = vMoment (R);
  n = NORMALIZE(vR(R));
  return SUM(TIMES(3.0 * SCALAR(m, n), n), TIMES(-1.0, m));
}
vector <double> vBsplit (double R) {
  double Rr = NORM(vR(R));
  double Rxy = sqrt(vR(R)[0] * vR(R)[0] + vR(R)[1] * vR(R)[1]);

  double costh = SCALAR(NORMALIZE(vR(R)), NORMALIZE(Globals::vOmega));
  double sinth = sqrt(1 - costh * costh);
  double cosphi = vR(R)[0] / Rxy;
  double sinphi = vR(R)[1] / Rxy;
  double phi = acos (cosphi);

  double psi1 = costh * cos(Globals::alpha) + sinth * sin(Globals::alpha) * cos(phi - Globals::PHI0 + Rr / Globals::RLC);

  double Br = (Globals::fr / (Rr * Rr * Globals::RLC)); // * tanh(psi1 / 0.1);
  double Bphi = -(Globals::fphi * Globals::fr * sinth / (Rr * Globals::RLC * Globals::RLC)); // * tanh(psi1 / 0.1);

  Br *= pow(Rr, 3);
  Bphi *= pow(Rr, 3);

  vector <double> temp(3);
  temp[0] = Br * sinth * cosphi - Bphi * sinphi;
  temp[1] = Br * sinth * sinphi + Bphi * cosphi;
  temp[2] = Br * costh;
  return temp;
}

vector <double> vB (double R) {
  if (Globals::fphi == 0 && Globals::fr == 0) {
    return vBdipole(R);
  } else {
    return SUM(vBsplit(R), vBdipole(R));
  }
}
vector <double> vb (double R) {
  return NORMALIZE(vB(R));
}

double theta_kb (double R) {
  vector <double> o(3);
  o[0] = sin(Globals::dzeta);
  o[1] = 0.0;
  o[2] = cos(Globals::dzeta);
  return ANGLE(vB(R), o);
}

vector <double> vBetaR (double R) {
  return TIMES(constants::R_star / constants::c, CROSS(Globals::vOmega, vR(R)));
}
vector <double> vUdr (double R) {
  vector <double> o(3);
  o[0] = sin(Globals::dzeta);
  o[1] = 0.0;
  o[2] = cos(Globals::dzeta);

  vector <double> vn;
  vector <double> vm;
  vn = NORMALIZE(SUM(o, TIMES(-SCALAR(o, vb(R)), vb(R))));
  vm = NORMALIZE(CROSS(vb(R), vn));

  vector <double> temp(3);
  temp[0] = SCALAR(vBetaR(R), vn);
  temp[1] = SCALAR(vBetaR(R), vm);
  if (temp[0] * temp[0] + temp[1] * temp[1] >= 1.0) {
    throw_error("ERROR: vUdr > 1.");
  }
  temp[2] = sqrt(1 - pow(temp[0], 2) - pow(temp[1], 2));
  return temp;
}
double delta (double R) {
  double vx = vUdr(R) [0];
  double vy = vUdr(R) [1];
  double sinth = sin(theta_kb(R));
  double costh = cos(theta_kb(R));
  double sign = sgn (- vy * costh / sqrt(pow(sinth - vx, 2) + pow(costh * vy , 2)));
  return sign * acos((sinth - vx) / sqrt(pow(sinth - vx, 2) + pow(costh * vy , 2)));
}
double BetaB (double R) {
  vector <double> o(3);
  o[0] = sin(Globals::dzeta);
  o[1] = 0.0;
  o[2] = cos(Globals::dzeta);

  vector <double> XX;
  vector <double> YY;
  XX = NORMALIZE(SUM(Globals::vOmega, TIMES(-SCALAR(o, Globals::vOmega), o)));
  YY = CROSS(o, XX);

  double bx = SCALAR(XX, vB(R));
  double by = SCALAR(YY, vB(R));
  // double bb = sqrt (bx * bx + by * by);
  // return acos (bx / bb) * sgn (by);
  return atan(by / bx);
}

double gFunc (double R) {
  double rperp, err;
  polint(Rtr, r_perp, (int) (1.5*Globals::RESCAPE/100.0)+2, R, &rperp, &err);
	double f = pow(rperp / sqrt(Globals::RLC), 2);
	double theta = ANGLE(vR(R), Globals::vOmega);
	double dtheta = 5.0 * constants::PI / 180.0;
	double gap = 1.0;
	if (Globals::alpha_deg > 80)
    gap = (1. - exp(-pow(constants::PI * 0.5 - theta, 2) / (2.0 * dtheta * dtheta)));
	return (pow(f, 2.5) * exp(-f * f) / (pow(f, 2.5) + pow(Globals::f0, 2.5))) * gap;
}
double Ne (double R) {
  double nGJ = SCALAR(Globals::vOmega, vB(R)) * (Globals::B0 / pow(NORM(vR(R)), 3)) / (2 * constants::PI * constants::c * constants::e);
  return Globals::lambda * gFunc (R) * nGJ;
}

double omegaB (double R) {
  return -constants::e * NORM(vB(R)) * (Globals::B0 / pow(NORM(vR(R)), 3)) / (constants::me * constants::c);
}
double omegaW (double R) {
  double vx = vUdr(R) [0];
  double vz = vUdr(R) [2];
  double sinth = sin(theta_kb(R));
  double costh = cos(theta_kb(R));
  return Globals::omega * (1 - sinth * vx - costh * vz);
}
double omegaP (double R) {
  return sqrt(4 * constants::PI * constants::e * constants::e * fabs(Ne(R)) / constants::me);
}

double Q (double R) {
  double vx = vUdr(R) [0];
  double vy = vUdr(R) [1];
  double vz = vUdr(R) [2];
  double sinth = sin(theta_kb(R));
  double costh = cos(theta_kb(R));
  return Globals::lambda * omegaB(R) * Globals::omega * (pow(sinth - vx, 2) + pow(vy * costh, 2)) / (2 * pow(Globals::gamma0, 3) * pow(omegaW(R), 2) * (costh * (1 - vx * vx - vy * vy) - vz * (1.0 - sinth * vx)));
}

double fDist (double gamma) {
  return ((6.0 * Globals::gamma0) / (pow(2.0, 1.0/6.0) * constants::PI)) * (pow(gamma, 4) / (2.0 * pow(gamma, 6) + pow(Globals::gamma0, 6)));
}

double gammaU (double R) {
  double vx = vUdr(R) [0];
  double vy = vUdr(R) [1];
  return pow(1 - vx * vx - vy * vy, -0.5);
}

double INTEGRAL (double gamma, double R) {
  double cA = pow(gammaU(R) * omegaW(R) / omegaB(R), 2);
  return -(pow(2, 2.0 / 3.0)*(-2*sqrt(3)*atan((2*pow(2, 1.0 / 3.0)*pow(gamma,2) - pow(Globals::gamma0,2))/
          (sqrt(3)*pow(Globals::gamma0,2))) - 2*log(pow(2, 1.0 / 3.0)*pow(gamma,2) + pow(Globals::gamma0,2)) +
          log(pow(2, 2.0 / 3.0)*pow(gamma,4) - pow(2, 1.0 / 3.0)*pow(gamma,2)*pow(Globals::gamma0,2) +
          pow(Globals::gamma0,4))) - pow(2, 1.0 / 3.0)*pow(Globals::gamma0,2)*cA*
          (2*sqrt(3)*atan((2*pow(2, 1.0 / 3.0)*pow(gamma,2) - pow(Globals::gamma0,2))/(sqrt(3)*pow(Globals::gamma0,2))) -
          2*log(pow(2, 1.0 / 3.0)*pow(gamma,2) + pow(Globals::gamma0,2)) +
          log(pow(2, 2.0 / 3.0)*pow(gamma,4) - pow(2, 1.0 / 3.0)*pow(gamma,2)*pow(Globals::gamma0,2) +
          pow(Globals::gamma0,4))) - 2*pow(Globals::gamma0,4)*pow(cA,2)*
          (log(2*pow(gamma,6) + pow(Globals::gamma0,6)) - 3*log(fabs((1 - gamma*sqrt(cA))*(1 + gamma*sqrt(cA))))))/
          (2.*pow(2, 1.0 / 6.0)*constants::PI*pow(Globals::gamma0,3)*(2 + pow(Globals::gamma0,6)*pow(cA,3)));
}
double Lambda (double R) {
  double vx = vUdr(R) [0];
  double vy = vUdr(R) [1];
  double sinth = sin(theta_kb(R));
  double costh = cos(theta_kb(R));
  double avrg = INTEGRAL(1000000.0, R) - INTEGRAL(0.0, R);
  return (-1.0 / 2.0) * pow(omegaP(R) * gammaU(R) / omegaW(R), 2) * avrg * (pow(sinth - vx, 2) + pow(vy * costh, 2));
}

double dtau (double R) {
  return pow(omegaP(R), 2) * fDist (fabs(omegaB(R)) / (omegaW(R) * gammaU(R)));
}

/* NEW STUFF HERE /> */

void CreateMas(int N){
  r_perp = new double [N];
  Rtr = new double [N];
  for (int i = 0; i < N; i++) {
    r_perp[i] = 0.0;
    Rtr[i] = 0.0;
  }
}

void DelMas(){
  delete[]r_perp;
  delete[]Rtr;
}

double r_perpFromR (double R1, double R2){
	string path = "my_output/";
	string name0 = "rperp";
	ofstream output0(path + name0 + ".dat");

	double rth = 0.0;

	CreateMas( (int)(1.5*Globals::RESCAPE/10.0) + 1);
  vector <double> r(3);
  vector <double> m(3);
  int i = 0;
  m = vMoment(R2); //take fix moment for a phase
  while(R1 <= R2){ //take point
  	r = vR(R2);
  	rth = sin(psi_m(R2))/sqrt(NORM(r));
  	while(NORM(r) > 1.0){ //integrate along B
  		r = SUM(r, TIMES(-0.01, Bxyz(r, m)));
  		//cout << r[0] << " " << r[1] << " " << r[2] << " " << NORM(r) <<"\n" << endl;
  	}
  	r_perp[i] = ANGLE(r, m);
  	Rtr[i] = R2;
  	output0 << R2 << " " << ANGLE(r, m) << " " << rth << "\n";
  	R2 = R2 - 100.0;
  	i+=1;
  }
  output0.close();
  return 0;
}

vector <double> Bxyz(vector <double> r, vector <double> m){
	vector <double> n(3);
	n = NORMALIZE(r);
  return SUM(TIMES(3.0 * SCALAR(m, n), n), TIMES(-1.0, m));
}

void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
// Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns a value y, and an error estimate dy. If P(x) is the polynomial of degree N − 1 such that P(xai) = yai,i = 1,...,n, then the returned value y = P(x).
{
	int i,m,ns=1;
	double den, dif, dift, ho, hp, w;
	double *c,*d;
	dif = fabs(x-xa[1]);
	c = new double [n];
	d = new double [n];
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift; }
		c[i]=ya[i];
    	d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0){
            	cout << "Error in routine polint\n";
            	break;
       		}
			den=w/den;
			d[i]=hp*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	delete [] d;
	delete [] c;
}
