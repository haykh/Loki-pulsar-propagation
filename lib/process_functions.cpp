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
double* der;


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


string path = "my_output/";
string name0 = "gfunc";
ofstream output0(path + name0 + ".dat");
double gFunc (double R) {
	//int k = (int)(4.0*Globals::RESCAPE/100.0) - 1;
	//double res = 4.0*Globals::RESCAPE - k*100.0;
	//int n = k + (int)(res/10.0)+1;
	int n = 61;
  double rperp;
  //int n = (int)(4.0*Globals::RESCAPE/100.0)+10;
  splint(Rtr, r_perp, der, n, R, &rperp);
  //cout << "ok\n";
	double f = rperp*rperp *Globals::RLC;
	double theta = ANGLE(vR(R), Globals::vOmega);
	double dtheta = 5.0 * constants::PI / 180.0;
	double gap = 1.0;
	if (Globals::alpha_deg > 80)
    	gap = (1. - exp(-pow(constants::PI * 0.5 - theta, 2) / (2.0 * dtheta * dtheta)));
    double fd = pow(sin(psi_m(R)), 2) * Globals::RLC / NORM(vR(R)); 
	//output0 << R << " " << rperp << " " << sin(psi_m(R))/sqrt(NORM(vR(R))) << "\n";
	output0 << R << " " << (pow(fd, 2.5) * exp(-fd * fd) / (pow(fd, 2.5) + pow(Globals::f0, 2.5))) * gap << " " << (pow(f, 2.5) * exp(-f * f) / (pow(f, 2.5) + pow(Globals::f0, 2.5))) * gap <<"\n";
    
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
  der = new double [N];
  for (int i = 0; i < N; i++) {
    r_perp[i] = 0.0;
    Rtr[i] = 0.0;
    der[i] = 0.0;
  }
}

void DelMas(){
  delete[]r_perp;
  delete[]Rtr;
  delete[]der;
}

double r_perpFromR (double R1, double R2){
	//string path = "my_output/";
	string name1 = "rperp";
	ofstream output1(path + name1 + ".dat");
	double rth = 0.0;
	
	//int k = (int)(4.0*Globals::RESCAPE/100.0) - 1;
	//double res = 4.0*Globals::RESCAPE - k*100.0;
	//int n = k + (int)(res/10.0)+1;
	//int n = (int)(4.0*Globals::RESCAPE/100.0) - 1 + (int)(4.0*Globals::RESCAPE - ((int)(4.0*Globals::RESCAPE/100.0) - 1)*100.0/10.0);
	//int n = (int)(4.0*Globals::RESCAPE/100.0)+10;
	//R2 = 4.0*Globals::RESCAPE;
	R2 = 5100.0;
	int n = 61;
	CreateMas(n);
	//cout << "ok \n";
  vector <double> r(3);
  vector <double> m(3);
  int i = n-1;
  while(R1 <= R2){ //take point
  	r = vR(R2);
  	m = vMoment(R2); //take fix moment for a phase
  	rth = sin(psi_m(R2))/sqrt(NORM(r));
  	while(NORM(r) > 1.0){ //integrate along B
  		r = SUM(r, TIMES(-0.01, Bxyz(r, m)));
  		//cout << r[0] << " " << r[1] << " " << r[2] << " " << NORM(r) <<"\n" << endl;
  	}
  	r_perp[i] = ANGLE(r, m);
  	Rtr[i] = R2;
  	//R2-=100.0;
  	//cout << Rtr[i] << " " << rth << " " << r_perp[i] << "\n";
  	//output1 << Rtr[i] << " " << rth << " " << r_perp[i] << "\n";
  	if(R2 >100.0)
  		R2 -= 100.0;
  	else
  		R2 -= 10.0;
  	i-=1;
  }
  spline(Rtr, r_perp, n, (r_perp[1]-r_perp[0])/(Rtr[1]-Rtr[0]), (r_perp[n-1]-r_perp[n-2])/(Rtr[n-1]-Rtr[n-2]), der);
  R2 = 5100.0;
  double rperp = 0.0;
  
  while(R2 > 0){
  	splint(Rtr, r_perp, der, n, R2, &rperp);
  	output1 << R2 << " " << rperp << " " << sin(psi_m(R2))/sqrt(NORM(vR(R2))) << "\n";
  	R2 -= 50.0;
  }
  output1.close();
  return 0;
}

vector <double> vBdipoleXYZ(vector <double> r, vector <double> m){
	vector <double> n(3);
	n = NORMALIZE(r);
  return SUM(TIMES(3.0 * SCALAR(m, n), n), TIMES(-1.0, m));
}

vector <double> vBsplitXYZ (vector <double> r) {
  double Rr = NORM(r);
  double Rxy = sqrt(r[0] * r[0] + r[1] * r[1]);

  double costh = SCALAR(NORMALIZE(r), NORMALIZE(Globals::vOmega));
  double sinth = sqrt(1 - costh * costh);
  double cosphi = r[0] / Rxy;
  double sinphi = r[1] / Rxy;
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

vector <double> Bxyz(vector <double> r, vector <double> m){
	
	if (Globals::fphi == 0 && Globals::fr == 0) {
    	return vBdipoleXYZ(r,m);
	} else {
		return SUM(vBsplitXYZ(r), vBdipoleXYZ(r,m));
	}
}


void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
//Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., yi = f(xi), with x1 < x2 < . . . < xN ,
//and given values yp1 and ypn for the first derivative of the interpolating function at points 1 and n, 
//respectively, this routine returns an array y2[1..n] that contains the second derivatives of the interpolating 
//function at the tabulated points xi. If yp1 and/or ypn are equal to 1 × 1030 or larger, the routine is signaled
//to set the corresponding boundary condition for a natural spline, with zero second derivative on that boundary.
{
	int i,k;
	double p,qn,sig,un,*u;
	u = new double [n-1]; 
	if (yp1 > 0.99e30)
        y2[0]=u[0]=0.0;
    else {
		y2[0] = -0.5;
		//The lower boundary condition is set either to be “nat- ural”
		//or else to have a specified first derivative.
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1); 
	}
	for (i=1;i<n-1;i++) { //This is the decomposition loop of the tridiagonal al-
	//gorithm. y2 and u are used for temporary storage of the decomposed factors.
 		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]); 
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
    	qn=un=0.0;
	else {
		qn=0.5;
		//The upper boundary condition is set either to be “natural”
		//or else to have a specified first derivative.
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[1-2]-y[n-2])/(x[n-1]-x[n-2])); 
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--) //This is the backsubstitution loop of the tridiagonal
		y2[k]=y2[k]*y2[k+1]+u[k]; //algorithm. 
	delete [] u;
}

void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
//Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai’s in order), and given the 
//array y2a[1..n], which is the output from spline above, and given a value of x, this routine returns a 
//cubic-spline interpolated value y.
{
	//void nrerror(char error_text[]); 
	int klo,khi,k;
	double h,b,a;
	//We will find the right place in the table by means of bisection. This is optimal if sequential calls to 
	//this routine are at random values of x. If sequential calls are in order, and closely spaced, one would do 
	//better to store previous values of klo and khi and test if they remain appropriate on the next call.
	klo=0;
	khi=n-1;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) 
			khi=k; 
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) 
		cout << "error: Bad xa input to routine splint" << "\n";
		//nrerror("Bad xa input to routine splint"); 
		//The xa’s must be distinct.
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h; 
	//Cubic spline polynomial is now evaluated. 
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}