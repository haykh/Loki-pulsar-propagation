#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include "pc_dens.h"
#include "aux/read_write.h"
#include "aux/functions.h"
#include "constants.h"
#include "initialize.h"

double* r_perp;
double* Rtr;
double* der;

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
  		r = SUM(r, TIMES(-0.01, vB_XYZ(r, m)));
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

double interpolate (double xData[], double yData[], int size, double x, bool extrapolate) {
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
