#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include "aux/read_write.h"
#include "aux/functions.h"
#include "constants.h"
#include "process_functions.h"
#include "b_field.h"
#include "initialize.h"
#include "pc_dens.h"
using namespace std;

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
