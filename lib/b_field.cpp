#include <math.h>
#include <vector>

#include "initialize.h"
#include "constants.h"
#include "aux/functions.h"
#include "process_functions.h"


vector <double> vBdipole_XYZ(vector <double> r, vector <double> m){
	vector <double> n(3);
	n = NORMALIZE(r);
  return SUM(TIMES(3.0 * SCALAR(m, n), n), TIMES(-1.0, m));
}

vector <double> vBsplit_XYZ (vector <double> r) {
  double Rr = NORM(r);
  double Rxy = sqrt(r[0] * r[0] + r[1] * r[1]);

  double costh = SCALAR(NORMALIZE(r), NORMALIZE(Globals::vOmega));
  double sinth = sqrt(1 - costh * costh);
  double cosphi = r[0] / Rxy;
  double sinphi = r[1] / Rxy;
  double phi = acos (cosphi);

  // double psi1 = costh * cos(Globals::alpha) + sinth * sin(Globals::alpha) * cos(phi - Globals::PHI0 + Rr / Globals::RLC);

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

vector <double> vB_XYZ(vector <double> r, vector <double> m) {
	if (Globals::fphi == 0 && Globals::fr == 0) {
    	return vBdipole_XYZ(r,m);
	} else {
		return SUM(vBsplit_XYZ(r), vBdipole_XYZ(r,m));
	}
}

vector <double> vBdipole (double R) {
  return vBdipole_XYZ(vR(R), vMoment(R));
}

vector <double> vBsplit (double R) {
  return vBsplit_XYZ(vR(R));
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
