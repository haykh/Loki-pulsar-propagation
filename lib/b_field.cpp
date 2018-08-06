#include <math.h>
#include <vector>

#include "initialize.h"
#include "constants.h"
#include "aux/functions.h"
#include "process_functions.h"
#include "aux/NRutil.h"


std::vector <double> vBdipole_XYZ(std::vector <double> r, std::vector <double> m){
	std::vector <double> n(3);
	n = NORMALIZE(r);
  return SUM(TIMES(3.0 * SCALAR(m, n), n), TIMES(-1.0, m));
}

std::vector <double> vBsplit_XYZ (std::vector <double> r) {
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

  std::vector <double> temp(3);
  temp[0] = Br * sinth * cosphi - Bphi * sinphi;
  temp[1] = Br * sinth * sinphi + Bphi * cosphi;
  temp[2] = Br * costh;
  return temp;
}

std::vector <double> vB_XYZ(std::vector <double> r, std::vector <double> m) {
	if (Globals::fphi == 0 && Globals::fr == 0) {
    return vBdipole_XYZ(r, m);
	} else {
		return SUM(vBsplit_XYZ(r), vBdipole_XYZ(r,m));
	}
}

std::vector <double> vb_XYZ(std::vector <double> r, std::vector <double> m) {
  return NORMALIZE(vB_XYZ(r, m));
}

std::vector <double> vBdipole (double R) {
  return vBdipole_XYZ(vR(R), vMoment(R));
}

std::vector <double> vBsplit (double R) {
  return vBsplit_XYZ(vR(R));
}

std::vector <double> vB (double R) {
  return vB_XYZ(vR(R), vMoment(R));
}
std::vector <double> vb (double R) {
  return NORMALIZE(vB(R));
}
