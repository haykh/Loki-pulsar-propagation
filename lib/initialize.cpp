#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "aux/read_write.h"
#include "aux/functions.h"
#include "constants.h"
#include "initialize.h"
using namespace std;

// finding initial point of generation for dipole field b0 (along field line)

vector <double> b0 (double th, double ph, double PHI0) {
  vector <double> n0(3);
  n0[0] = sin(th) * cos(ph);
  n0[1] = sin(th) * sin(ph);
  n0[2] = cos(th);

  vector <double> m(3);
  m[0] = sin(Globals::alpha) * cos(Globals::PHI0);
  m[1] = sin(Globals::alpha) * sin(Globals::PHI0);
  m[2] = cos(Globals::alpha);

  return NORMALIZE(SUM(TIMES(3.0 * SCALAR(m, n0), n0), TIMES(-1.0, m)));
}
double func1 (double x, double y, double PHI0) {
  vector <double> o(3);
  o[0] = sin(Globals::dzeta);
  o[1] = 0.0;
  o[2] = cos(Globals::dzeta);
  return CROSS(b0(x, y, PHI0), o)[0];
}
double func2 (double x, double y, double PHI0) {
  vector <double> o(3);
  o[0] = sin(Globals::dzeta);
  o[1] = 0.0;
  o[2] = cos(Globals::dzeta);
  return CROSS(b0(x, y, PHI0), o)[1];
}
double DX (double (*func)(double, double, double), double x, double y, double PHI0) {
  double h = 0.00001;
  double fm2 = func(x - 2 * h, y, PHI0);
  double fp2 = func(x + 2 * h, y, PHI0);
  double fm1 = func(x - h, y, PHI0);
  double fp1 = func(x + h, y, PHI0);
  return (fm2 - 8 * fm1 + 8 * fp1 - fp2) / (12 * h);
}
double DY (double (*func)(double, double, double), double x, double y, double PHI0) {
  double h = 0.00001;
  double fm2 = func(x, y - 2 * h, PHI0);
  double fp2 = func(x, y + 2 * h, PHI0);
  double fm1 = func(x, y - h, PHI0);
  double fp1 = func(x, y + h, PHI0);
  return (fm2 - 8 * fm1 + 8 * fp1 - fp2) / (12 * h);
}

void findInitPoints (double PHI0) {
  double X = Globals::alpha, Y = Globals::PHI0;
  for(int i = 0; i < 100; i ++) {
      double f1x = DX(func1, X, Y, Globals::PHI0);
      double f2x = DX(func2, X, Y, Globals::PHI0);
      double f1y = DY(func1, X, Y, Globals::PHI0);
      double f2y = DY(func2, X, Y, Globals::PHI0);
      double f1 = func1(X, Y, Globals::PHI0);
      double f2 = func2(X, Y, Globals::PHI0);
      double dX = (f1y * f2 - f1 * f2y) / (f1x * f2y - f1y * f2x);
      double dY = (f1x * f2 - f1 * f2x) / (f1y * f2x - f2y * f1x);
      X += dX;
      Y += dY;
  }
  Globals::theta_em = X;
  Globals::phi_em = Y;
}

// initialize Globals >
double Globals::theta_em, Globals::phi_em,
      Globals::B12, Globals::B0,
      Globals::freqGHz, Globals::omega,
      Globals::Period, Globals::Omega,
      Globals::lambda, Globals::gamma0, Globals::f0,
      Globals::mode, Globals::fr, Globals::fphi,
      Globals::alpha_deg, Globals::beta_deg, Globals::alpha, Globals::beta, Globals::dzeta,
      Globals::PHI0,
      Globals::R_em, Globals::RLC, Globals::RESCAPE, Globals::ROMODE;
vector <double> Globals::vOmega;
string Globals::RUN_ID, Globals::input_name, Globals::out_path;
// < initialize Globals

void define_Globals() {
  Globals::B12 = read_from_file(Globals::input_name, "B12"); // Surface B-field in 10^12 Gs
  Globals::Period = read_from_file(Globals::input_name, "Period"); // Rotation period in sec
  Globals::freqGHz = read_from_file(Globals::input_name, "freqGHz"); // Radiation frequency in GHz

  Globals::lambda = read_from_file(Globals::input_name, "lambda"); // Plasma multiplicity // normally 10000
  Globals::gamma0 = read_from_file(Globals::input_name, "gamma0"); // Plasma mean gamma factor // normally 50
  Globals::f0 = read_from_file(Globals::input_name, "f0"); // Polar Cap gap width // normally 0.5
  Globals::R_em = read_from_file(Globals::input_name, "R_em"); // Emission radius in star radii

  Globals::mode = read_from_file(Globals::input_name, "mode"); // 1 = O-mode & 0 = X-mode
  Globals::fr = read_from_file(Globals::input_name, "fr"); // Split-monopole parameter
  Globals::fphi = read_from_file(Globals::input_name, "fphi"); // Split-monopole parameter

  Globals::alpha_deg = read_from_file(Globals::input_name, "alpha_deg");
  Globals::beta_deg = read_from_file(Globals::input_name, "beta_deg");

  Globals::alpha = Globals::alpha_deg * constants::PI / 180.0; // Inclination angle in radians
  Globals::beta = Globals::beta_deg * constants::PI / 180.0; // Line of sight angle in radians;
  Globals::dzeta = Globals::alpha - Globals::beta; // Minimum angle between the rotation axis and the line of sight

  Globals::B0 = Globals::B12 * 1.0e12; // Surface magnetic field in Gs
  Globals::Omega = 2.0 * constants::PI / Globals::Period; // Rotation frequency
  Globals::omega = 2.0 * constants::PI * Globals::freqGHz * 1.0e9; // Radiation circular frequency

  Globals::vOmega.push_back(0.);
  Globals::vOmega.push_back(0.);
  Globals::vOmega.push_back(Globals::Omega);

  Globals::RLC = (constants::c / Globals::Omega) / constants::R_star;
  Globals::RESCAPE = 1.0e3 * pow(Globals::lambda / 1.0e4, 1.0/3.0) * pow(Globals::gamma0 / 100.0, -6.0/5.0) * pow(Globals::B0 / 1.0e12, 2.0/5.0) * pow(Globals::freqGHz, -2.0/5.0) * pow(Globals::Period, -1.0/5.0);
  Globals::ROMODE = 1.0e2 * pow(Globals::lambda / 1.0e4, 1.0/3.0) * pow(Globals::gamma0 / 100.0, 1.0/3.0) * pow(Globals::B0 / 1.0e12, 1.0/3.0) * pow(Globals::freqGHz, -2.0/3.0) * pow(Globals::Period, -1.0/3.0);
}

void initialize(int argc, char* argv[]) {
  read_in_out(Globals::input_name, Globals::out_path, argc, argv);

  define_Globals();

  Globals::RUN_ID = read_from_file_str(Globals::input_name, "run_id", "my_run");
  cout << "RUN_ID: " << Globals::RUN_ID << "\n\n";

  cout << "\nR_esc: " << Globals::RESCAPE << endl;
  cout << "R_A: " << Globals::ROMODE << endl;
  cout << "R_lc: " << Globals::RLC << endl << endl;

  string MODE;
  if (Globals::mode == 0) MODE = "X-mode";
  else MODE = "O-mode";

  // Create output directory if doesn't exist
  struct stat st = {0};
  if (stat(Globals::out_path.c_str(), &st) == -1) {
    mkdir(Globals::out_path.c_str(), 0700);
  }

  ofstream outputData(Globals::out_path + "/" + Globals::RUN_ID + ".dat");
  outputData
      << "alpha = " << Globals::alpha_deg
      << "\nbeta = " << Globals::beta_deg
      << "\n\nPeriod = " << Globals::Period
      << "\nB12 = " << Globals::B12
      << "\nfGHz = " << Globals::freqGHz
      << "\n\nlambda = " << Globals::lambda
      << "\ngamma0 = " << Globals::gamma0
      << "\nf0 = " << Globals::f0
      << "\nR_em = " << Globals::R_em
      << "\n\n" + MODE
      << "\nfr = " << Globals::fr
      << "\nfphi = " << Globals::fphi
      << "\n\n\nR_LC = " << Globals::RLC
      << "\nR_escape = " << Globals::RESCAPE
      << "\nR_A = " << Globals::ROMODE;
  outputData.close();
}
