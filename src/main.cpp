#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <typeinfo>

#ifdef MPI
  #include "mpi.h"
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>

using namespace std;

#include "../lib/constants.h"
#include "../lib/read_write.h"

#include "../lib/functions.h"
#include "../lib/process_functions.h"
#include "../lib/initialize.h"

#include "../lib/integrator.h"
#include "../lib/RHS.h"

/*
  Default libraries
*/
#include "../lib/diffeqsolver.h"

void displayVector (vector <double> a) {
    cout << endl << a[0] << endl << a[1] << endl << a[2] << endl;
}

int main(int argc, char* argv[]) {
  #ifdef MPI
    cout << "mpi enabled\n";
  #endif
  initialize(argc, argv);

  ofstream output0(Globals::out_path + "/" + Globals::RUN_ID + "_0.dat");
  ofstream output1(Globals::out_path + "/" + Globals::RUN_ID + "_1.dat");

  // SIMULATION STARTS HERE />

  double phi_t_start = read_from_file(Globals::input_name, "phi_start");
  double phi_t_end = read_from_file(Globals::input_name, "phi_end");
  double phi_t_step = read_from_file(Globals::input_name, "phi_step");
  for (double phi_t = phi_t_start; phi_t <= phi_t_end; phi_t += phi_t_step) { // Phase switch
    cout << "PHI: " << phi_t << endl;
    Globals::PHI0 = phi_t * constants::PI / 180.0;
    findInitPoints (Globals::PHI0);

    double x1, x2, dep_vars[2];
    x1 = 0.0;
    x2 = 1.5 * Globals::RESCAPE;

    // Initial values />
    if (Globals::mode == 0) { // X-mode
      dep_vars[0] = BetaB(x1) + delta(x1) + constants::PI / 2.0;
      dep_vars[1] = Arcsinh(1.0 / Q(x1)) / 2.0;
    } else { // O-mode
      dep_vars[0] = BetaB(x1) + delta(x1);
      dep_vars[1] = Arcsinh(-1.0 / Q(x1)) / 2.0;
    }
    // </ Initial values

    double PA = dep_vars[0] * 180 / constants::PI;
    double tau = constants::PI * constants::R_star * integrate(dtau, x1, Globals::RLC) / (constants::c * Globals::omega);
    double II0 = gFunc(x1);
    double II = II0 * exp (-tau);

    double VV = II * tanh(2.0 * dep_vars[1]);
    output0 << phi_t << " " << II0 << " " << VV << " " << PA << endl;

    int nvar = 2, nok = 0, nbad = 0;
    double deps = 1.0, h1 = 1.0e-14, hmin = 1.0e-15;
    odeint(dep_vars, nvar, x1, x2, deps, h1, hmin, nok, nbad, RHS);

    VV = II * tanh(2.0 * dep_vars[1]);
    PA = dep_vars[0] * 180.0 / constants::PI;

    // cout << "\tI: " << II << "\n\tV: " << VV << "\n\tPA: " << -PA << endl << endl;
    output1 << phi_t << " " << II << " " << VV << " " << -PA << endl;
  }
  output0.close();
  output1.close();
	return 0;
}
