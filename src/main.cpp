#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <typeinfo>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>

using namespace std;

#include "../lib/constants.h"
#include "../lib/pulsar_parameters.h"

#include "../lib/functions.h"
#include "../lib/process_functions.h"
#include "../lib/auxiliary.h"

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
  // Reading flags >
  string input_name, out_path = "output";
  bool found_input = false;
  for (int i=0; i<argc; ++i) {
		if (typeid(argv[i]) == typeid(char*)) {
			if (strncmp(argv[i], "-i", 3) == 0) {
				if (i+1 < argc) {
					input_name = argv[i+1];
          found_input = true;
        }
				else {
          throw_error("ERROR: No input file given.");
				}
			}
      if (strncmp(argv[i], "-o", 3) == 0) {
				if (i+1 < argc) {
					out_path = argv[i+1];
        }
				else {
          throw_error("ERROR: No output path given.");
				}
			}
		}
	}
  if (!found_input) {
    throw_error("ERROR: No input file given.");
  } else {
    cout << "INPUT: " << input_name << "\n";
  }
  // < Reading flags

  define_Globals(input_name);

  string run_id = read_from_file_str(input_name, "run_id", "my_run");
  cout << "RUN_ID: " << run_id << "\n\n";

  cout << "\nR_esc: " << Globals::RESCAPE << endl;
  cout << "R_A: " << Globals::ROMODE << endl;
  cout << "R_lc: " << Globals::RLC << endl << endl;

  string MODE;
  if (Globals::mode == 0) MODE = "X-mode";
  else MODE = "O-mode";

  struct stat st = {0};
  if (stat(out_path.c_str(), &st) == -1) {
      mkdir(out_path.c_str(), 0700);
  }

  ofstream outputData(out_path + "/" + run_id + ".dat");
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

//
//    SIMULATION STARTS HERE
//

  ofstream output0(out_path + "/" + run_id + "[0].dat");
  ofstream output1(out_path + "/" + run_id + "[1].dat");

  /*ofstream plot(out_path + "betadelta.dat");
  PHI0 = 5.0 * constants::PI / 180.0;
  findInitPoints(PHI0);
  for (double r = 58.0; r <= 62.0; r += 0.1) {
  	plot << r << " " << BetaB(r) + delta(r) << "\n";
  	//plot << r << " " << cos(theta_kb(r)) << "\n";
  }
  plot.close();
  return 0;*/

  double phi_t_start = read_from_file(input_name, "phi_start");
  double phi_t_end = read_from_file(input_name, "phi_end");
  double phi_t_step = read_from_file(input_name, "phi_step");
  for (double phi_t = phi_t_start; phi_t <= phi_t_end; phi_t += phi_t_step) { // Phase switch
    cout << "PHI: " << phi_t << endl;
    Globals::PHI0 = phi_t * constants::PI / 180.0;
    findInitPoints (Globals::PHI0);

    double x1, x2, dep_vars[2];
    x1 = 0.0;
    // if (RESCAPE < 1000.0) x2 = 1500.0;
    x2 = 1.5 * Globals::RESCAPE;

    //freeM();
    r_perpFromR (x1, x2);
   // PlotRperp ();

    /*string name1 = "theory_r_perp";
	   ofstream output2(out_path + name1 + ".dat");
	    double RR = x2;
	     while(x1 < RR){
		     output2 << RR << " " << psi_m(RR)/sqrt(NORM(vR(RR)))<< endl;
		       RR = RR - 100.0;
	        }
	         output2.close();
           //cout << "r_perp for initial point from dipol: " << psi_m(x2)/sqrt(NORM(vR(x2))) << "\n" << endl;*/

    /*Initial values*/
    if (Globals::mode == 0) { // X-mode
      dep_vars[0] = BetaB(x1) + delta(x1) + constants::PI / 2.0;
      dep_vars[1] = Arcsinh(1.0 / Q(x1)) / 2.0;
    } else { // O-mode
      dep_vars[0] = BetaB(x1) + delta(x1);
      dep_vars[1] = Arcsinh(-1.0 / Q(x1)) / 2.0;
    }
    /*--------------*/

    double PA = dep_vars[0] * 180 / constants::PI;
    double tau = constants::PI * constants::R_star * integrate(dtau, x1, Globals::RLC) / (constants::c * Globals::omega);
    double gf = gFunc(x1);
    double II0 = gf;
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

    DelMas();
  }
  output0.close();
  output1.close();
	return 0;
}
