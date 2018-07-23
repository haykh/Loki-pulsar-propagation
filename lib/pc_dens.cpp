#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>

#include "pc_dens.h"
#include "aux/read_write.h"
#include "aux/functions.h"
#include "aux/NRutil.h"
#include "b_field.h"
#include "process_functions.h"
#include "constants.h"
#include "initialize.h"

#define SQR(x) ( (x)*(x) )

#ifdef INTBACK

  double* pcdens::rps;
  double* pcdens::Rs;
  int pcdens::N;

  /*std::vector <double> euler_3d (std::vector<double> vr0,
                                std::vector <double> (*vDeriv)(std::vector <double>, std::vector <double>),
                                std::vector <double> arg2,
                                double stepsize) {
    double h = stepsize;

    std::vector <double> vXYZ = vr0;
    while (NORM(vXYZ) > 1.0)
      vXYZ = SUM(vXYZ, TIMES(-h, vDeriv(vXYZ, arg2)));
    return vXYZ;
  }*/

  std::vector <double> euler_3d (std::vector<double> vr0,
                                std::vector <double> (*vDeriv)(std::vector <double>, double),
                                double arg2,
                                double stepsize) {
    double h = stepsize;

    std::vector <double> vXYZ = vr0;
    while (NORM(vXYZ) > 1.0)
      vXYZ = SUM(vXYZ, TIMES(-h, vDeriv(vXYZ, arg2)));
    return vXYZ;
  }

  bool check_rpFromR (string &filename) {
    stringstream filename_stream;
    // filename is: <fr>_<fphi>_<phi0>_<Rem>_<alpha>_<beta>_<P>.rpFromR
    //    > where <P> is period in milliseconds
    //    > <alpha> and <beta> are in degrees * 100
    //    > <phi0> is in degrees * 10
    //    > <fr> & <fphi> are * 10
    //    > `m` stands for minus
    //    example: 7_0_m5_3000_4500_m100_500.rpFromR
    if (Globals::fr < 0)
      filename_stream << "m";
    filename_stream << int(fabs(Globals::fr) * 10) << "_";
    if (Globals::fphi < 0)
      filename_stream << "m";
    filename_stream << int(fabs(Globals::fphi) * 10) << "_";
    if (Globals::PHI0 < 0)
      filename_stream << "m";
    filename_stream << int(fabs(Globals::PHI0) * 10 * 180.0 / constants::PI) << "_";
    filename_stream << int(Globals::R_em) << "_";
    if (Globals::alpha_deg < 0)
      filename_stream << "m";
    filename_stream << int(fabs(Globals::alpha_deg) * 100) << "_";
    if (Globals::beta_deg < 0)
      filename_stream << "m";
    filename_stream << int(fabs(Globals::beta_deg) * 100) << "_";
    filename_stream << int(Globals::Period * 1000) << ".rpFromR";
    filename = "rpFromR_data/" + filename_stream.str();

    struct stat buffer;
    return (stat (filename.c_str(), &buffer) == 0);
  }

  void fill_rpFromR (string filename, double Rmax) {
  	double rup = Rmax;
    std::vector <double> vm_, rups_;
    double euler_step = 0.01 * (Globals::RLC / 5000.0);

    while (rup >= 0.0) {
      rups_.push_back(rup);
      rup -= (rup > (Globals::RLC / 5.)) ? (Globals::RLC / 50.) : ((rup > (Globals::RLC / 50.)) ? (Globals::RLC / 500.) : (Globals::RLC / 1000.));
    }
    rups_.push_back(0.0);

    pcdens::N = rups_.size();
    pcdens::rps = new double [pcdens::N];
    pcdens::Rs = new double [pcdens::N];

    for (int i = 0; i < pcdens::N; i++) {
      rup = rups_[pcdens::N - i - 1];
      #ifndef MPI
        print_progress(rup / Rmax, 26, "\t");
      #endif
    	vm_ = vMoment(rup);
      pcdens::Rs[i] = rup;
      //pcdens::rps[i] = sin(ANGLE(euler_3d (vR(rup), vb_XYZ, vm_, euler_step), vm_));
      // CHANGE THIS:
      pcdens::rps[i] = sin(ANGLE(euler_3d (vR(rup), vb_XYZ, rup, euler_step), vm_));
      cout << rup << " " << pcdens::rps[i] << endl;
    }
    user_cout("");

    struct stat st = {0};
    std::string foldername = filename.substr (0, filename.find("/"));
    if (stat(foldername.c_str(), &st) == -1) {
      mkdir(foldername.c_str(), 0700);
    }
    ofstream output0(filename);
    for (int i = 0; i < pcdens::N; i++) {
      output0 << pcdens::Rs[i] << " " << pcdens::rps[i] << "\n";
    }
    output0.close();
  }

  void read_rpFromR (string filename) {
    std::vector <double> Rs_, rps_;
    ifstream data_file(filename);
    if (data_file.is_open()) {
      while (!data_file.eof()) {
        double x_;
        data_file >> x_;
        Rs_.push_back(x_);
        data_file >> x_;
        rps_.push_back(x_);
      }
    } else {
      throw_error("ERROR: Unable to open `rperFromR` file.");
    }
    pcdens::N = Rs_.size() - 1;
    pcdens::rps = new double [pcdens::N];
    pcdens::Rs = new double [pcdens::N];
    for (int i = 0; i < pcdens::N; i++) {
      pcdens::Rs[i] = Rs_[i];
      pcdens::rps[i] = rps_[i];
    }
  }

  void rpFromR (double Rmax) {
    std::string fname;
    if (check_rpFromR(fname)) {
      user_cout("\t`rpFromR` file found. Reading...");
      read_rpFromR (fname);
    } else {
      user_cout("\tNo `rpFromR` file found. Computing and writing...");
      fill_rpFromR (fname, Rmax);
    }
  }

#endif

double gFunc (double R) {
  #ifdef INTBACK
    double rp = interpolate(pcdens::Rs, pcdens::rps, pcdens::N, R, false);
  #else
    double rp = sin(psi_m(R)) / sqrt(NORM(vR(R)));
  #endif
  double f = SQR(rp) * Globals::RLC;
	double theta = ANGLE(vR(R), Globals::vOmega);
	double dtheta = 5.0 * constants::PI / 180.0;
	double gap = 1.0;
	if (Globals::alpha_deg > 80)
    gap = (1. - exp(-pow(constants::PI * 0.5 - theta, 2) / (2.0 * dtheta * dtheta)));
	return (pow(f, 2.5) * exp(-f * f) / (pow(f, 2.5) + pow(Globals::f0, 2.5))) * gap;
}
