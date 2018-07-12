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

using namespace std;

#include "../lib/constants.h"
#include "../lib/aux/read_write.h"
#include "../lib/b_field.h"

#include "../lib/aux/functions.h"
#include "../lib/process_functions.h"
#include "../lib/pc_dens.h"
#include "../lib/initialize.h"

#include "../lib/aux/integrator.h"
#include "../lib/RHS.h"
#include "../lib/aux/diffeqsolver.h"

void displayVector (vector <double> a) {
  cout << endl << a[0] << endl << a[1] << endl << a[2] << endl;
}

int main(int argc, char* argv[]) {
  initialize(argc, argv);

  string fname0_ = Globals::out_path + "/" + Globals::RUN_ID + "_0.out";
  string fname1_ = Globals::out_path + "/" + Globals::RUN_ID + "_1.out";
  #ifndef MPI
    ofstream output0(fname0_);
    ofstream output1(fname1_);
  #else
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &mpi::size);
    MPI_Comm_rank (MPI_COMM_WORLD, &mpi::rank);
  #endif

  double phi_t_start = read_from_file(Globals::input_name, "phi_start");
  double phi_t_end = read_from_file(Globals::input_name, "phi_end");
  double phi_t_step = read_from_file(Globals::input_name, "phi_step");

  #ifdef MPI
    int total_steps = (int)((phi_t_end - phi_t_start) / phi_t_step) + 1;
    vector <vector <double> > steps_rank(mpi::size);
    int rnk = 0;
    for (double phi_t = phi_t_start; phi_t <= phi_t_end; phi_t += phi_t_step) {
      steps_rank[rnk % mpi::size].push_back(phi_t);
      rnk ++;
    }
    if (mpi::rank == 0) mpi::offset = 0;
    else mpi::offset = mpi::rank * (steps_rank[mpi::rank - 1].size() * 4) * sizeof(double);
    double *buffer0 = new double [steps_rank[mpi::rank].size() * 4];
    double *buffer1 = new double [steps_rank[mpi::rank].size() * 4];
    int buff_step = 0;
  #endif

  #ifndef MPI
    for (double phi_t = phi_t_start; phi_t <= phi_t_end; phi_t += phi_t_step) { // Phase switch
      cout << "PHI: " << phi_t << endl;
  #else
    for (auto phi_t : steps_rank[mpi::rank]) {
      cout << "rank: #" << mpi::rank << " | PHI: " << phi_t << endl;
  #endif
      Globals::PHI0 = phi_t * constants::PI / 180.0;
      findInitPoints (Globals::PHI0);

      // Initial values & limits />
      double x1, x2, dep_vars[2];
      x1 = 10.0;
      x2 = 2. * Globals::RESCAPE;

      #ifdef INTBACK
        cout << "Calculating r_perp(R) array...\n";
        rpFromR (Globals::RLC);
        cout << "r_perp(R) array done.\n";
      #endif

      if (Globals::mode == 0) { // X-mode
        dep_vars[0] = BetaB(x1) + delta(x1) + constants::PI / 2.0;
        dep_vars[1] = Arcsinh(1.0 / Q(x1)) / 2.0;
      } else { // O-mode
        dep_vars[0] = BetaB(x1) + delta(x1);
        dep_vars[1] = Arcsinh(-1.0 / Q(x1)) / 2.0;
      }
      // </ Initial values & limits

      double PA = dep_vars[0] * 180 / constants::PI;
      double tau = constants::PI * constants::R_star * integrate(dtau, x1, Globals::RLC) / (constants::c * Globals::omega);
      double II0 = gFunc(x1);
      double II = II0 * exp (-tau);

      double VV = II * tanh(2.0 * dep_vars[1]);
      #ifndef MPI
        output0 << phi_t << " " << II0 << " " << VV << " " << PA << endl;
      #else
        buffer0[buff_step] = phi_t;
        buffer0[buff_step+1] = II0;
        buffer0[buff_step+2] = VV;
        buffer0[buff_step+3] = PA;
      #endif

      #ifndef MPI
        cout << "Solving ODE...\n";
      #endif
      odeint(dep_vars, 2, x1, x2, 1.0, 1e-14, 1e-15, 0, 0, RHS);
      #ifndef MPI
        cout << "ODE done.\n";
      #endif

      VV = II * tanh(2.0 * dep_vars[1]);
      PA = dep_vars[0] * 180.0 / constants::PI;
      #ifndef MPI
        cout << "\tI: " << II << "\n\tV: " << VV << "\n\tPA: " << -PA << endl << endl;
        output1 << phi_t << " " << II << " " << VV << " " << -PA << endl;
      #else
        buffer1[buff_step] = phi_t;
        buffer1[buff_step+1] = II;
        buffer1[buff_step+2] = VV;
        buffer1[buff_step+3] = PA;
        buff_step += 4;
      #endif

      #ifdef INTBACK
        delete[] pcdens::rps;
        delete[] pcdens::Rs;
      #endif
  }
  #ifndef MPI
    output0.close();
    output1.close();
  #else
    MPI_Status status;

    MPI_File output0;
    MPI_File_open(MPI_COMM_WORLD, fname0_.c_str(),
                MPI_MODE_CREATE|MPI_MODE_WRONLY,
                MPI_INFO_NULL, &output0);
    MPI_File_write_at(output0, mpi::offset, buffer0, steps_rank[mpi::rank].size() * 4, MPI_DOUBLE, &status);
  	MPI_File_close (&output0);

    MPI_File output1;
    MPI_File_open(MPI_COMM_WORLD, fname1_.c_str(),
        MPI_MODE_CREATE|MPI_MODE_WRONLY,
        MPI_INFO_NULL, &output1);
    MPI_File_write_at(output1, mpi::offset, buffer1, steps_rank[mpi::rank].size() * 4, MPI_DOUBLE, &status);
  	MPI_File_close (&output1);

  	MPI_Finalize ();
  #endif
	return 0;
}
