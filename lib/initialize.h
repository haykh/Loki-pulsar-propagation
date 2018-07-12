#pragma INITIALIZE

namespace Globals {
  extern double theta_em, phi_em,
        B12, B0,
        freqGHz, omega,
        Period, Omega,
        lambda, gamma0, f0,
        mode, fr, fphi,
        alpha_deg, beta_deg, alpha, beta, dzeta,
        PHI0,
        R_em, RLC, RESCAPE, ROMODE;
  extern std::vector <double> vOmega;
  extern std::string RUN_ID, input_name, out_path;
}

void findInitPoints (double PHI0);

void define_Globals();

void initialize(int argc, char* argv[]);

#ifdef MPI
  namespace mpi {
    extern int rank, size, offset;
  }
#endif
