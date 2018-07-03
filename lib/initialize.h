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
  extern vector <double> vOmega;
}

void findInitPoints (double PHI0);

void define_Globals(string input_name);
