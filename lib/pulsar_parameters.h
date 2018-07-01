#pragma PULSAR_PARAMETERS
using namespace std;

void throw_error(const string msg);
double read_from_file (string input_name, const string param);
string read_from_file_str (string input_name, const string param);
void define_Globals(string input_name);

namespace Globals {
  extern double B12;
  extern double B0;

  extern double freqGHz;
  extern double omega;

  extern double Period;
  extern double Omega;

  extern double lambda;
  extern double gamma0;
  extern double f0;

  extern double mode;
  extern double fr;
  extern double fphi;
  extern double BMULT;

  extern double alpha_deg;
  extern double beta_deg;

  extern double alpha;
  extern double beta;
  extern double dzeta;

  extern double PHI0;

  extern vector <double> vOmega;

  extern double R_em;

  extern double RLC;
  extern double RESCAPE;
  extern double ROMODE;
}
