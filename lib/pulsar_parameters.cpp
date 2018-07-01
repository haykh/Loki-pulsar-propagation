#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>

#include <sstream>
#include <algorithm>
#include <iterator>
using namespace std;

#include "constants.h"
#include "pulsar_parameters.h"

void throw_error(const string msg) {
  cout << msg << endl;
  exit (EXIT_FAILURE);
}
bool is_number(const string& s) {
    if (s.empty())
        return false;
    bool trigger = false;
    if (s[0] == '.') {
        return false;
    }
    for (int i = 0; i < s.size(); i ++) {
        char c = s[i];
        if (c == '.' && !trigger) {
            trigger = !trigger;
            continue;
        } else if (c == '.' && trigger) {
            return false;
        }
        if (i == 0 && c == '-') {
            continue;
        }
        if (!std::isdigit(c))
          return false;
    }
    return true;
}
double read_from_file (string input_name, const string param) {
    ifstream infile(input_name);
    string str;
    stringstream msg;
    if (infile.is_open()) {
        while (getline(infile, str)) {
          // cout << iter << ": " << str << "\n";
          istringstream iss(str);
          vector<string> words{istream_iterator<string>{iss}, istream_iterator<string>{}};
          if(words.size() < 1) continue;
          // cout << words[0] << "\n";
          if (words[0] == param) {
              if(is_number(words[1])) {
                  return atof(words[1].c_str());
              } else {
                  msg << "ERROR: Converting string to number for parameter '" << param << "' in the input file.";
                  throw_error(msg.str());
              }
              break;
          }
        }
    } else {
        throw_error("ERROR: Cannot open the input file.");
    }
    msg << "ERROR: Cannot find the given parameter '" << param << "' in the input file.";
    throw_error(msg.str());
    return 0;
}
string read_from_file_str (string input_name, const string param) {
    ifstream infile(input_name);
    string str;
    stringstream msg;
    if (infile.is_open()) {
        while (getline(infile, str)) {
          istringstream iss(str);
          vector<string> words{istream_iterator<string>{iss}, istream_iterator<string>{}};
          if(words.size() < 1) continue;
          if (words[0] == param) {
            if(!is_number(words[1])) {
              return words[1].c_str();
            } else {
              msg << "ERROR: param '" << param << "' is not a string in the input file.";
              throw_error(msg.str());
            }
            break;
          }
        }
    } else {
      throw_error("ERROR: Cannot open the input file.");
    }
    msg << "ERROR: Cannot find the given parameter '" << param << "' in the input file.";
    throw_error(msg.str());
    return "NaN";
}

// initialize Globals >
double Globals::B12;
double Globals::B0;

double Globals::freqGHz;
double Globals::omega;

double Globals::Period;
double Globals::Omega;

double Globals::lambda;
double Globals::gamma0;
double Globals::f0;

double Globals::mode;
double Globals::fr;
double Globals::fphi;
double Globals::BMULT;

double Globals::alpha_deg;
double Globals::beta_deg;

double Globals::alpha;
double Globals::beta;
double Globals::dzeta;

double Globals::PHI0;

vector <double> Globals::vOmega;

double Globals::R_em;

double Globals::RLC;
double Globals::RESCAPE;
double Globals::ROMODE;
// < initialize Globals

void define_Globals(string input_name) {
  Globals::B12 = read_from_file(input_name, "B12"); // Surface B-field in 10^12 Gs
  Globals::Period = read_from_file(input_name, "Period"); // Rotation period in sec
  Globals::freqGHz = read_from_file(input_name, "freqGHz"); // Radiation frequency in GHz

  Globals::lambda = read_from_file(input_name, "lambda"); // Plasma multiplicity // normally 10000
  Globals::gamma0 = read_from_file(input_name, "gamma0"); // Plasma mean gamma factor // normally 50
  Globals::f0 = read_from_file(input_name, "f0"); // Polar Cap gap width // normally 0.5
  Globals::R_em = read_from_file(input_name, "R_em"); // Emission radius in star radii

  Globals::mode = read_from_file(input_name, "mode"); // 1 = O-mode & 0 = X-mode
  Globals::fr = read_from_file(input_name, "fr"); // Split-monopole parameter
  Globals::fphi = read_from_file(input_name, "fphi"); // Split-monopole parameter

  Globals::alpha_deg = read_from_file(input_name, "alpha_deg");
  Globals::beta_deg = read_from_file(input_name, "beta_deg");

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
