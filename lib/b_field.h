#pragma B_FIELD

void ReadFile();

//std::vector <double> vB_XYZ(std::vector <double> r, std::vector <double> m);
//std::vector <double> vb_XYZ(std::vector <double> r, double R);

//std::vector <double> vB (double R);
std::vector <double> vb (double R);

std::vector <double> vB_XYZ(std::vector <double> r, double R);
std::vector <double> vB(double R);

std::vector <double> vb_XYZ(std::vector <double> r, double R);
std::vector <double> vb (double R);
void bnorm();
//double norm;

std::vector <double> vB_dipole(std::vector <double> r, std::vector <double> m);
std::vector <double> vb_dipole(std::vector <double> r, std::vector <double> m);

void PrintField();