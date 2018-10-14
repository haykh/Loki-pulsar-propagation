#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>

#include <math.h>

#include "H5Cpp.h"
#include "aux/read_write.h"
#include "b_field.h"

using std::cout;
using std::endl;
using namespace H5;

#include "initialize.h"
#include "constants.h"
#include "aux/functions.h"
#include "process_functions.h"
#include "aux/NRutil.h"

//ofstream plot("b_field.dat");

const H5std_string	FILE_NAME("dipole.hdf5");
const H5std_string	DATASET_NAME0("bx");
const H5std_string	DATASET_NAME1("by");
const H5std_string	DATASET_NAME2("bz");

const int 	DIM0 = 257;	               // dataset dimensions
const int 	DIM1 = 257;
const int 	DIM2 = 257;

double *bx = new double[DIM0*DIM1*DIM2];    // buffer for data to write
double *by = new double[DIM0*DIM1*DIM2];
double *bz = new double[DIM0*DIM1*DIM2];

void ReadFile(){
	int i, j, k;
    for (j = 0; j < DIM0; j++){
		for (i = 0; i < DIM1; i++){
			for (k = 0; k < DIM2; k++){
				bx[(i*DIM0+j)*DIM0+k] = 0.0;
				by[(i*DIM0+j)*DIM0+k] = 0.0;
				bz[(i*DIM0+j)*DIM0+k] = 0.0;
			}
		}
	}
	
	try
    {
	// Turn off the auto-printing when failure occurs so that we can
	// handle the errors appropriately
	//Exception::dontPrint();

	// Open an existing file and dataset.
	H5File file(FILE_NAME, H5F_ACC_RDWR);
	
	//Open dataset0 for bx
	DataSet dataset0 = file.openDataSet(DATASET_NAME0);
	DataSpace dataspace0 = dataset0.getSpace();
	hsize_t dims0[3];
    int rank0 = dataspace0.getSimpleExtentDims(dims0);
    DataSpace memspace0(rank0, dims0);
    dataset0.read(bx, PredType::NATIVE_DOUBLE, memspace0, dataspace0);
    
    //Open dataset0 for by
    DataSet dataset1 = file.openDataSet(DATASET_NAME1);
	DataSpace dataspace1 = dataset1.getSpace();
	hsize_t dims1[3];
    int rank1 = dataspace1.getSimpleExtentDims(dims1);
    DataSpace memspace1(rank1, dims1);
    dataset1.read(by, PredType::NATIVE_DOUBLE, memspace1, dataspace1);
    
    //Open dataset0 for bz
	DataSet dataset2 = file.openDataSet(DATASET_NAME2);
	DataSpace dataspace2 = dataset2.getSpace();
	hsize_t dims2[3];
    int rank2 = dataspace0.getSimpleExtentDims(dims2);
    DataSpace memspace2(rank2, dims2);
    dataset2.read(bz, PredType::NATIVE_DOUBLE, memspace2, dataspace2);

    }  // end of try block

    // catch failure caused by the H5File operations
    catch(FileIException error)
    {
		throw_error("h5 error");
    }

    // catch failure caused by the DataSet operations
    catch(DataSetIException error)
    {
    	throw_error("h5 error");
    }
    //cout << bx[(200*DIM0+200)*DIM0+200] << " " << by[(200*DIM0+200)*DIM0+200] << " " << bz[(200*DIM0+200)*DIM0+200] << "\n";
}

std::vector <double> vB_XYZ(std::vector <double> r, double R) {
	std::vector <double> b(3);
	double def = 1.0;
	Point3	p;
	double R_star_n = 1.0;
    
    /*p.x =  r[0] * cos(Globals::PHI0 + R/Globals::RLC) + r[1] * sin(Globals::PHI0 + R/Globals::RLC);
    p.y = -r[0] * sin(Globals::PHI0 + R/Globals::RLC) + r[1] * cos(Globals::PHI0 + R/Globals::RLC);
    p.x = R_star_n*p.z+(DIM0-1.0)/2.0;
    p.y = R_star_n*p.y+(DIM1-1.0)/2.0;
    p.z = R_star_n*r[2]+(DIM2-1.0)/2.0;*/
    
    Globals::PHI0 = 20.0;
    p.x =  r[0] * cos(Globals::PHI0 + R/Globals::RLC) + r[1] * sin(Globals::PHI0 + R/Globals::RLC);
    p.y = -r[0] * sin(Globals::PHI0 + R/Globals::RLC) + r[1] * cos(Globals::PHI0 + R/Globals::RLC);
    p.x = R_star_n*p.x+(DIM0-1.0)/2.0;
    p.y = R_star_n*p.y+(DIM1-1.0)/2.0;
    p.z = R_star_n*r[2]+(DIM2-1.0)/2.0;

    /*Point3 rad;
    rad.x = R_star_n*sin(Globals::alpha);
    rad.y = 0.0;
    rad.z = R_star_n*cos(Globals::alpha);
    rad.x =  rad.x * cos(Globals::PHI0 + R/Globals::RLC) + rad.y * sin(Globals::PHI0 + R/Globals::RLC);
    rad.y = -rad.x * sin(Globals::PHI0 + R/Globals::RLC) + rad.y * cos(Globals::PHI0 + R/Globals::RLC);
    
    rad.x += (DIM0-1.0)/2.0;
    rad.y += (DIM0-1.0)/2.0;
    rad.z += (DIM0-1.0)/2.0;
    std::vector <double> b_norm(3);
    b_norm[0] = trilinear(&rad, bx, DIM0, DIM1, DIM2, def);
	b_norm[1] = trilinear(&rad, by, DIM0, DIM1, DIM2, def);
	b_norm[2] = trilinear(&rad, bz, DIM0, DIM1, DIM2, def);
    double norm = NORM(b_norm);*/
    double norm = 1.0;
    
    int i, j ,k;
    i = (int)p.z;
    k = (int)p.y;
    j = (int)p.x;
    cout << bx[(i*DIM0+j)*DIM0+k] << " " << by[(i*DIM0+j)*DIM0+k] << " " << bz[(i*DIM0+j)*DIM0+k] << "\n";
    
	b[0] = trilinear(&p, bx, DIM0, DIM1, DIM2, def)/norm;
	b[1] = trilinear(&p, by, DIM0, DIM1, DIM2, def)/norm;
	b[2] = trilinear(&p, bz, DIM0, DIM1, DIM2, def)/norm;
	cout << b[0] << " " << b[1] << " " << b[2] << "\n";
	// exit(1);
	return b;
}

std::vector <double> vB(double R) {
	return vB_XYZ(vR(R), R);
}

std::vector <double> vb_XYZ(std::vector <double> r, double R) {
  return NORMALIZE(vB_XYZ(r, R));
}

std::vector <double> vb (double R) {
  return NORMALIZE(vB(R));
}

std::vector <double> vB_dipole(std::vector <double> r, std::vector <double> m) {
  vector <double> n;
  n = NORMALIZE(r);
  return SUM(TIMES(3.0 * SCALAR(m, n), n), TIMES(-1.0, m));
}

std::vector <double> vb_dipole(std::vector <double> r, std::vector <double> m) {
	return NORMALIZE(vB_dipole(r, m));
}

void PrintField(){
	ofstream outputb("b_field.dat");
	int i, j, k;
	i = 10 + 128;
	j = 10 + 128;
	int l;
	std::vector <double> r(3);
	r[0] = -128;
	r[1] = 10.0;
	r[2] = 10.0;
	std::vector <double> m(3);
	double R = 0.0;
	double PHI = 20.0;
	m[0] = sin(Globals::alpha) * cos(PHI + R / Globals::RLC);
	m[1] = sin(Globals::alpha) * sin(PHI + R / Globals::RLC);
	m[2] = cos(Globals::alpha);
	double Bx, By, Bz;
    double norm = 0.0;
    std::vector <double> vR(3);
    vR[0] = -120.0;
    vR[1] = -23.0;
    vR[2] = 120.0;
	cout << vb_XYZ(vR, R)[0] << " " << vb_XYZ(vR, R)[1] << " " << vb_XYZ(vR, R)[2] << "\n";
	cout << vb_dipole(vR, m)[0] << " " << vb_dipole(vR, m)[1] << " " << vb_dipole(vR, m)[2] << "\n";
	for(l = 0; l < DIM0; l++){
		//outputb << vB_XYZ(r, R)[0]/NORM(vB_XYZ(r, R)) << " " << vb_dipole(r, m)[0] << "\n";
		outputb << bx[(i*DIM0+j)*DIM0+l] << " " << vb_dipole(r, m)[0] << "\n";
		r[0] += 1.0;
	}
	outputb.close();
}