#include <iostream>
#include <string>
#include "H5Cpp.h"
#include "aux/read_write.h"
#include "b_field.h"

using std::cout;
using std::endl;
using namespace H5;

#include <math.h>
#include <vector>

#include "initialize.h"
#include "constants.h"
#include "aux/functions.h"
#include "process_functions.h"
#include "aux/NRutil.h"

const H5std_string	FILE_NAME("fhrs.002");
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
}

std::vector <double> vB_XYZ(std::vector <double> r, double R) {
	std::vector <double> b(3);
	double def = 1.0; //?
	Point3	p;
    p.x = 25.0*r[0]+(DIM0-1.0)/2.0;
    p.y = 25.0*r[1]+(DIM0-1.0)/2.0;
    p.z = 25.0*r[2]+(DIM0-1.0)/2.0;
    p.x = p.x * cos(R/Globals::RLC) - p.y * sin(R/Globals::RLC);
    p.y = p.x * sin(R/Globals::RLC) + p.y * cos(R/Globals::RLC);
	b[0] = trilinear(&p, bx, DIM0, DIM1, DIM2, def);
	b[1] = trilinear(&p, by, DIM0, DIM1, DIM2, def);
	b[2] = trilinear(&p, bz, DIM0, DIM1, DIM2, def);
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