#ifndef FIELDMAP_h
#define FIELDMAP_h

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <string>

using namespace std;

// Read the AESOP B field map from a file (text or binary) and linearly interpolate it in 3-D
// Robert P. Johnson       May 3, 2018
class FieldMap {
	int nFldPnt;
	int nEntered;
	string FileName;
	string Type;
	float *fMapX;     // X component of the field in Tesla
	float *fMapY;     // Y component
	float *fMapZ;     // Z component
	double zOffset;   // Location of the center of the magnet in z
	double dx, dy, dz;
	double x0, y0, z0;
	double TriLinear(int i, int j, int k, double xd, double yd, double zd, float *f);

public:
	FieldMap(string FileName = " ", string Type = "none", int nPoints = 81);
	double *GetField(double x, double y, double z);
	void setOriginStep(double x0 = -200., double dx = 5., double zOffset = -106.3409);
	void addEntry(float Bx, float By, float Bz);
	void writeBinaryFile(string fName);
	float *getEntry(int i);
	void reset() {
		nEntered = 0;
	}
	~FieldMap();
};

#endif