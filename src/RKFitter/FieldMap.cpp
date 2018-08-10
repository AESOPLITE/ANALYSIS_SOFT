#include "FieldMap.h"

// Load the AESOP-Lite B-field map and interpolate it
// If no filename is provided or Type is not "binary" or "text", then the map can be filled by the methods
// setOriginStep and addEntry.
// If Type is "text" the map will be (slowly) read from the original text file.  In that case, the writeBinaryFile method
// can be called to create a (machine-dependent) binary file. Then use Type="binary" to do a fast binary read.
// Note: for the 5mm map, nPoints = 81, dz=5, and x0=-200
FieldMap::FieldMap(string FileName, string Type, int nPoints)
{
	cout << ">>>>>>>> FieldMap constructor called for " << nPoints << " points per dimension <<<<<<<<<<<" << endl;
	nFldPnt = nPoints;
	this->FileName = FileName;
	this->Type = Type;
	fMapX = NULL;
	fMapY = NULL;
	fMapZ = NULL;
	size_t found = FileName.find("1mm");
	double xUnit = 1.0;
	double bUnit = 10000.;
	if (found != string::npos) {
		xUnit = 1000.;
		bUnit = 1.;
	}
	if (FileName == " " || (Type != "binary" && Type != "text")) {
		cout << "FieldMap: no input file is to be read. Field map must be filled by function call." << endl;
		nEntered = 0;
		int size = nFldPnt*nFldPnt*nFldPnt;
		fMapX = new float[size];
		cout << "FieldMap: making an array for Bx of size " << size << endl;
		fMapY = new float[size];
		cout << "FieldMap: making an array for By of size " << size << endl;
		fMapZ = new float[size];
		cout << "FieldMap: making an array for Bz of size " << size << endl;
		return;
	}
	zOffset = -106.3409;
	cout << "FieldMap: z position of the magnet center is set to " << zOffset << endl;
	nEntered = nFldPnt*nFldPnt*nFldPnt;
	cout << "Loading the Field Map from " << FileName << "; zOffset=" << zOffset << endl;

	if (Type == "binary") { // The fast, but machine-dependent option.  See code below to write out a suitable binary file.
		cout << "    Reading the field map from a binary file." << endl;
		int size = nFldPnt*nFldPnt*nFldPnt;
		fMapX = new float[size];
		cout << "FieldMap: making an array for Bx of size " << size << endl;
		fMapY = new float[size];
		cout << "FieldMap: making an array for By of size " << size << endl;
		fMapZ = new float[size];
		cout << "FieldMap: making an array for Bz of size " << size << endl;
		ifstream infile;
		infile.open(FileName, ios::binary | ios::in);
		int nPrint = 0;
		if (infile.is_open()) {
			infile.read((char*)&x0, sizeof(double));
			infile.read((char*)&y0, sizeof(double));
			infile.read((char*)&z0, sizeof(double));
			infile.read((char*)&dx, sizeof(double));
			infile.read((char*)&dy, sizeof(double));
			infile.read((char*)&dz, sizeof(double));
			for (int i = 0; i < nFldPnt; i++) {
				for (int j = 0; j < nFldPnt; j++) {
					for (int k = 0; k < nFldPnt; k++) {
						infile.read((char*)(fMapX + i*nFldPnt*nFldPnt + j*nFldPnt + k), sizeof(float));
						infile.read((char*)(fMapY + i*nFldPnt*nFldPnt + j*nFldPnt + k), sizeof(float));
						infile.read((char*)(fMapZ + i*nFldPnt*nFldPnt + j*nFldPnt + k), sizeof(float));
						if (nPrint < 20) {
							nPrint++;
							cout << "i=" << i << " j=" << j << " k=" << k << " Bx=" << fMapX[i*nFldPnt*nFldPnt + j*nFldPnt + k]
								<< " By=" << fMapY[i*nFldPnt*nFldPnt + j*nFldPnt + k] << " Bz=" << fMapZ[i*nFldPnt*nFldPnt + j*nFldPnt + k] << endl;
						}
					}
				}
			}
			infile.close();
			cout << "FieldMap: x0=" << x0 << "  dx=" << dx << endl;
			cout << "FieldMap: y0=" << y0 << "  dy=" << dy << endl;
			cout << "FieldMap: z0=" << z0 << "  dz=" << dz << endl;
		}
		else {
			cout << "failed to open the field map" << endl;
		}
	}
	else { // Option "text" for reading the field map from a text file. Machine-independent but quite slow.
		cout << "    Reading the field map from a text file." << endl;
		ifstream mapFile(FileName);
		if (mapFile.is_open()) {
			string line;
			getline(mapFile, line);
			cout << "FieldMap line 1: " << line << endl;
			getline(mapFile, line);
			cout << "FieldMap line 2: " << line << endl;
			int lnCnt = 0;
			int size = nFldPnt*nFldPnt*nFldPnt;
			fMapX = new float[size];
			cout << "FieldMap: making an array for Bx of size " << size << endl;
			fMapY = new float[size];
			cout << "FieldMap: making an array for By of size " << size << endl;
			fMapZ = new float[size];
			cout << "FieldMap: making an array for Bz of size " << size << endl;
			for (int i = 0; i < nFldPnt; i++) {
				for (int j = 0; j < nFldPnt; j++) {
					for (int k = 0; k < nFldPnt; k++) {
						if (mapFile.good()) {
							getline(mapFile, line);
							//cout << line << endl;
							lnCnt++;
							float x, y, z, Bx, By, Bz;
							sscanf(line.c_str(), "%f %f %f %f %f %f", &x, &y, &z, &Bx, &By, &Bz);
							//cout << "Bx= " << Bx << " By= " << By << " Bz= " << Bz << endl;
							fMapX[i*nFldPnt*nFldPnt + j*nFldPnt + k] = (float)Bx/bUnit;
							fMapY[i*nFldPnt*nFldPnt + j*nFldPnt + k] = (float)By/bUnit;
							fMapZ[i*nFldPnt*nFldPnt + j*nFldPnt + k] = (float)Bz/bUnit;
							x *= xUnit;
							y *= xUnit;
							z *= xUnit;
							if (i == 0) x0 = x;
							if (j == 0) y0 = y;
							if (k == 0) z0 = z;
							if (i == 1) dx = x - x0;
							if (j == 1) dy = y - y0;
							if (k == 1) dz = z - z0;
						}
						else {
							cout << "Field map " << FileName << " is incomplete " << endl;
							return;
						}
					}
				}
				//cout << "FieldMap: " << i << " " << fMapX[i*nFldPnt*nFldPnt] << " " << fMapY[i*nFldPnt*nFldPnt] << " " << fMapZ[i*nFldPnt*nFldPnt] << endl;
			} 
			cout << "FieldMap: read " << lnCnt << " data lines from the field map text file " << FileName << endl; 
			cout << "FieldMap: x0=" << x0 << "  dx=" << dx << endl;
			cout << "FieldMap: y0=" << y0 << "  dy=" << dy << endl;
			cout << "FieldMap: z0=" << z0 << "  dz=" << dz << endl;
		}  
	}  
}

void FieldMap::setOriginStep(double x0, double dx, double zOffset) {
	this->x0 = x0;
	this->y0 = x0;
	this->z0 = x0;
	this->dx = dx;
	this->dy = dx;
	this->dz = dx;
	this->zOffset = zOffset;
	cout << "FieldMap::setOriginStep: x0=" << this->x0 << "  dx=" << this->dx << endl;
	cout << "FieldMap::setOriginStep: y0=" << y0 << "  dy=" << dy << endl;
	cout << "FieldMap::setOriginStep: z0=" << z0 << "  dz=" << dz << endl;
	cout << "FieldMap::setOriginStep: z position of the magnet center is set to " << zOffset << endl;
}

void FieldMap::addEntry(float Bx, float By, float Bz) {
	if (nEntered >= nFldPnt*nFldPnt*nFldPnt) {
		cout << "FieldMap::addEntry: error---no more space to enter field points" << endl;
		return;
	}
	fMapX[nEntered] = Bx;
	fMapY[nEntered] = By;
	fMapZ[nEntered] = Bz;
	nEntered++;
	if (nEntered == nFldPnt*nFldPnt*nFldPnt) cout << "FieldMap::addEntry: the B field map is now full" << endl;
}

float *FieldMap::getEntry(int i) {
	float *B = new float[3];
	if (i > nFldPnt) return B;
	B[0] = fMapX[i];
	B[1] = fMapY[i];
	B[2] = fMapZ[i];
	return B;
}

 double *FieldMap::GetField(double x, double y, double z) {
	if (nEntered != nFldPnt*nFldPnt*nFldPnt) cout << "FieldMap::GetField WARNING: the field map is still incomplete!" << endl;
	double zmag = z - zOffset;
	double xmag = x;
	double ymag = y;
	int i = (int)floor((xmag - x0) / dx);
	if (i < 0) {
		i = 0;
		xmag = x0;
	}
	if (i > nFldPnt - 2) {
		i = nFldPnt - 2;
		xmag = x0 + nFldPnt*dx;
	}
	int j = (int)floor((ymag - y0) / dy);
	if (j < 0) {
		j = 0;
		ymag = y0;
	}
	if (j > nFldPnt - 2) {
		j = nFldPnt - 2;
		ymag = y0 + nFldPnt*dy;
	}
	int k = (int)floor((zmag - z0) / dz);
	if (k < 0) {
		k = 0;
		zmag = z0;
	}
	if (k > nFldPnt - 2) {
		k = nFldPnt - 2;
		zmag = z0 + nFldPnt*dz;
	}
	double xd = (xmag - (x0+i*dx)) / dx;
	double yd = (ymag - (y0+j*dy)) / dy;
	double zd = (zmag - (z0+k*dz)) / dz;
	double *B = new double[3];
	B[0] = TriLinear(i, j, k, xd, yd, zd, fMapX);
	B[1] = TriLinear(i, j, k, xd, yd, zd, fMapY);
	B[2] = TriLinear(i, j, k, xd, yd, zd, fMapZ);
	return B;
}

double FieldMap::TriLinear(int i, int j, int k, double xd, double yd, double zd, float *f) {
	double c00 = f[i*nFldPnt*nFldPnt + j*nFldPnt + k] * (1.0 - xd) + f[(i+1)*nFldPnt*nFldPnt + j*nFldPnt + k] * xd; // interpolate in x
	double c01 = f[i*nFldPnt*nFldPnt + j*nFldPnt + k + 1] * (1.0 - xd) + f[(i+1)*nFldPnt*nFldPnt + j*nFldPnt + k + 1] * xd;
	double c10 = f[i*nFldPnt*nFldPnt + (j+1)*nFldPnt + k] * (1.0 - xd) + f[(i+1)*nFldPnt*nFldPnt + (j+1)*nFldPnt + k] * xd;
	double c11 = f[i*nFldPnt*nFldPnt + (j+1)*nFldPnt + k+1] * (1.0 - xd) + f[(i+1)*nFldPnt*nFldPnt + (j+1)*nFldPnt + k+1] * xd;
	double c0 = c00 * (1.0 - yd) + c10 * yd; // interpolate in y
	double c1 = c01 * (1.0 - yd) + c11 * yd;
	double c = c0 * (1.0 - zd) + c1 * zd; // interpolate in z
	return c;
}

void FieldMap::writeBinaryFile(string fName) { // Make a binary field map file that can be read much more quickly
	ofstream outfile;
	outfile.open(fName, ios::binary | ios::out);
	outfile.write((char*)&x0, sizeof(double));
	outfile.write((char*)&y0, sizeof(double));
	outfile.write((char*)&z0, sizeof(double));
	outfile.write((char*)&dx, sizeof(double));
	outfile.write((char*)&dy, sizeof(double));
	outfile.write((char*)&dz, sizeof(double));
	for (int i = 0; i < nFldPnt; i++) {
		for (int j = 0; j < nFldPnt; j++) {
			for (int k = 0; k < nFldPnt; k++) {
				outfile.write((char*)(fMapX + i*nFldPnt*nFldPnt + j*nFldPnt + k), sizeof(float));
				outfile.write((char*)(fMapY + i*nFldPnt*nFldPnt + j*nFldPnt + k), sizeof(float));
				outfile.write((char*)(fMapZ + i*nFldPnt*nFldPnt + j*nFldPnt + k), sizeof(float));
			}
		}
	}
	outfile.close();
}

FieldMap::~FieldMap()
{
	cout << ">>>>>>>> FieldMap destructor called <<<<<<<<<<<<<" << endl;
	delete[] fMapX;
	fMapX = NULL;
	delete[] fMapY;
	fMapY = NULL;
	delete[] fMapZ;
	fMapZ = NULL;
}
