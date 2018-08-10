#ifndef RKFITTER_h
#define RKFITTER_h

//#define DLIB   // dlib is needed for algorithms 1 and 2

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
# include <cstdlib>
# include <iomanip>
# include <ctime>

#include <dlib/optimization.h>
#include <dlib/global_optimization.h>


#include "TkrData.h"
#include "FieldMap.h"
#include "RungeKutta4.h"

using namespace dlib;
using namespace std;

// Fit AESOP-Lite tracks to a model based on Runga-Kutta integration through the B-field map
// Robert P. Johnson       May 3, 2018
// July 6: including multiple scattering in the chi^2 calculation

class RKfitter
{
	typedef dlib::matrix<double, 0, 1> column_vector;

	bool verbose;
	FieldMap *fM;
	TkrData *tD;
	double *a;  // Track parameters
	double chiSquared;
	double *xIntercept;
	double *yIntercept;
	double *zIntercept;
	double reqmin;
	double *step;
	double stepSize;
	double sigma;
	double ynewlo;
	int icount;
	int numres;
	int ifault;
	int maxCalls;
	int *hits;
	bool multScat; 
	double SiThickness;
	double X0;    // For calculating radiation lengths in silicon
	double z0;    // Starting point in z for the track
	double **Cx;  // Covariance matrix of measurements, including multiple scattering
	void nelmin(int n, double start[], double xmin[], double *ynewlo, double reqmin, double step[], int konvge, int kcount,
		        int *icount, int *numres, int *ifault);
	void timestamp();
	RungeKutta4 *rk4;
	double chi2pm(int i, double di, int j, double dj);
	void hessian(double h[5]);
	double H[5][5];
	int invert(double **a, int n);
	double **C;
	int algorithm;
	int linearFit(double x[], double y[], int nhits, int ibp, double *a, double *b, double *c, double *xsq);

public:
	RKfitter(bool verbose, double z0, FieldMap *fM, TkrData *tD, bool multScat=true, double stepSize=5.0, int alg=0);
	double chi2(double a[]);
	double chi2_2(double S, double C);
#ifdef DLIB
	static double chi2b(column_vector b); 
	static double chi2bend(double S, double C);
// Note, to use algorithms 1 and 2 DLIB needs to be defined here and also compiled and linked.
// Plus, you must define chi2b and chi2bend as follows, in the top-level program, with rfk_global
// defined as a global variable pointer.  Then rfk_global must be set to point to the one and only
// instance of RKfitter. That is so these functions can be passed successfully to dlib routines.
/*
	RKfitter *rkf_global;
	double RKfitter::chi2b(column_vector b) {   // Version to be used by dlib
		int N = b.size();
		double a[5];
		for (int i = 0; i < N; i++) a[i] = b(i);
		//cout << "chi2b: a= " << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << " " << a[4] << endl;
		return rkf_global->chi2(a);
	}
	double RKfitter::chi2bend(double S, double C) {   // 2-D version to be used by dlib
		return rkf_global->chi2_2(S, C);
	}
*/
#endif
	double chi2m(double a[]);  // version with multiple scattering
	double chi2nm(double a[]); // version without multiple scattering
	int fitIt(bool genStartGuess, double guess[5], std::vector<int> hitSelection);
	void getIntercept(int Layer, double r[3]);
	double chiSqr() { return chiSquared; }
	void tkrParams(double prm[5]) {
		for (int i = 0; i < 5; i++) prm[i] = a[i];
	}
	void covariance(double Cov[5][5]) {
		for (int i = 0; i < 5; i++) {
			for (int j = 0; j < 5; j++) {
				Cov[i][i] = C[i][j];
			}
		}
	}
	void errors(double e[5]) {
		for (int i = 0; i < 5; i++) {
			e[i] = sqrt(C[i][i]);
		}
	}
	void print(string s);
	~RKfitter();
};

#endif