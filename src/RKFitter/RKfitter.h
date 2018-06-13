#ifndef RKFITTER_h
#define RKFITTER_h

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
# include <cstdlib>
# include <iomanip>
# include <ctime>

#include "TkrData.h"
#include "FieldMap.h"
#include "RungeKutta4.h"

using namespace std;

// Fit AESOP-Lite tracks to a model based on Runga-Kutta integration through the B-field map
// Robert P. Johnson       May 3, 2018
class RKfitter
{
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
	double z0;    // Starting point in z for the track
	void nelmin(int n, double start[], double xmin[], double *ynewlo, double reqmin, double step[], int konvge, int kcount,
		        int *icount, int *numres, int *ifault);
	void timestamp();
	RungeKutta4 *rk4;
	double chi2pm(int i, double di, int j, double dj);
	void hessian(double h[5]);
	double H[5][5];
	int invert(double **a, int n);
	double **C;
	int linearFit(double x[], double y[], int nhits, int ibp, double *a, double *b, double *c, double *xsq);

public:
	RKfitter(bool verbose, double z0, FieldMap *fM, TkrData *tD);
	double chi2(double a[]);
	int fitIt(bool genStartGuess, double guess[5], vector<int> hitSelection);
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