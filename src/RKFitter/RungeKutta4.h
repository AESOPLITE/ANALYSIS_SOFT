#ifndef RUNGEKUTTA4_h
#define RUNGEKUTTA4_h

#include <random>
#include "FieldMap.h"

using namespace std;

// Perform a 4th order Runge-Kutta numerical integration of the Lorentz force through a B field map
// Robert P. Johnson       May 3, 2018
class RungeKutta4 {
	double h, h2;   // The step size along the arc and its square
	double alpha;
	FieldMap *fM;
	double *f(double x[3], double p[3]);
	double *xA;
	double *yA;
	double *zA;
	double *sA;
	double *pX;
	double *pY;
	double *pZ;
	int nScat;
	double *zScat;
	double *scatAng;
	int aSize;  // Current size of the arrays
	int nStep;
	double X0;
	int IndexNow;

public:
	RungeKutta4(double dz, FieldMap *fM);
	RungeKutta4(double dz, FieldMap *fM, int nScat, double zScat[]);
	double *Integrate(double Q, double r0[3], double p0[3], double s, double deltaZ);
	double *Integrate(double Q, double r0[3], double p0[3], double s, double deltaZ, double ranNorm[]);
	double *getX(int i);
	double *getP(int i);
	double *getP();
	double *getX(double z, bool *flag);
	double getS(int i);
	double getScat(int i);
	int getNsteps();
	~RungeKutta4();
};

#endif