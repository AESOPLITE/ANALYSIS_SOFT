#ifndef RUNGEKUTTA4_h
#define RUNGEKUTTA4_h

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
	int aSize;  // Current size of the arrays
	int nStep;

public:
	RungeKutta4(double dz, FieldMap *fM);
	double *Integrate(double Q, double r0[3], double p0[3], double s, double deltaZ);
	double *getX(int i);
	double *getX(double z, bool *flag);
	double getS(int i);
	int getNsteps();
	~RungeKutta4();
};

#endif RUNGEKUTTA4_h