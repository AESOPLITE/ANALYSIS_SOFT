#include "RKfitter.h"

// Calculate the chi^2 of the fit for a given track. This is what we try to minimize.
double RKfitter::chi2(double a[]) {
	//cout << "RKfitter::chi2: entering with a=" << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << " " << a[4] << endl;
	// estimate how far to integrate, to cover all layers
	double arg = 1.0 - a[2] * a[2] - a[3] * a[3];
	if (arg <= 0.) {
		return 9.9e9;
	}
	double Delta_z = z0 - (tD->zLayer[tD->nLayers - 1]);
	double ctz = -sqrt(arg);
	double s = -Delta_z / ctz + 200.0;
	if (s > 500. || s < 0.) s = 500.;

	// Set up the input parameters needed by the integrator
	double Q;
	if (a[4] < 0.) Q = -1.0; else Q = 1.0;
	double p = Q / a[4];
	double r0[3] = { a[0], a[1], z0 };
	double p0[3] = { p*a[2], p*a[3], p*ctz };

	// Integrate all the way through the instrument
	if (verbose) {
		cout << "RKfitter::chi2: r0=" << r0[0] << " " << r0[1] << " " << r0[2];
		cout << "   p0=" << p0[0] << " " << p0[1] << " " << p0[2] << "  s=" << s << endl;
	}
	double *xEnd= rk4->Integrate(Q, r0, p0, s, abs(Delta_z)+40.);
	if (verbose) cout << "completed integration at xEnd=" << xEnd[0] << " " << xEnd[1] << " " << xEnd[2]  << endl;
	delete[] xEnd;

	// Loop over layers and add up the chi^2
	double result = 0.;
	for (int lyr = 0; lyr < tD->nLayers; lyr++) {
		if (tD->hits[lyr].size() == 0) continue;  // skip if no hits on the layer
		if (hits[lyr] < 0) continue;              // allow user to skip layer
		bool flag;
		if (tD->orientation[lyr] == 'n') {
			double *rInterp = rk4->getX(tD->zLayer[lyr], &flag);
			double xInterp = rInterp[0];
			if (!flag) {
				xInterp = 220.;   // force a large chi^2 if the track doesn't even intercept the layer
				xIntercept[lyr] = -9999.;
				yIntercept[lyr] = -9999.;
			} else {
				xIntercept[lyr] = rInterp[0];
				yIntercept[lyr] = rInterp[1];
			}
			double xMeas = tD->hits[lyr].at(hits[lyr]);
			double incr = pow((xInterp - xMeas) / sigma, 2);
			result += incr;	
			delete[] rInterp;
			//cout << "   Layer " << lyr << " flag=" << flag << " xInterp=" << xInterp << " xMeas=" << xMeas << " chi2_inc=" << incr << endl;
		}
		else {
			double *rInterp = rk4->getX(tD->zLayer[lyr], &flag);
			double yInterp = rInterp[1];
			if (!flag) {
				yInterp = 220.;
				xIntercept[lyr] = -9999.;
				yIntercept[lyr] = -9999.;
			} else {
				xIntercept[lyr] = rInterp[0];
				yIntercept[lyr] = rInterp[1];
			}
			double yMeas = tD->hits[lyr].at(hits[lyr]);
			double incr = pow((yInterp - yMeas) / sigma, 2);
			result += incr;
			delete[] rInterp;
			//cout << "   Layer " << lyr << " flag=" << flag << " yInterp=" << yInterp << " yMeas=" << yMeas << " chi2_inc=" << incr << endl;
		}
	}
	if (verbose) {
		cout << "RKfitter::chi2: with sigma = " << sigma << "  a=";
		for (int i = 0; i < 5; i++) cout << a[i] << " ";
		cout << "  chi^2=" << result << endl;
	}
	return result;
}

// Calculate numerically the 2nd derivative matrix of the chi^2 function
void RKfitter::hessian(double hin[5]) {
	double fxy = chi2(a);
	double h[5];
	for (int i = 0; i < 5; i++) {
		h[i] = hin[i] / 2.;
		double delta;
		do {
			h[i] *= 2.;
			delta = (chi2pm(i, h[i], i, 0.)-fxy)/fxy;
		} while (abs(delta) < 0.00001);
		//cout << "RKfitter::hessian: h[" << i << "]=" << h[i] << endl;
	}
	for (int i = 0; i < 5; i++) {
		double fxpy = chi2pm(i, h[i], i, 0.);
		double fxmy = chi2pm(i, -h[i], i, 0.);
		H[i][i] = (fxpy - 2.*fxy + fxmy) / (h[i] * h[i]);
		for (int j = i + 1; j < 5; j++) {
			double fxpyp = chi2pm(i, h[i], j, h[j]);
			double fxmym = chi2pm(i, -h[i], j, -h[j]);
			double fxpym = chi2pm(i, h[i], j, -h[j]);
			double fxmyp = chi2pm(i, -h[i], j, h[j]);
			H[i][j] = (fxpyp - fxpym - fxmyp + fxmym) / (4.0 * h[i] * h[j]);
		}
	}
	for (int i=1; i<5; i++) {
		for (int j = 0; j < i; j++) {
			H[i][j] = H[j][i];
		}
	}
}

// Utility routine to facilitate the hessian calculation
double RKfitter::chi2pm(int i, double di, int j, double dj) {
	double b[5];
	for (int k = 0; k < 5; k++) {
		if (k == i) b[k] = a[k] + di;
		else if (k == j) b[k] = a[k] + dj;
		else b[k] = a[k];
	}
	return chi2(b);
}

RKfitter::RKfitter(bool verbose, double z0,  FieldMap *fM, TkrData *tD) {
	this->verbose = verbose;
	this->z0 = z0;
	this->fM = fM;
	this->tD = tD;
	hits = new int[tD->nLayers];
	xIntercept = new double[tD->nLayers];
	yIntercept = new double[tD->nLayers];
	zIntercept = new double[tD->nLayers];
	for (int lyr = 0; lyr < tD->nLayers; lyr++) {
		zIntercept[lyr] = tD->zLayer[lyr];
	}
	maxCalls = 800;   // Maximum function calls allowed in the minimization search
	reqmin = 0.0001;    // convergence check parameter
	step = new double[5];
	step[0] = 0.3;    // initial step for x position in the minimization search
	step[1] = 0.3;    // initial step for y position
	step[2] = 0.005;  // initial step size for x direction cosine
	step[3] = 0.005;  // initial step size for y direction cosine
	step[4] = 5.0;    // initial step size for 1/p as a percentage
	stepSize = 5.0;   // Runge-Kutta integration step size (comparable to the B field map precision)
	rk4 = new RungeKutta4(stepSize, fM);
	sigma = tD->stripPitch / sqrt(12.0);  // Assumed resolution of the AESOP tracker
	a = new double[5];
	C = new double*[5];
	for (int i = 0; i < 5; i++) {
		C[i] = new double[5];
	}
}

// Execute the brute-force track fit using a non-derivative simplex minimization method
void RKfitter::fitIt(bool genStartGuess, double guess[5], vector<int> hitSelection) {
	// hitSelection specifies which hit to use from each layer of the TkrData structure
	// guess[5] is the externally supplied starting guess for the track
	// if genStartGuess is true, then the program will generate an initial guess from a linear fit

	if (hitSelection.size() != tD->nLayers) cout << "RKfitter:fitIt, wrong number of hits specified, need " << tD->nLayers << endl;
	for (int lyr = 0; lyr < tD->nLayers; lyr++) hits[lyr] = hitSelection[lyr];
	double *temp = new double[5];

	// estimate some reasonable values for the initial step in each of the parameters
	// and copy the initial guess, since the minimization routine may overwrite it
	double initStep[5];
	for (int i = 0; i < 5; i++) {
		temp[i] = guess[i];
		initStep[i] = step[i];
	}
	initStep[4] = initStep[4] * a[4] / 100.;

	// if requested, cook up an initial guess here, instead of using the one supplied
	if (verbose) cout << "RKfitter::fitIt: supplied initial guess= " << temp[0] << " " << temp[1] << " " << temp[2] << " " << temp[3] << " " << temp[4] << endl;
	if (genStartGuess) {
		if (verbose) cout << "Generate the starting guess from line and parabola fits to " << tD->nLayers << " tracker layers" << endl;
		double *zn = new double[tD->nLayers];
		double *zb = new double[tD->nLayers];
		double *x = new double[tD->nLayers]; 
		double *y = new double[tD->nLayers];
		int Nn = 0;
		int Nb = 0;
		for (int lyr = 0; lyr < tD->nLayers; lyr++) {
			if (verbose) cout << "    Layer " << lyr;
			if (tD->orientation[lyr] == 'n') {
				zn[Nn] = tD->zLayer[lyr];
				x[Nn] = tD->hits[lyr].at(hits[lyr]);
				if (verbose) cout << " Nn=" << Nn << "  z=" << zn[Nn] << "  x=" << x[Nn] << endl;
				Nn++;
			}
			else {
				zb[Nb] = tD->zLayer[lyr];
				y[Nb] = tD->hits[lyr].at(hits[lyr]);
				if (verbose) cout << " Nb=" << Nb << "  z=" << zb[Nb] << "  y=" << y[Nb] << endl;
				Nb++;
			}
		}  
		double a, b, c, d, e, xsqn, xsqb;
		if (linearFit(zn, x, Nn, 0, &a, &b, &c, &xsqn) == 0) {
			if (linearFit(zb, y, Nb, 1, &c, &d, &e, &xsqb) == 0) {
				if (verbose) cout << "RKfitter::fitIt: linear parameters: " << a << " " << b << " " << c << " " << d << " " << e << endl;
				double clight = 2.99793e8; // Speed of light in m/s
				double *Bf = fM->GetField(temp[0], temp[1], tD->zLayer[2]);
				double B = sqrt(Bf[0] * Bf[0] + Bf[1] * Bf[1] + Bf[2] * Bf[2]);
				double alpha = 1.0e12 / (clight * B); 
				//cout << "RKfitter::fitIt: B= " << Bf[0] << " " << Bf[1] << " " << Bf[2] << endl;
				if (verbose) cout << "Rkfitter::fitIt: B field in center of detector = " << B << " tesla,    alpha=" << alpha << endl;
				temp[0] = -(a + b*z0);
				temp[1] = -(c + z0*(d + z0*e));
				double slp2 = d + 2.*e*z0;
				temp[2] = -b / sqrt(1. + b*b + slp2*slp2);
				temp[3] = -slp2 / sqrt(1. + b*b + slp2*slp2);
				temp[4] = -1.25*alpha*2.0*e*sqrt(1.0 - temp[2] * temp[2]);
				if (verbose) cout << "RKfitter::fitIt: initial guess=" << temp[0] << " " << temp[1] << " " << temp[2] << " " << temp[3] << " " << temp[4] << endl;
				delete[] Bf;
			} else {
				cout << "RKfitter::fitIt: parabola fit failed when getting the starting guess" << endl;
			}
		}
		else {
			cout << "RKfitter::fitIt: line fit failed when getting the starting guess" << endl;
		}
		delete[] zn;
		delete[] zb;
		delete[] x;
		delete[] y;
	}

	int nVar = 5;     // Number of parameters being fit
	int Konvge = 5;   // How often to check for convergence

	// call the canned simplex minimization routine
	nelmin(nVar, temp, a, &ynewlo, reqmin, initStep, Konvge, maxCalls, &icount, &numres, &ifault);
	//cout << "RKfitter::fitIt: from nelmin, newlo=" << ynewlo << " count=" << icount << " numres=" << numres << " error=" << ifault << endl;

	// Propagation of errors:
	double h[5] = { 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001*abs(a[4]) };
	hessian(h);
	//double T[5][5];
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			C[i][j] = 0.5*H[i][j];
			//if (i != j) C[i][j] = 0.;
			//T[i][j] = C[i][j];
		}
	}

	invert(C,5);  // Invert the hessian matrix (times -0.5) to get the covariance matrix

/*
    // Check the matrix inversion
	double U[5][5];
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			U[i][j] = 0.;
			for (int k = 0; k < 5; k++) {
				U[i][j] += C[i][k] * T[k][j];
			}
		}
	}
	if (verbose) {
		cout << "     Unit Matrix:" << endl;
		for (int i = 0; i < 5; i++) {
			cout << "         ";
			for (int j = 0; j < 5; j++) {
				cout << U[i][j] << "   ";
			}
			cout << endl;
		}
	}
*/
	// Calculate the chi^2 once more in order to get the final residuals
	chiSquared = chi2(a);
	if (verbose) {
		cout << "RKfitter::fitIt: final fit chi^2 = " << chiSquared << endl;
		for (int lyr = 0; lyr < tD->nLayers; lyr++) {
			double r[3];
			getIntercept(lyr, r);
			double residual;
			if (tD->orientation[lyr] == 'n') {
				residual = r[0] - tD->hits[lyr].at(0);
			} else {
				residual = r[1] - tD->hits[lyr].at(0);
			}
			cout << "    layer " << lyr << "  residual = " << residual << endl;
		}
	}

	delete[] temp;
}

// Return the position at each silicon plane interpolated from the Runge-Kutta integration
void RKfitter::getIntercept(int Layer, double r[3]) {
	r[0] = xIntercept[Layer];
	r[1] = yIntercept[Layer];
	r[3] = zIntercept[Layer];
}

void RKfitter::print(string s) {
	cout << endl;
	cout << "RKfitter fit results for hits ";
	for (int i = 0; i < tD->nLayers; i++) cout << hits[i] << " ";
	cout << endl;
	cout << "    Fitted track parameters at z=" << z0 << ": x=" << a[0] << " y=" << a[1] << " ctx=" << a[2] << " cty=" << a[3] << " Q/p=" << a[4] << endl;
	double p = abs(1 / a[4]);
	cout << "    Fitted momentum = " << p << " GeV" << endl;
	cout << "    Minimum chi^2=" << ynewlo << " in " << icount << " function calls and " << numres << " restarts, with " << ifault << " errors" << endl;
	cout << "    Hessian matrix:" << endl;
	for (int i = 0; i < 5; i++) {
		cout << "         ";
		for (int j = 0; j < 5; j++) {
			cout << H[i][j] << "   ";
		}
		cout << endl;
	}
	cout << "     Covariance Matrix:" << endl;
	for (int i = 0; i < 5; i++) {
		cout << "         ";
		for (int j = 0; j < 5; j++) {
			cout << C[i][j] << "   ";
		}
		cout << endl;
	}
	cout << "     Errors: ";
	for (int i = 0; i < 5; i++) {
		cout << sqrt(C[i][i]) << "  ";
	}
	cout << endl << endl;
}

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
void RKfitter::invert(double **a, int n) { // code slightly modified from Numerical Recipes in C "gaussj"
	int i, icol, irow, j, k, l, ll;
	double big, dum, pivinv, temp;
	int *indxc = new int[n];
	int *indxr = new int[n];
	int *ipiv = new int[n];
	for (j = 0; j<n; j++) ipiv[j] = 0;
	for (i = 0; i<n; i++) {
		big = 0.0;
		for (j = 0; j<n; j++)
			if (ipiv[j] != 1)
				for (k = 0; k<n; k++) {
					if (ipiv[k] == 0) {
						if (abs(a[j][k]) >= big) {
							big = abs(a[j][k]);
							irow = j;
							icol = k;
						}
					}
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l = 0; l<n; l++) SWAP(a[irow][l], a[icol][l]);
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if (a[icol][icol] == 0.0) throw("gaussj: Singular Matrix");
		pivinv = 1.0 / a[icol][icol];
		a[icol][icol] = 1.0;
		for (l = 0; l<n; l++) a[icol][l] *= pivinv;
		for (ll = 0; ll<n; ll++)
			if (ll != icol) {
				dum = a[ll][icol];
				a[ll][icol] = 0.0;
				for (l = 0; l<n; l++) a[ll][l] -= a[icol][l] * dum;
			}
	}
	for (l = n - 1; l >= 0; l--) {
		if (indxr[l] != indxc[l])
			for (k = 0; k<n; k++)
				SWAP(a[k][indxr[l]], a[k][indxc[l]]);
	}

	delete[] ipiv;
	delete[] indxr;
	delete[] indxc;
}
#undef SWAP

// Line or parabola fit translated from John Clem's FORTRAN program "find_track"
int RKfitter::linearFit(double x[], double y[], int nhits, int ibp, double *a, double *b, double *c, double *xsq) {

	// x, y    input : tracker data in either bending or non - bending
	// nhits   input : number of tracker points in either bending or non - bending
	// ibp     input : indicator flag 1 = bending plane data, otherwise non - bending
	// a, b, c output least square fit to either straight line y = a + bx or 2nd order polynomial y = a + b*x + cx ^ 2
	// xsq     output Total Chi Square

	// Initialization
	double sum = 0.0;
	double sumx = 0.0;
	double sumxy = 0.0;
	double sumx2y = 0.0;
	double sumy = 0.0;
	double sumx2 = 0.0;
	double sumx3 = 0.0;
	double sumx4 = 0.0;

	// Matrix Elements
	for (int ii = 0; ii < nhits; ii++) {
		sum = sum + 1.0;
		sumx = sumx + x[ii];
		sumx2 = sumx2 + x[ii] * x[ii];
		sumx3 = sumx3 + x[ii] * x[ii] * x[ii];
		sumx4 = sumx4 + x[ii] * x[ii] * x[ii] * x[ii];
		sumy = sumy + y[ii];
		sumxy = sumxy + x[ii] * y[ii];
		sumx2y = sumx2y + x[ii] * x[ii] * y[ii];
	}

	if (ibp == 1) { // second order y = a + b*x + cx ^ 2

		double det = sum*sumx2*sumx4 + sumx*sumx3*sumx2 + sumx2*sumx*sumx3 - (sumx2*sumx2*sumx2 + sumx3*sumx3*sum + sumx4*sumx*sumx);
		double a1 = sumy*sumx2*sumx4 + sumx*sumx3*sumx2y + sumx2*sumxy*sumx3 - (sumx2y*sumx2*sumx2 + sumx3*sumx3*sumy + sumx4*sumxy*sumx);
		double a2 = sum*sumxy*sumx4 + sumy*sumx3*sumx2 + sumx2*sumx*sumx2y - (sumx2*sumxy*sumx2 + sumx2y*sumx3*sum + sumx4*sumx*sumy);
		double a3 = sum*sumx2*sumx2y + sumx*sumxy*sumx2 + sumy*sumx*sumx3 - (sumx2*sumx2*sumy + sumx3*sumxy*sum + sumx2y*sumx*sumx);
		if (det == 0.) {
			cout << "linearFit: zero determinate in the parabola fit.  Abort!" << endl;
			return -1;
		}
		*a = a1 / det;
		*b = a2 / det;
		*c = a3 / det;
		*xsq = 0.0;
		for (int ii = 0; ii < nhits; ii++) {
			*xsq += pow(y[ii] - *a - *b*x[ii] - *c*x[ii] * x[ii], 2);
		}
		return 0;
	}

	// linear y = a + bx
	double det = sum*sumx2 - sumx*sumx;
	if (det == 0.) {
		cout << "linearFit: zero determinate in the line fit.  Abort!" << endl;
		return -1;
	}
	*a = (sumx2*sumy - sumx*sumxy) / det;
	*b = (sum*sumxy - sumx*sumy) / det;
	*c = 0.;
	*xsq = 0.0;
	for (int ii = 0; ii < nhits; ii++) {
		*xsq += pow(y[ii] - *a - *b*x[ii], 2);
	}
	return 0;
}

RKfitter::~RKfitter()
{
	//cout << ">>>>>> RKfitter cleaning up <<<<<<<<<<<<" << endl;
	delete rk4;
	delete[] a;
	delete[] step;
	delete[] hits;
	delete[] xIntercept;
	delete[] yIntercept;
	delete[] zIntercept;
	for (int i = 0; i < 5; i++) delete[] C[i];
	delete[] C;
}

// Nelder-Mead algorithm from the Web, made here into a member function for convenience
void RKfitter::nelmin(int n, double start[], double xmin[],	double *ynewlo, double reqmin, double step[], int konvge, int kcount,
	                  int *icount, int *numres, int *ifault)

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NELMIN minimizes a function using the Nelder-Mead algorithm.
	//
	//  Discussion:
	//
	//    This routine seeks the minimum value of a user-specified function.
	//
	//    Simplex function minimisation procedure due to Nelder+Mead(1965),
	//    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
	//    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
	//    25, 97) and Hill(1978, 27, 380-2)
	//
	//
	//    This routine does not include a termination test using the
	//    fitting of a quadratic surface.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license. 
	//
	//  Modified:
	//
	//    27 February 2008
	//
	//  Author:
	//
	//    Original FORTRAN77 version by R ONeill.
	//    C++ version by John Burkardt.
	//
	//  Reference:
	//
	//    John Nelder, Roger Mead,
	//    A simplex method for function minimization,
	//    Computer Journal,
	//    Volume 7, 1965, pages 308-313.
	//
	//    R ONeill,
	//    Algorithm AS 47:
	//    Function Minimization Using a Simplex Procedure,
	//    Applied Statistics,
	//    Volume 20, Number 3, 1971, pages 338-345.
	//
	//  Parameters:
	//
	//    Input, double FN ( double x[] ), the name of the routine which evaluates
	//    the function to be minimized.  (I removed this argument and embedded calls to member function chi2---R.P.J.)
	//
	//    Input, int N, the number of variables.
	//
	//    Input/output, double START[N].  On input, a starting point
	//    for the iteration.  On output, this data may have been overwritten.
	//
	//    Output, double XMIN[N], the coordinates of the point which
	//    is estimated to minimize the function.
	//
	//    Output, double YNEWLO, the minimum value of the function.
	//
	//    Input, double REQMIN, the terminating limit for the variance
	//    of function values.
	//
	//    Input, double STEP[N], determines the size and shape of the
	//    initial simplex.  The relative magnitudes of its elements should reflect
	//    the units of the variables.
	//
	//    Input, int KONVGE, the convergence check is carried out 
	//    every KONVGE iterations.
	//
	//    Input, int KCOUNT, the maximum number of function 
	//    evaluations.
	//
	//    Output, int *ICOUNT, the number of function evaluations 
	//    used.
	//
	//    Output, int *NUMRES, the number of restarts.
	//
	//    Output, int *IFAULT, error indicator.
	//    0, no errors detected.
	//    1, REQMIN, N, or KONVGE has an illegal value.
	//    2, iteration terminated because KCOUNT was exceeded without convergence.
	//
{
	double ccoeff = 0.5;
	double del;
	double dn;
	double dnn;
	double ecoeff = 2.0;
	double eps = 0.001;
	int i;
	int ihi;
	int ilo;
	int j;
	int jcount;
	int l;
	int nn;
	double *p;
	double *p2star;
	double *pbar;
	double *pstar;
	double rcoeff = 1.0;
	double rq;
	double x;
	double *y;
	double y2star;
	double ylo;
	double ystar;
	double z;
	//
	//  Check the input parameters.
	//
	if (reqmin <= 0.0)
	{
		*ifault = 1;
		return;
	}

	if (n < 1)
	{
		*ifault = 1;
		return;
	}

	if (konvge < 1)
	{
		*ifault = 1;
		return;
	}

	p = new double[n*(n + 1)];
	pstar = new double[n];
	p2star = new double[n];
	pbar = new double[n];
	y = new double[n + 1];

	*icount = 0;
	*numres = 0;

	jcount = konvge;
	dn = (double)(n);
	nn = n + 1;
	dnn = (double)(nn);
	del = 1.0;
	rq = reqmin * dn;
	//
	//  Initial or restarted loop.
	//
	for (; ; )
	{
		for (i = 0; i < n; i++)
		{
			p[i + n*n] = start[i];
		}
		y[n] = chi2(start);
		*icount = *icount + 1;

		for (j = 0; j < n; j++)
		{
			x = start[j];
			start[j] = start[j] + step[j] * del;
			for (i = 0; i < n; i++)
			{
				p[i + j*n] = start[i];
			}
			y[j] = chi2(start);
			*icount = *icount + 1;
			start[j] = x;
		}
		//                    
		//  The simplex construction is complete.
		//                    
		//  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
		//  the vertex of the simplex to be replaced.
		//                
		ylo = y[0];
		ilo = 0;

		for (i = 1; i < nn; i++)
		{
			if (y[i] < ylo)
			{
				ylo = y[i];
				ilo = i;
			}
		}
		//
		//  Inner loop.
		//
		for (; ; )
		{
			if (kcount <= *icount)
			{
				break;
			}
			*ynewlo = y[0];
			ihi = 0;

			for (i = 1; i < nn; i++)
			{
				if (*ynewlo < y[i])
				{
					*ynewlo = y[i];
					ihi = i;
				}
			}
			//
			//  Calculate PBAR, the centroid of the simplex vertices
			//  excepting the vertex with Y value YNEWLO.
			//
			for (i = 0; i < n; i++)
			{
				z = 0.0;
				for (j = 0; j < nn; j++)
				{
					z = z + p[i + j*n];
				}
				z = z - p[i + ihi*n];
				pbar[i] = z / dn;
			}
			//
			//  Reflection through the centroid.
			//
			for (i = 0; i < n; i++)
			{
				pstar[i] = pbar[i] + rcoeff * (pbar[i] - p[i + ihi*n]);
			}
			ystar = chi2(pstar);
			*icount = *icount + 1;
			//
			//  Successful reflection, so extension.
			//
			if (ystar < ylo)
			{
				for (i = 0; i < n; i++)
				{
					p2star[i] = pbar[i] + ecoeff * (pstar[i] - pbar[i]);
				}
				y2star = chi2(p2star);
				*icount = *icount + 1;
				//
				//  Check extension.
				//
				if (ystar < y2star)
				{
					for (i = 0; i < n; i++)
					{
						p[i + ihi*n] = pstar[i];
					}
					y[ihi] = ystar;
				}
				//
				//  Retain extension or contraction.
				//
				else
				{
					for (i = 0; i < n; i++)
					{
						p[i + ihi*n] = p2star[i];
					}
					y[ihi] = y2star;
				}
			}
			//
			//  No extension.
			//
			else
			{
				l = 0;
				for (i = 0; i < nn; i++)
				{
					if (ystar < y[i])
					{
						l = l + 1;
					}
				}

				if (1 < l)
				{
					for (i = 0; i < n; i++)
					{
						p[i + ihi*n] = pstar[i];
					}
					y[ihi] = ystar;
				}
				//
				//  Contraction on the Y(IHI) side of the centroid.
				//
				else if (l == 0)
				{
					for (i = 0; i < n; i++)
					{
						p2star[i] = pbar[i] + ccoeff * (p[i + ihi*n] - pbar[i]);
					}
					y2star = chi2(p2star);
					*icount = *icount + 1;
					//
					//  Contract the whole simplex.
					//
					if (y[ihi] < y2star)
					{
						for (j = 0; j < nn; j++)
						{
							for (i = 0; i < n; i++)
							{
								p[i + j*n] = (p[i + j*n] + p[i + ilo*n]) * 0.5;
								xmin[i] = p[i + j*n];
							}
							y[j] = chi2(xmin);
							*icount = *icount + 1;
						}
						ylo = y[0];
						ilo = 0;

						for (i = 1; i < nn; i++)
						{
							if (y[i] < ylo)
							{
								ylo = y[i];
								ilo = i;
							}
						}
						continue;
					}
					//
					//  Retain contraction.
					//
					else
					{
						for (i = 0; i < n; i++)
						{
							p[i + ihi*n] = p2star[i];
						}
						y[ihi] = y2star;
					}
				}
				//
				//  Contraction on the reflection side of the centroid.
				//
				else if (l == 1)
				{
					for (i = 0; i < n; i++)
					{
						p2star[i] = pbar[i] + ccoeff * (pstar[i] - pbar[i]);
					}
					y2star = chi2(p2star);
					*icount = *icount + 1;
					//
					//  Retain reflection?
					//
					if (y2star <= ystar)
					{
						for (i = 0; i < n; i++)
						{
							p[i + ihi*n] = p2star[i];
						}
						y[ihi] = y2star;
					}
					else
					{
						for (i = 0; i < n; i++)
						{
							p[i + ihi*n] = pstar[i];
						}
						y[ihi] = ystar;
					}
				}
			}
			//
			//  Check if YLO improved.
			//
			if (y[ihi] < ylo)
			{
				ylo = y[ihi];
				ilo = ihi;
			}
			jcount = jcount - 1;

			if (0 < jcount)
			{
				continue;
			}
			//
			//  Check to see if minimum reached.
			//
			if (*icount <= kcount)
			{
				jcount = konvge;

				z = 0.0;
				for (i = 0; i < nn; i++)
				{
					z = z + y[i];
				}
				x = z / dnn;

				z = 0.0;
				for (i = 0; i < nn; i++)
				{
					z = z + pow(y[i] - x, 2);
				}

				if (z <= rq)
				{
					break;
				}
			}
		}
		//
		//  Factorial tests to check that YNEWLO is a local minimum.
		//
		for (i = 0; i < n; i++)
		{
			xmin[i] = p[i + ilo*n];
		}
		*ynewlo = y[ilo];

		if (kcount < *icount)
		{
			*ifault = 2;
			break;
		}

		*ifault = 0;

		for (i = 0; i < n; i++)
		{
			del = step[i] * eps;
			xmin[i] = xmin[i] + del;
			z = chi2(xmin);
			*icount = *icount + 1;
			if (z < *ynewlo)
			{
				*ifault = 2;
				break;
			}
			xmin[i] = xmin[i] - del - del;
			z = chi2(xmin);
			*icount = *icount + 1;
			if (z < *ynewlo)
			{
				*ifault = 2;
				break;
			}
			xmin[i] = xmin[i] + del;
		}

		if (*ifault == 0)
		{
			break;
		}
		//
		//  Restart the procedure.
		//
		for (i = 0; i < n; i++)
		{
			start[i] = xmin[i];
		}
		del = eps;
		*numres = *numres + 1;
	}
	delete[] p;
	delete[] pstar;
	delete[] p2star;
	delete[] pbar;
	delete[] y;

	return;
}
//****************************************************************************80

void RKfitter::timestamp(void)

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

	static char time_buffer[TIME_SIZE];
	const struct tm *tm;
	size_t len;
	time_t now;

	now = time(NULL);
	tm = localtime(&now);

	len = strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

	cout << time_buffer << "\n";

	return;
# undef TIME_SIZE
}

