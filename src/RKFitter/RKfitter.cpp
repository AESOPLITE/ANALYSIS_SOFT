#include "RKfitter.h"

// Calculate the chi^2 of the fit for a given track. This is what we try to minimize. 
double RKfitter::chi2_2(double S, double C) {
	if (algorithm != 0) icount++;
	double b[5];
	b[0] = a[0];
	b[1] = a[1];
	b[2] = a[2];
	b[3] = S;
	b[4] = C;
	if (multScat) return chi2m(b);
	return chi2nm(b);
}

double RKfitter::chi2(double a[]) {
	if (algorithm != 0) icount++;
	//if (verbose) cout << "Entering RKfitter::chi2 with multScat=" << multScat << " bobyqa=" << bobyqa << " count=" << icount << endl;
	if (multScat) return chi2m(a);
	return chi2nm(a);
}

double RKfitter::chi2nm(double a[]) {  // Version without multiple scattering
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
	//if (verbose) {
	//	cout << "RKfitter::chi2nm: r0=" << r0[0] << " " << r0[1] << " " << r0[2];
	//	cout << "   p0=" << p0[0] << " " << p0[1] << " " << p0[2] << "  s=" << s << endl;
	//}
	double *xEnd = rk4->Integrate(Q, r0, p0, s, abs(Delta_z));
	if (verbose) cout << "completed integration at xEnd=" << xEnd[0] << " " << xEnd[1] << " " << xEnd[2] << endl;
	delete[] xEnd;

	// Loop over layers and add up the chi^2
	double result = 0.;
	for (int lyr = 0; lyr < tD->nLayers; lyr++) {
		//if (verbose) cout << "lyr=" << lyr << " # hits=" << tD->hits[lyr].size() << " orient=" << tD->orientation[lyr] << endl;
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
			}
			else {
				xIntercept[lyr] = rInterp[0];
				yIntercept[lyr] = rInterp[1];
			}
			double xMeas = tD->hits[lyr].at(hits[lyr]);
			double incr = pow((xInterp - xMeas) / sigma, 2);
			result += incr;
			delete[] rInterp;
			//if (verbose) cout << "   Layer " << lyr << " flag=" << flag << " xInterp=" << xInterp << " xMeas=" << xMeas << " chi2_inc=" << incr << endl;
		}
		else {
			double *rInterp = rk4->getX(tD->zLayer[lyr], &flag);
			double yInterp = rInterp[1];
			if (!flag) {
				yInterp = 220.;
				xIntercept[lyr] = -9999.;
				yIntercept[lyr] = -9999.;
			}
			else {
				xIntercept[lyr] = rInterp[0];
				yIntercept[lyr] = rInterp[1];
			}
			double yMeas = tD->hits[lyr].at(hits[lyr]);
			double incr = pow((yInterp - yMeas) / sigma, 2);
			result += incr;
			delete[] rInterp;
			//if (verbose) cout << "   Layer " << lyr << " flag=" << flag << " yInterp=" << yInterp << " yMeas=" << yMeas << " chi2_inc=" << incr << endl;
		}
	}
	if (verbose) {
		cout << "RKfitter::chi2nm: with sigma = " << sigma << "  a=";
		for (int i = 0; i < 5; i++) cout << a[i] << " ";
		cout << "  chi^2=" << result << endl;
	}
	return result;
}

// Calculate the chi^2 of the fit for a given track. This is what we try to minimize. Multiple scattering version.
double RKfitter::chi2m(double a[]) {
    //if (verbose) cout << "RKfitter::chi2m: entering with a=" << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << " " << a[4] << endl;
	// estimate how far to integrate, to cover all layers

	double arg = 1.0 - a[2] * a[2] - a[3] * a[3];
	if (arg <= 0.) {
		return 9.9e9;
	}
	double Delta_z = z0 - (tD->zLayer[tD->nLayers - 1]) + 20.;
	double ctz = -sqrt(arg);
	double s = -Delta_z / ctz + 150.0;
	if (s > 500. || s < 0.) s = 500.;
	//cout << "RKfitter::chi2m: z0=" << z0 << " Delta_z=" << Delta_z << " and s=" << s << endl;

	// Set up the input parameters needed by the integrator
	double Q;
	if (a[4] < 0.) Q = -1.0; else Q = 1.0;
	double p = Q / a[4];
	double r0[3] = { a[0], a[1], z0 };
	double p0[3] = { p*a[2], p*a[3], p*ctz };

	// Integrate all the way through the instrument
	//if (verbose) {
	//	cout << "RKfitter::chi2m: r0=" << r0[0] << " " << r0[1] << " " << r0[2];
	//	cout << "   p0=" << p0[0] << " " << p0[1] << " " << p0[2] << "  s=" << s << endl;
	//}
	double *xEnd= rk4->Integrate(Q, r0, p0, s, abs(Delta_z));
	//if (verbose) cout << "completed integration at xEnd=" << xEnd[0] << " " << xEnd[1] << " " << xEnd[2]  << endl;
	delete[] xEnd;

	// Calculate the covariance matrix of the measurements, including multiple scattering

	double *zLyr = new double[tD->nLayers];
	double *theta0 = new double[tD->nLayers];
	double *residual = new double[tD->nLayers];
	for (int i = 0; i < tD->nLayers; i++) {
		//if (verbose) cout << "RKfitter::chi2m: lyr=" << i << " # hits=" << tD->hits[i].size() << " orient=" << tD->orientation[i] << endl;
		bool flag;

		double *rInterp = rk4->getX(tD->zLayer[i], &flag);  // Find the intersection with each silicon layer

		Cx[i][i] = sigma*sigma;  // Including here the measurement uncertainties in the diagonal elements
		double *pInt = rk4->getP();
		double pmom = sqrt(pInt[0] * pInt[0] + pInt[1] * pInt[1] + pInt[2] * pInt[2]);
		double ct = abs(pInt[2]) / pmom;
		double radLens = SiThickness / X0 / ct;
		theta0[i] = (sqrt(radLens) * 0.0136 / pmom) * (1.0 + 0.038 * log(radLens));
		if (!flag) {
			xIntercept[i] = -9999.;
			yIntercept[i] = -9999.;
		}
		else {
			xIntercept[i] = rInterp[0];
			yIntercept[i] = rInterp[1];
		}
		zLyr[i] = rInterp[2];
		//if (verbose) cout << " layer " << i << " p=" << pmom << " ct=" << ct << " theta0=" << theta0[i] << endl;
		// Calculate the diagonal elements
		for (int k = 0; k < i; k++) { // include the effect of scattering in all higher layers
			double leverArm = abs(zLyr[i] - zLyr[k]);
			double dperp = leverArm * theta0[k];  // lever arm times scattering angle
			Cx[i][i] += dperp*dperp;  // add all the contributions in quadrature
			//if (verbose) cout << "     from layer " << k << " lever=" << leverArm << "  theta=" << theta0[k] << "  dperp=" << dperp << " sigma=" << sqrt(Cx[i][i]) << endl;
		}

		// Calculate the off-diagonal elements
		for (int j = 0; j < i; j++) {
			Cx[i][j] = 0;
			for (int k = 0; k < j; k++) {
				Cx[i][j] += abs(zLyr[i] - zLyr[k])*abs(zLyr[j] - zLyr[k]) * theta0[k] * theta0[k];
			}
			Cx[j][i] = Cx[i][j];
		}

		// Calculate the residual
		if (hits[i] >= 0) {
			if (tD->orientation[i] == 'n') {
				double xInterp = rInterp[0];
				if (!flag) {
					xInterp = 220.;   // force a large chi^2 if the track doesn't even intercept the layer
				}
				double xMeas = tD->hits[i].at(hits[i]);
				residual[i] = (xInterp - xMeas);
			}
			else {
				double yInterp = rInterp[1];
				if (!flag) {
					yInterp = 220.;
				}
				double yMeas = tD->hits[i].at(hits[i]);
				residual[i] = (yInterp - yMeas);
			}
			//if (verbose) cout << "RKfitter::chi2m: layer " << i << " residual=" << residual[i] << endl;
		}
	}
	delete[] zLyr;
	delete[] theta0;

	// Invert the covariance matrix
	/*
	double **T; 
	if (verbose) {
		cout << "RKfitter::chi2m: covariance matrix:" << endl;
		T = new double*[tD->nLayers];
		for (int i = 0; i < tD->nLayers; i++) {
			T[i] = new double[tD->nLayers];
			for (int j = 0; j < tD->nLayers; j++) {
				T[i][j] = Cx[i][j];
				cout << "  " << Cx[i][j];
			}
			cout << endl;
		}
	}
	*/
	int err = invert(Cx, tD->nLayers);
	if (err != 0) cout << "RKfitter::chi2m: singular covariance matrix encountered" << endl;

	// Check the matrix inversion
	/*
	if (verbose) {
		double **U = new double*[tD->nLayers];
		for (int i = 0; i < tD->nLayers; i++) {
			U[i] = new double[tD->nLayers];
			for (int j = 0; j < tD->nLayers; j++) {
				U[i][j] = 0.;
				for (int k = 0; k < tD->nLayers; k++) {
					U[i][j] += Cx[i][k] * T[k][j];
				}
			}
		}
		cout << "  RKfitter::chi2m: Testing the matrix inversion. Unit Matrix:" << endl;
		for (int i = 0; i < tD->nLayers; i++) {
			cout << "         ";
			for (int j = 0; j < tD->nLayers; j++) {
				cout << U[i][j] << "   ";
			}
			cout << endl;
			delete[] T[i];
			delete[] U[i];
		}
	}
	*/

	// Loop over layers and add up the chi^2
	double result = 0.;
	for (int i = 0; i < tD->nLayers; i++) {
		if (hits[i] < 0) continue;
		for (int j = 0; j < tD->nLayers; j++) {
			if (hits[j] < 0) continue;
			result += residual[i] * Cx[i][j] * residual[j];
			if (isnan(result)) {
				result = 9.9e22;
				break;
			}
		}
	}
	delete[] residual;

	if (verbose) {
		cout << "RKfitter::chi2: with sigma = " << sigma << " mult scat=" << multScat << "  a=";
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

RKfitter::RKfitter(bool verbose, double z0,  FieldMap *fM, TkrData *tD, bool multScat, double stepSize, int alg) {
	// verbose - set true to get lots of printout
	// z0 - starting point in z for the integration
	// FieldMap - magnetic field map
	// TkrData - data to be fit
	// multScat - set true to include multiple scattering in the chi^2 calculation
	// stepSize - Runge-Kutta integration step size (normally set to 5mm)
	// alg - minimization algorithm
	//       0 = Nelder-Mead simplex algorithm (no derivatives)
	//       1 = dlib bobyqa algorithm (doesn't work well; not sure why)
	//       2 = dlib gradient algorithm (numerical derivatives); goes to Nelder-Mean if chi^2 is bad
	//       Algorithm 2 is faster than 0 but only slightly better in resolution
	if (verbose) cout << "Entering the RKfitter constructor with z0=" << z0 << endl;
	this->verbose = verbose;
	this->z0 = z0;
	this->fM = fM;
	this->tD = tD;
	this->multScat = multScat;    
	this->algorithm = alg;
	hits = new int[tD->nLayers];
	for (int i = 0; i < tD->nLayers; i++) hits[i] = 0;
	xIntercept = new double[tD->nLayers];
	yIntercept = new double[tD->nLayers];
	zIntercept = new double[tD->nLayers];
	for (int lyr = 0; lyr < tD->nLayers; lyr++) {
		zIntercept[lyr] = tD->zLayer[lyr];
	}
	maxCalls = 1000;   // Maximum function calls allowed in the minimization search
	reqmin = 0.0001;    // convergence check parameter
	step = new double[5];
	step[0] = 1.;    // initial step for x position in the minimization search
	step[1] = 1.;    // initial step for y position
	step[2] = 0.04;  // initial step size for x direction cosine
	step[3] = 0.04;  // initial step size for y direction cosine
	step[4] = 20.0;    // initial step size for 1/p as a percentage
	this->stepSize = stepSize;   // Runge-Kutta integration step size (comparable to the B field map precision)
	double rho = 2.329;           // Density of silicon in g/cm^2
	SiThickness = 0.4;
	X0 = (21.82 / rho) * 10.0;

	rk4 = new RungeKutta4(stepSize, fM);
	sigma = tD->stripPitch / sqrt(12.0);  // Assumed resolution of the AESOP tracker
	a = new double[5];
	C = new double*[5];
	for (int i = 0; i < 5; i++) {
		C[i] = new double[5];
	}
	Cx = new double*[tD->nLayers];
	if (multScat) {
		for (int i = 0; i < tD->nLayers; i++) {
			Cx[i] = new double[tD->nLayers];
		}
	}
}

// Execute the brute-force track fit using a canned minimization method
int RKfitter::fitIt(bool genStartGuess, double guess[5], std::vector<int> hitSelection) {
	// hitSelection specifies which hit to use from each layer of the TkrData structure
	// guess[5] is the externally supplied starting guess for the track
	// if genStartGuess is true, then the program will generate an initial guess from a linear fit

	if (verbose) cout << "RKfitter::fitIt: guess=" << guess[0] << " " << guess[1] << " " << guess[2] << " " << guess[3] << " " << guess[4] << endl;
	if (hitSelection.size() != tD->nLayers) {
		cout << "RKfitter:fitIt, wrong number of hits specified, need " << tD->nLayers << endl;
		return -1;
	}
	for (int lyr = 0; lyr < tD->nLayers; lyr++) hits[lyr] = hitSelection.at(lyr);
	double *temp = new double[5];

	// estimate some reasonable values for the initial step in each of the parameters
	// and copy the initial guess, since the minimization routine may overwrite it
	double initStep[5];
	for (int i = 0; i < 5; i++) {
		temp[i] = guess[i];
		initStep[i] = step[i];
	}
	initStep[4] = initStep[4] * guess[4] / 100.;
	//cout << "RKfitter::fitIt: initial step sizes= ";
	//for (int i = 0; i < 5; i++) cout << initStep[i] << " ";
	//cout << endl;

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
			if (tD->hits[lyr].size() < 1) continue;
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
		if (linearFit(zn, x, Nn, 0, &a, &b, &c, &xsqn) == 0) {  // Line fit
			if (linearFit(zb, y, Nb, 1, &c, &d, &e, &xsqb) == 0) {  // Parabola fit
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

	// Call the minimization routine
#ifdef DLIB
	column_vector start = { temp[0],temp[1],temp[2],temp[3],temp[4] };
	column_vector lower_bound = { -120., -120., -1., -1., -1.e100 };
	column_vector upper_bound = { 120., 120., 1., 1., 1.e100 };
#endif
	icount = 0;
	numres = 0;
	ifault = 0;
	switch(algorithm) {
	case 0:
		nelmin(nVar, temp, a, &ynewlo, reqmin, initStep, Konvge, maxCalls, &icount, &numres, &ifault);
		//cout << "RKfitter::fitIt: from nelmin, newlo=" << ynewlo << " count=" << icount << " numres=" << numres << " error=" << ifault << endl;
		//if (ifault != 0) cout << "ifault=" << ifault << " returned by nelmin" << endl;
		break;
	case 1:
		if (verbose) cout << "Calling bobqa with parameters = " << temp[0] << " " << temp[1] << " " << temp[2] << " " << temp[3] << " " << temp[4] << endl;
#ifdef DLIB
		try {
			dlib::find_min_bobyqa(&chi2b,
				start,
				9,    // number of interpolation points
				lower_bound,
				upper_bound,
				0.45,    // initial trust region radius
				1.e-4,  // stopping trust region radius
				maxCalls    // max number of objective function evaluations
				);
		} catch (error e)	{
			cout << "RKfitter: An exception occurred in find_min_bobyqa: " << e.info << endl;
			ifault++;
		}
		if (verbose) cout << "bobyqa solution: \n" << start << "for " << icount << " function calls." << endl;
		for (int i = 0; i < nVar; i++) {
			a[i] = start(i);
			//cout << "i=" << i << " a=" << a[i] << endl;
		}
#else
		cout << "RKfitter: dlib must be compiled and linked in order to use algorithm 1" << endl;
#endif
		ynewlo = chi2(a);
/*		try {
			auto result = dlib::find_min_global(
				&chi2bend,
				{ max(a[3] - 0.0025,-1.), a[4] * 0.99 },
				{ min(a[3] + 0.0025,1.), a[4] * 1.01 },
				max_function_calls(100)
			);
			double tmp[5];
			tmp[0] = a[0];
			tmp[1] = a[1];
			tmp[2] = a[2];
			tmp[3] = result.x(0);
			tmp[4] = result.x(1);
			double c2 = chi2(tmp);
			if (c2 < ynewlo) {
				ynewlo = c2;
				a[3] = tmp[3];
				a[4] = tmp[4];
			}
			if (verbose) cout << "global solution: \n" << result.x << "for " << icount << " function calls. chi2=" << c2 << endl;
		}
		catch (error e) {
			cout << "RKfitter: An exception occurred in find_min_global: " << e.info << endl;
			ifault++;
		} */
		break;
	case 2:
		bool erroc = false;
#ifdef DLIB
		try {
			find_min_using_approximate_derivatives(bfgs_search_strategy(),
				objective_delta_stop_strategy(1e-2),
				&chi2b, start, -1);
		}
		catch (error e) {
			cout << "RKfitter: An exception occurred in find_min_using_approximate_derivatives: " << e.info << endl;
			ifault++;
			erroc = true;
		}
		for (int i = 0; i < nVar; i++) {
			a[i] = start(i);
			//cout << "i=" << i << " a=" << a[i] << endl;
		}
#else
		cout << "dlib must be compiled and linked in order to use algorithm 2 " << endl;
#endif
		ynewlo = chi2(a);
		if (ynewlo > 25. || erroc) {
			//cout << "RKfitter::fitIt: chi2=" << ynewlo << ", switching to Nelder-Mead method." << endl;
			nelmin(nVar, temp, a, &ynewlo, reqmin, initStep, Konvge, maxCalls, &icount, &numres, &ifault);
		}
		break;
	}

	// Propagation of errors:
	double h[5] = { 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001*abs(a[4]) };
	hessian(h);
	//double T[5][5];
	int err = 0;
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			C[i][j] = 0.5*H[i][j];
			if (isnan(C[i][j])) {
				err = 2;
				break;
			}
			//if (i != j) C[i][j] = 0.;
			//T[i][j] = C[i][j];
		}
		if (err != 0) break;
	}

	if (err == 0) err = invert(C,5);  // Invert the hessian matrix (times 0.5) to get the covariance matrix
	if (err != 0) cout << "RKfitter::fitIt: bad or singular covariance matrix encountered. Error=" << err << endl;

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
			if (hits[lyr] < 1) continue;
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
	return ifault;
}

// Return the position at each silicon plane interpolated from the Runge-Kutta integration
void RKfitter::getIntercept(int Layer, double r[3]) {
	r[0] = xIntercept[Layer];
	r[1] = yIntercept[Layer];
	r[2] = zIntercept[Layer];
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
int RKfitter::invert(double **a, int n) { // code slightly modified from Numerical Recipes in C "gaussj"
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
		if (a[icol][icol] == 0.0) {
			cout << "gaussj: singular matrix encountered; aborting the matrix inversion" << endl;
			return -1;
		}
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
	return 0;
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
	if (multScat) for (int i = 0; i < tD->nLayers; i++) delete[] Cx[i];
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

