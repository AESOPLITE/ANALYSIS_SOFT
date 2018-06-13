#include "RungeKutta4.h"

RungeKutta4::RungeKutta4(double dz, FieldMap *fM) {  // Version with no scattering
	// dz is the step size
	// fM is the magnetic field map
	this->fM = fM;
	h = dz;
	h2 = h*h;
	nStep = 0;
	aSize = 150;
	xA = new double[aSize];
	yA = new double[aSize];
	zA = new double[aSize];
	sA = new double[aSize];
	double rho = 2.329;            // Density of silicon in g/cm^2
	X0 = (21.82 / rho) * 10.0;
	//	for (int i = 0; i < 100; i++) {
	//		double num = distribution(generator);
	//		cout << " random number = " << num << endl;
	//	}
	this->nScat = 0;
	this->zScat = new double;
	scatAng = new double;
}

RungeKutta4::RungeKutta4(double dz, FieldMap *fM, int nScat, double zScat[]) {
	// dz is the step size
	// fM is the magnetic field map
	this->fM = fM;
	h = dz;
	h2 = h*h;
	nStep = 0;
	aSize = 150;
	xA = new double[aSize];
	yA = new double[aSize];
	zA = new double[aSize];
	sA = new double[aSize];
	double rho = 2.329;            // Density of silicon in g/cm^2
	X0 = (21.82 / rho) * 10.0;
//	for (int i = 0; i < 100; i++) {
//		double num = distribution(generator);
//		cout << " random number = " << num << endl;
//	}
	this->nScat = nScat;
	if (nScat > 0) {
		this->zScat = new double[nScat];
		scatAng = new double[nScat];
		for (int i = 0; i < nScat; i++) this->zScat[i] = zScat[i];
	}
}

// Do a Runge-Kutta 4th-order integration through the field map starting from r0 and p0 and for distance s
// First, a version that can be called when no random scattering is performed (see the constructors)
double *RungeKutta4::Integrate(double Q, double r0[3], double p0[3], double s, double deltaZ) {
	double *rN;
	rN = NULL;
	return Integrate(Q, r0, p0, s, deltaZ, rN);
}

double *RungeKutta4::Integrate(double Q, double r0[3], double p0[3], double s, double deltaZ, double ranNorm[]) {
	// Q is the charge sign
	// r0 is the initial point in mm
	// p0 is the initial momentum in GeV/c
	// s is the distance to propagate (approximate to distance dx)
	// deltaZ is the maximum absolute value change in z.  Set very large to ignore.
	// ranNorm is a set of normally distributed random numbers for the scattering in each plane
	alpha = Q * 2.99792458e-4;
	double *r = new double[6];
	for (int i = 0; i < 3; i++) {
		r[i] = r0[i];
		r[i + 3] = p0[i];
	}
	nStep = (int)(s / h) + 1;
	if (nStep > aSize) {
		delete[] xA;
		delete[] yA;
		delete[] zA;
		delete[] sA;
		aSize = nStep;
		xA = new double[aSize];
		yA = new double[aSize];
		zA = new double[aSize];
		sA = new double[aSize];
		cout << "RungeKutta4::integrate: INCREASING ARRAY SIZE TO " << aSize << endl;
	}
	double sNow = 0.;
	int iScat = 0;
	for (int step = 0; step < nStep; step++) {
		double ri[3] = { r[0], r[1], r[2] };
		double pi[3] = { r[3], r[4], r[5] };
		double *k1 = f(ri, pi);
		double r1[3] = { r[0] + h2 * k1[0], r[1] + h2 * k1[1], r[2] + h2 * k1[2] };
		double p1[3] = { r[3] + h2 * k1[3], r[4] + h2 * k1[4], r[5] + h2 * k1[5] };
		double *k2 = f(r1, p1);
		double r2[3] = { r[0] + h2 * k2[0], r[1] + h2 * k2[1], r[2] + h2 * k2[1] };
		double p2[3] = { r[3] + h2 * k2[3], r[4] + h2 * k2[4], r[5] + h2 * k2[5] };
		double *k3 = f(r2, p2);
		double r3[3] = { r[0] + h * k3[0], r[1] + h * k3[1], r[2] + h * k3[2] };
		double p3[3] = { r[3] + h * k3[3], r[4] + h * k3[4], r[5] + h * k3[5] };
		double *k4 = f(r3, p3);
		for (int i = 0; i < 6; i++) {
			r[i] = r[i] + h * (k1[i] / 6. + k2[i] / 3. + k3[i] / 3. + k4[i] / 6.);
		}
		xA[step] = r[0];
		yA[step] = r[1];
		zA[step] = r[2];
		sNow = sNow + h;
		sA[step] = sNow;
		delete[] k1;
		delete[] k2;
		delete[] k3;
		delete[] k4;
		if (abs(r[2] - r0[2]) >= deltaZ) {
			nStep = step + 1;
			break;
		}
		if (iScat < nScat) {
			if (step > 0) {
				//cout << step << " " << zScat[iScat] << " " << zA[step] << " " << zA[step - 1] << endl;
				if ((zScat[iScat] > zA[step - 1] && zScat[iScat] <= zA[step]) || (zScat[iScat] < zA[step - 1] && zScat[iScat] >= zA[step])) {
					double pmom = sqrt(r[3] * r[3] + r[4] * r[4] + r[5] * r[5]);
					double ct = abs(r[5]) / pmom;
					double radLens = 0.4 / X0 / ct;
					double theta0 = (sqrt(radLens) * 0.0136 / pmom) * (1.0 + 0.038 * log(radLens));
					double thetax = theta0*ranNorm[iScat];
					double thetay = theta0*ranNorm[iScat];
					double t[3];
					for (int i = 0; i < 3; i++) t[i] = r[i + 3] / pmom;
					double u[3] = { t[1] / sqrt(t[0] * t[0] + t[1] * t[1]) , -t[0] / sqrt(t[0] * t[0] + t[1] * t[1]), 0. };
					double v[3] = { -t[2] * u[1] , t[2] * u[0], t[0] * u[1] - t[1] * u[0] };
					double tp[3];
					tp[0] = sin(thetax);
					tp[1] = sin(thetay);
					tp[2] = sqrt(1.0 - tp[0] * tp[0] - tp[1] * tp[1]);
					//cout << "p=" << pmom << " ct=" << ct << " theta0=" << theta0 << " thetax=" << thetax << " thetay=" << thetay << endl;
					//cout << " t=" << t[0] << " " << t[1] << " " << t[2] << endl;
					//cout << " u=" << u[0] << " " << u[1] << " " << u[2] << endl;
					//cout << " v=" << v[0] << " " << v[1] << " " << v[2] << endl;
					//cout << " tp=" << tp[0] << " " << tp[1] << " " << tp[2] << endl;
					double dot = 0;
					for (int i = 0; i < 3; i++) {
						double tRot = u[i] * tp[0] + v[i] * tp[1] + t[i] * tp[2];
						r[3 + i] = tRot*pmom;  // Rotated momentum vector
						dot += t[i] * tRot;
						//cout << "   i=" << i << " p[i]=" << r[3 + i] << " tRot=" << tRot << endl;
					}
					scatAng[iScat] = acos(dot);
					//cout << "Layer " << iScat << "; zA = " << zA[step] << " " << zA[step-1] << "   scattering angle = " << scatAng[iScat] << endl;
					iScat++;
				}
			}
		}
	}
	for (int i = iScat; i < nScat; i++) scatAng[iScat] = -99.;
	return r;
}

int RungeKutta4::getNsteps() {
	return nStep;
}

// Get the position on the trajectory at point i
double *RungeKutta4::getX(int i) {
	double *r = new double[3];
	r[0] = xA[i];
	r[1] = yA[i];
	r[2] = zA[i];
	return r;
}

// Interpolate the position on the trajectory at height z
double *RungeKutta4::getX(double z, bool *flag) {
	*flag = false;
	int sgn;
	if (zA[1] > zA[0]) sgn = 1; else sgn = -1;  // Going up or down?
	double dz = zA[1] - zA[0];
	if (dz == 0.) dz = sgn*h;
	int imax = nStep;
	int i0 = (int)floor((z - zA[0]) / dz);
	if (i0 >= imax - 1) i0 = imax - 2;
	if (i0 < 0) i0 = 0;
	/*
	cout << "zA= ";
	for (int i = 0; i < nStep; i++) {
		cout << zA[i] << " ";
	}
	cout << endl;
	cout << " z=" << z << " dz=" << dz << " zA[0]=" << zA[0] << " imax=" << imax << " sgn=" << sgn << endl;
	*/
	int idx = i0;

	if (sgn < 0) {
		//cout << "getX:" << i0 << " " << zA[i0] << " " << z << endl;
		if (zA[i0] < z) {
			for (int i = i0; i > 0; i--) {
				//cout << "   < i=" << i << " zA=" << zA[i] << endl;
				if (zA[i - 1] > z) {
					idx = i;
					*flag = true;
					break;
				}
			}
		}
		else {
			for (int i = i0 + 1; i < imax; i++) {
				//cout << "   > i=" << i << " zA=" << zA[i] << endl;
				if (zA[i] < z) {
					idx = i;
					*flag = true;
					break;
				}
			}
		}
	}
	if (!flag || z < zA[idx] || z > zA[idx + sgn]) { // well, rats! Let's try a dumb, slow method. . .
		idx = -1;
		for (int i = 0; i < nStep-1; i++) {
			if (z < zA[i + 1] && z > zA[i]) {
				sgn = 1;
				idx = i;
				*flag = true;
				break;
			}
			if (z < zA[i] && z > zA[i + 1]) {
				sgn = -1;
				idx = i+1;
				*flag = true;
				break;
			}
		}
		if (!flag) {
			double *r = new double[3];
			r[0] = 999.;
			r[1] = 999.;
			r[2] = 999.;
			return r;
		}
		//cout << "RungeKutta4:getX: dumb method was successfully used to find the interpolation point." << endl;
	}
	double del = (z - zA[idx]) / (zA[idx + sgn] - zA[idx]);
	double *r = new double[3];
	//cout << "idx=" << idx << endl;
	r[0] = xA[idx] + (xA[idx + sgn] - xA[idx])*del;   // Linear interpolation
	r[1] = yA[idx] + (yA[idx + sgn] - yA[idx])*del;
	r[2] = z;
	//cout << "   z=" << z << " del=" << del << "  zA[idx]=" << zA[idx] << " zA[idx+sgn]=" << zA[idx + sgn] << endl;
	double zint = zA[idx] + (zA[idx + sgn] - zA[idx])*del;
	//cout << "zint=" << zint << endl;

	return r;
}

double RungeKutta4::getS(int i) {
	return sA[i];
}

double *RungeKutta4::f(double x[3], double p[3]) { // Return all the derivatives for Lorentz force
	double *B= fM->GetField(x[0], x[1], x[2]);
	double *d = new double[6];
	double pmag = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
	d[0] = p[0] / pmag; // dx/ds (assuming the electron moves at the speed of light)
	d[1] = p[1] / pmag;
	d[2] = p[2] / pmag;
	d[3] = alpha * (d[1] * B[2] - d[2] * B[1]); // dp/ds
	d[4] = alpha * (d[2] * B[0] - d[0] * B[2]);
	d[5] = alpha * (d[0] * B[1] - d[1] * B[0]);
	delete[] B;
	return d;
}

double RungeKutta4::getScat(int i) {
	if (i < nScat) {
		return scatAng[i];
	}
	else {
		return -99.;
	}
}

RungeKutta4::~RungeKutta4()
{
	//cout << "************** RungeKutta4 cleaning up arrays ************" << endl;
	delete[] xA;
	xA = NULL;
	delete[] yA;
	yA = NULL;
	delete[] zA;
	zA = NULL;
	delete[] sA;
	sA = NULL;
	delete[] scatAng;
	scatAng = NULL;
	delete[] zScat;
	zScat = NULL;
}
