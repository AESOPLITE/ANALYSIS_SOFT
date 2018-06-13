// RKTrackFit.cpp : Defines the entry point for the console application.
//

#define _USE_MATH_DEFINES

#include <cstdio>
#include <chrono>
#include <vector>
#include <cmath>
#include <string>
#include <random>

#include "FieldMap.h"
#include "RungeKutta4.h"
#include "TkrData.h"
#include "RKfitter.h"
#include "Histogram.h"

using namespace std;

// Program to test fitting AESOP-Lite hits to a track model built on Runge-Kutta integration through the field map
// Robert P. Johnson       May 3, 2018
int main()
{
	cout << "RKTrackFit: Entering Main" << endl;
	typedef std::chrono::high_resolution_clock myclock;
	myclock::time_point beginning = myclock::now();

	// obtain a seed from the timer
	myclock::duration d = myclock::now() - beginning;
	unsigned seed = d.count();

	std::mt19937 generator(seed);
	double sigma = 0.228 / sqrt(12.0);
	cout << "     Smearing the measurements by a sigma of " << sigma << endl;
	normal_distribution<double> distribution(0., sigma);

	// Load the B field map. Reading the binary file is far faster than reading the text version (but is machine dependent)
	string fN = "C:\\Users\\Robert\\Documents\\Work\\AESOP\\HTIDS_2014\\Recon\\TrackFit\\RKTrackFit\\fieldmap5mm.bin";
	FieldMap *fM = new FieldMap(fN, "binary", 81);

/*
	// Test building the field map by function call
	int tsize = 81 * 81 * 81;
	float *BX = new float[tsize];
	float *BY = new float[tsize];
	float *BZ = new float[tsize];
	for (int i = 0; i < tsize; i++) {
		float *Bfield = fM->getEntry(i);
		BX[i] = Bfield[0];
		BY[i] = Bfield[1];
		BZ[i] = Bfield[2];
		delete[] Bfield;
	}
	delete fM;
	fM = new FieldMap();
	fM->setOriginStep();
	for (int i = 0; i < tsize; i++) {
		fM->addEntry(BX[i], BY[i], BZ[i]);
	}
	delete[] BX;
	delete[] BY;
	delete[] BZ;
*/

	// When the field map is read from a text file, we can write it back out in binary with the following line of code:
	//fM->writeBinaryFile("C:\\Users\\Robert\\Documents\\Work\\AESOP\\HTIDS_2014\\Recon\\TrackFit\\RKTrackFit\\fieldmap5mm.bin");

	// Get the tracker geometry
	double zLayer[7] = { -14.60, -34.60, -94.60, -154.6, -174.6, -194.6, -214.6 };
	char oLayer[7] = { 'n','b','b','b','n','b','n' };    // Orientation, nonbending or bending
	int nLyrs = 7;
	double ladderOffsetLeft[7] = { 8.889, 8.866, 8.904, 8.871, 8.886, 8.873, 8.884 };
	double ladderOffsetRight[7] = { 98.472, 98.447, 98.478, 98.465, 98.460, 98.466, 98.466 };

	// Create a data list to pass to the fitting routine
	TkrData *Td = new TkrData();
	for (int lyr = 0; lyr < nLyrs; lyr++) {
		Td->addLyr(oLayer[lyr], zLayer[lyr], ladderOffsetLeft[lyr], ladderOffsetRight[lyr]);
	}
	Td->print("testing testing");

	// Generate a track to fit
	double c = 2.99793e8; // Speed of light in m/s
	//double alpha = 1.0e12 / (c * B); // Convert from pt in GeV to curvature in mm
	double xStarti[3] = { 0., 0., 0. };
	double thetaI = 175.0;
	double phiI = 45.0*M_PI / 180.;
	double dcos[3];

	double momentum = 0.05;    // GeV
	double Q = 1.;

	Histogram hop(200, -50., 0.5, "Q/momentum", "1/GeV", "tracks");
	Histogram hopr(200, -10., 0.1, "Q/momentum error", "sigmas", "tracks");
	Histogram hxerr(200, -5., 0.05, "x position error", "mm", "tracks");
	Histogram hyerr(200, -5., 0.05, "y position error", "mm", "tracks");
	Histogram hresn(100, -1., 0.02, "Residual in non-bending plane", "mm", "hits");
	Histogram hresb(100, -1., 0.02, "Residual in bending plane", "mm", "hits");
	Histogram hchi2(100, 0., 25., "Fit chi^2", "chi^2", "tracks");
	Histogram hTchi2(100, 0., 400., "MC truth chi^2", "chi^2", "tracks");
	Histogram hIchi2(100, 0., 400., "Initial-guess chi^2", "chi^2", "tracks");
	Histogram hscat(100, 0., 0.002, "MC scattering angles", "rad", "layers");

	default_random_engine generator2;
	uniform_real_distribution<double> distribution2(0.0, 1.0);
	default_random_engine generator3;
	normal_distribution<double> distribution3(0., 1.0);
	int nIterations = 500;
	int itrPrnt = 0;
	int nScatters = 6;
	double *ranNorm = new double[nScatters+1];
	FILE *fPict;
	double xStart[5];
	for (int iter = 0; iter < nIterations; iter++) {
		// Randomly move around the starting point and direction event to event, but not the momentum
		xStart[0] = xStarti[0] + 20.0*(distribution2(generator2)-0.5);  
		xStart[1] = xStarti[1] + 20.0*(distribution2(generator2)-0.5);
		xStart[2] = xStarti[2];
		double phi = phiI + 2.0*M_PI*distribution2(generator2);
		double theta = (thetaI + (180.0-thetaI)*distribution2(generator2))*M_PI/180.0;
		dcos[0] = cos(phi)*sin(theta);
		dcos[1] = sin(phi)*sin(theta);
		dcos[2] = cos(theta);
		double K = Q / momentum;
		double pStart[3] = { momentum*dcos[0], momentum*dcos[1], momentum*dcos[2] };

		double truth[5];
		truth[0] = xStart[0];
		truth[1] = xStart[1];
		truth[2] = dcos[0];
		truth[3] = dcos[1];
		truth[4] = K;

		// Introduce random Gaussian errors into the initial "guess"
		double guess[5];
		guess[0] = xStart[0] + 0.25*distribution3(generator3);
		guess[1] = xStart[1] + 0.25*distribution3(generator3);
		guess[2] = dcos[0] + 0.004*distribution3(generator3);
		guess[3] = dcos[1] + 0.004*distribution3(generator3);
		guess[4] = K*(1.0 + 0.025*distribution3(generator3));
		if (iter < 10) cout << "Iteration " << iter << " x=" << xStart[0] << " y=" << xStart[1] << " phi=" << phi << " theta=" << theta << endl;

		double dz = 5.0;
		for (int i = 0; i < nScatters; i++) {
			ranNorm[i] = distribution3(generator3);
			//cout << "ran = " << ranNorm[i] << endl;
		};
		RungeKutta4 *rk4;
		double *xEnd;
		if (nScatters > 0) {
			rk4 = new RungeKutta4(dz, fM, nScatters, zLayer);
			xEnd = rk4->Integrate(Q, xStart, pStart, 300., 9999., ranNorm);
		}
		else {
			rk4 = new RungeKutta4(dz, fM);
			xEnd = rk4->Integrate(Q, xStart, pStart, 300., 9999.);
		}
		
		delete[] xEnd;

		// Write out a gnuplot file for plotting the track and the simulated hits
		if (iter == itrPrnt) {
			fPict = fopen("C:\\Users\\Robert\\Documents\\Work\\AESOP\\HTIDS_2014\\Recon\\TrackFit\\MCtrack.gp", "w");
			if (fPict == NULL) {
				cout << "RKTrackFit: could not open the gnuplot output file" << endl;
				return -1;
			}
			fprintf(fPict, "set xlabel 'X'\n");
			fprintf(fPict, "set ylabel 'Y'\n");
			fprintf(fPict, "$runga0 << EOD\n");
			for (int i = 0; i < rk4->getNsteps(); i++) {
				double *xPlt = rk4->getX(i);
				fprintf(fPict, "%10.6f %10.6f %10.6f\n", xPlt[0], xPlt[1], xPlt[2]);
				delete[] xPlt;
			}
			fprintf(fPict, "EOD\n");
			fprintf(fPict, "$hits0 << EOD\n");
		}

		// Generate a "measurement" at each of the 7 layers
		Td->clrHits();
		for (int lyr = 0; lyr < nLyrs; lyr++) {
			double *xMC;
			bool flag;
			xMC = rk4->getX(zLayer[lyr], &flag);
			if (iter == itrPrnt) cout << "    Layer " << lyr << " flag=" << flag << " xMC=" << xMC[0] << " " << xMC[1] << " " << xMC[2] << endl;
			if (!flag) continue;
			double grn = distribution(generator);  // Gaussian smearing of the measurement
			if (oLayer[lyr] == 'n') {
				xMC[0] += grn;    // non-bending view measures x (B-field is in x direction)
				Td->addHit(lyr, xMC[0]);
			}
			else {
				xMC[1] += grn;    // bending view measures y
				Td->addHit(lyr, xMC[1]);
			}
			if (iter == itrPrnt) fprintf(fPict, "%10.6f %10.6f %10.6f\n", xMC[0], xMC[1], xMC[2]);
			delete[] xMC;
		}
		if (iter == itrPrnt) {
			fprintf(fPict, "EOD\n");
			fprintf(fPict, "splot $runga0 u 1 : 2 : 3 with lines lw 3, $hits0 u 1 : 2 : 3 with points pt 6 ps 2\n");
			fclose(fPict);
		}
		for (int i = 0; i < nScatters; i++) {
			double ang = rk4->getScat(i);
			//cout << " Scat ang=" << ang << endl;
			hscat.entry(ang);
		}
		delete rk4;

		if (iter < 1) Td->print("test");  // Prints the data in a nice format

		RKfitter *rkf = new RKfitter(iter<0, xStart[2], fM, Td);
		double initChi2 = rkf->chi2(guess);
		double truthChi2 = rkf->chi2(truth);
		vector<int> hits = { 0, 0, 0, 0, 0, 0, 0 };
		int errCode = rkf->fitIt(false, guess, hits);
		if (iter == itrPrnt) rkf->print("test fit");
		//cout << "Initial chi^2 = " << initChi2 << "  Truth chi^2 target=" << truthChi2 << " Error=" << errCode << "  Fit chi^2=" << rkf->chiSqr() << endl;

		double a[5];
		rkf->tkrParams(a);
		double e[5];
		rkf->errors(e);
		hop.entry(a[4]);
		hopr.entry((a[4]-K)/e[4]);
		hxerr.entry((a[0] - xStart[0]) / e[0]);
		hyerr.entry((a[1] - xStart[1]) / e[1]);
		hchi2.entry(rkf->chiSqr());
		hIchi2.entry(initChi2);
		hTchi2.entry(truthChi2);
		for (int lyr = 0; lyr < Td->nLayers; lyr++) {
			double r[3];
			rkf->getIntercept(lyr, r);
			if (Td->orientation[lyr] == 'n') {
				double residual = r[0] - Td->hits[lyr].at(0);
				hresn.entry(residual);
			} else {
				double residual = r[1] - Td->hits[lyr].at(0);
				hresb.entry(residual);
			}
		}
		cout << "Event " << iter << " chi^=" << rkf->chiSqr() << " Error=" << errCode << " Track=" << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << " " << a[4] << endl;
		delete rkf;
	}

	cout << "output histograms" << endl;
	string path = "C:\\Users\\Robert\\Documents\\Work\\AESOP\\HTIDS_2014\\Recon\\TrackFit\\";
	FILE *fp;
	string fileName = path + "inverseP.gp";
	fp = fopen(fileName.c_str(), "w");
	if (fp != NULL) {
		hop.plot(fp, true, " ", " ");
		fclose(fp);
	}
	fileName = path + "nonBendResid.gp";
	fp = fopen(fileName.c_str(), "w");
	if (fp != NULL) {
		hresn.plot(fp, true, " ", " ");
		fclose(fp);
	}
	fileName = path + "bendResid.gp";
	fp = fopen(fileName.c_str(), "w");
	if (fp != NULL) {
		hresb.plot(fp, true, " ", " ");
		fclose(fp);
	}
	fileName = path + "chi2.gp";
	fp = fopen(fileName.c_str(), "w");
	if (fp != NULL) {
		hchi2.plot(fp, true, " ", " ");
		fclose(fp);
	}
	fileName = path + "initChi2.gp";
	fp = fopen(fileName.c_str(), "w");
	if (fp != NULL) {
		hIchi2.plot(fp, true, " ", " ");
		fclose(fp);
	}
	fileName = path + "truthChi2.gp";
	fp = fopen(fileName.c_str(), "w");
	if (fp != NULL) {
		hTchi2.plot(fp, true, " ", " ");
		fclose(fp);
	}
	fileName = path + "momError.gp";
	fp = fopen(fileName.c_str(), "w");
	if (fp != NULL) {
		hopr.plot(fp, true, " ", " ");
		fclose(fp);
	}
	fileName = path + "xError.gp";
	fp = fopen(fileName.c_str(), "w");
	if (fp != NULL) {
		hxerr.plot(fp, true, " ", " ");
		fclose(fp);
	}
	fileName = path + "yError.gp";
	fp = fopen(fileName.c_str(), "w");
	if (fp != NULL) {
		hyerr.plot(fp, true, " ", " ");
		fclose(fp);
	}
	fileName = path + "scatAng.gp";
	fp = fopen(fileName.c_str(), "w");
	if (fp != NULL) {
		hscat.plot(fp, true, " ", " ");
		fclose(fp);
	}


	delete fM;
	delete Td;
	cout << "RKTrackFit:  ....Done....   Normal return code." << endl;

    return 0;
}

