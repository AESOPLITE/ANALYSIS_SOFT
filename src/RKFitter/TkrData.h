#ifndef TKRDATA_h
#define TKRDATA_h

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

// Summary of the tracker data to be used by the Runge-Kutta track fit
// Robert P. Johnson       May 3, 2018
class TkrData
{
public:
  int eventNumber;
	int nLayers;
	std::vector<double> zLayer;             // the z position of each SSD layer
	std::vector<char> orientation;          // orientation of each layer 'n' or 'b'
	std::vector<double> ladderOffsetLeft;
	std::vector<double> ladderOffsetRight;
	std::vector<std::vector<double>> hits;       // a vector of hit coordinate values for each layer
	std::vector<std::vector<double>> xHitMC;
	std::vector<std::vector<double>> yHitMC;
	std::vector<std::vector<double>> xHitPR;
	std::vector<std::vector<double>> yHitPR;
	double stripPitch;

	TkrData(int eventNumber);
	void addLyr(char orientation, double z, double offsetLeft, double offsetRight);
	void addHit(int layer, double value);
	void addHit(int layer, double value, double xMC, double yMC, double xPR, double yPR);
	void clrHits();
	void print(string s);
	~TkrData();
};

#endif
