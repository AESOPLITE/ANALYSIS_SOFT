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
	int nLayers;
	vector<double> zLayer;             // the z position of each SSD layer
	vector<char> orientation;          // orientation of each layer 'n' or 'b'
	vector<double> ladderOffsetLeft;
	vector<double> ladderOffsetRight;
	vector<vector<double>> hits;       // a vector of hit coordinate values for each layer
	double stripPitch;

	TkrData();
	void addLyr(char orientation, double z, double offsetLeft, double offsetRight);
	void addHit(int layer, double value);
	void clrHits();
	void print(string s);
	~TkrData();
};

#endif TKRDATA_h