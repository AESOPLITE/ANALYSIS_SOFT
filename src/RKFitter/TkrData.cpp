#include "TkrData.h"

TkrData::TkrData()
{
	nLayers = 0;
	stripPitch = 0.228;
	zLayer.reserve(7);
	orientation.reserve(7);
	ladderOffsetLeft.reserve(7);
	ladderOffsetRight.reserve(7);
	hits.reserve(7);
}

void TkrData::addLyr(char orient, double z, double offsetLeft, double offsetRight) {
	nLayers++;
	zLayer.push_back(z);
	orientation.push_back(orient);
	ladderOffsetLeft.push_back(offsetLeft);
	ladderOffsetRight.push_back(offsetRight);
	vector<double> hitList;
	hits.push_back(hitList);
}

void TkrData::print(string s) {
	cout << "Tracker data " << s << "  Number of layers = " << nLayers << endl;
	for (int lyr = 0; lyr < nLayers; lyr++) {
		cout << "    Layer " << lyr << " orientation is " << orientation[lyr];
		cout << "  z = " << zLayer[lyr] << "  Offsets= " << ladderOffsetLeft[lyr] << " " << ladderOffsetRight[lyr] << endl;
		cout << "          Hits = ";
		for (unsigned int hit = 0; hit < hits[lyr].size(); hit++) {
			cout << hits[lyr].at(hit) << " ";
		}
		cout << endl;
	}
}

void TkrData::addHit(int layer, double value) {
	if (layer < nLayers) {
		hits[layer].push_back(value);
	}
	else {
		cout << "TkrData::addLyr: layer " << layer << " is not yet defined" << endl;
	}
}

void TkrData::clrHits() {
	for (int layer = 0; layer < nLayers; layer++) {
		hits[layer].clear();
	}
}

TkrData::~TkrData()
{
	zLayer.clear();
	orientation.clear();
	ladderOffsetLeft.clear();
	ladderOffsetRight.clear();
	for (int lyr = 0; lyr < nLayers; lyr++) hits.at(lyr).clear();
	hits.clear();
}
