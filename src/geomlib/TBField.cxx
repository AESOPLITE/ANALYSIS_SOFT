//*************************************************************************
//* =====================
//*  TBField Class
//* =====================
//*
//* (Description)
//*   A singleton class to hold information of detector system
//*   used in Kalman filter classes.
//* (Requires)
//*     TObject
//* (Provides)
//*     class TBField
//* (Update Recored)
//*   2005/02/23  A.Yamaguchi  	Original version.
//*   2005/08/14  K.Fujii        Removed CalcTable(), GetMeasLayerTable(),
//*                              GetPhiTable(), and GetDir() and added
//*                              Transport() to do their functions.
//*   2010/04/06  K.Fujii        Modified Transport() to allow a 1-dim hit,
//*                              for which pivot is at the expected hit.
//*
//*************************************************************************

#include <memory>      // from STL
#include <cmath>
#include <iostream>    // from STL

#include <TFile.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH3D.h>
#include <map> 
#include "TBField.h"   // from KalTrackLib

ClassImp(TBField)

const double OffsetMag = 106.3409;		//offset of center of magnet to bottom of T3 in mm

Bool_t   TBField::fUseUniformBfield = kFALSE;
//Bool_t   TBField::fUseUniformBfield = kTRUE;
Double_t TBField::fFieldCoeff       = 1.;


///////////////////////////////////////////////////////////////////////
/////////////////////////V2////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
/*
TH3D* TBField::fieldmapX = new TH3D("field map Bx" , "Bx [G]; x[mm]; y[mm]; z[mm]", 81 , -200, 200, 81, -200, 200, 81, -200, 200);
TH3D* TBField::fieldmapY = new TH3D("field map By" , "By [G]; x[mm]; y[mm]; z[mm]", 81 , -200, 200, 81, -200, 200, 81, -200, 200);
TH3D* TBField::fieldmapZ = new TH3D("field map Bz" , "Bz [G]; x[mm]; y[mm]; z[mm]", 81, -200, 200, 81, -200, 200, 81, -200, 200);

using namespace std;

//_________________________________________________________________________
//  ----------------------------------
//   Ctors and Dtor
//  ----------------------------------

TBField::TBField()
{
}

TBField::~TBField()
{
}

//_________________________________________________________________________
//  ----------------------------------
//   Utility Methods
//  ----------------------------------
//_________________________________________________________________________



//function to be called once to set the vector field map
bool TBField::SetMagField() {
	
	//first write field map into root filecd

	 TFile*file=new TFile("fieldmap.root","RECREATE");
	
//Fill fieldmap.root following FLUKA coordinates, do coordinate transform later	
     TNtuple *T= new TNtuple("fieldmap","fieldmap","x:y:z:bx:by:bz");
     T->ReadFile("/home/smechbal/ANALYSIS_SOFT/fieldmap.txt");
	 T->Write();	
	file->Close();

	//now open file
    TFile*file_map=new TFile("fieldmap.root","READ");
	TNtuple*ntuple=(TNtuple*)file_map->Get("fieldmap");

	// create 3D histograms	

    int nentries = ntuple->GetEntries();

	Float_t x;
	Float_t y;
	Float_t z;
	Float_t bx;
	Float_t by;
	Float_t bz;
	ntuple->SetBranchAddress("x",&x);
	ntuple->SetBranchAddress("y",&y);
	ntuple->SetBranchAddress("z",&z);
	ntuple->SetBranchAddress("bx",&bx);
	ntuple->SetBranchAddress("by",&by);
	ntuple->SetBranchAddress("bz",&bz);

	for(int ientry=0; ientry < nentries ;ientry++) {
		ntuple->GetEntry(ientry);

			    fieldmapX->SetBinContent(
				fieldmapX->GetXaxis()->FindBin(x),
				fieldmapX->GetYaxis()->FindBin(y),
				fieldmapX->GetZaxis()->FindBin(z),
				bx);
			    fieldmapY->SetBinContent(
				fieldmapY->GetXaxis()->FindBin(x),
				fieldmapY->GetYaxis()->FindBin(y),
				fieldmapY->GetZaxis()->FindBin(z),
				by);
			    fieldmapZ->SetBinContent(
				fieldmapZ->GetXaxis()->FindBin(x),
				fieldmapZ->GetYaxis()->FindBin(y),
				fieldmapZ->GetZaxis()->FindBin(z),
				bz);

	}
	file_map->Close();
	return true;
}


void TBField::Get(const double&x, const double&y, const double&z, double& BxFLUKA, double& ByFLUKA, double& BzFLUKA) {

//interpolate and convert from [G] to [T]
	//cout << "Interpolate at x = " << x << "  y = " << y << " z = " << z << endl;
	BxFLUKA = (fieldmapX->Interpolate(x,y,z))/10000;
	ByFLUKA = (fieldmapY->Interpolate(x,y,z))/10000;
	BzFLUKA = (fieldmapZ->Interpolate(x,y,z))/10000;
	int binx = fieldmapX->GetXaxis()->FindBin(x);
	int biny = fieldmapX->GetYaxis()->FindBin(y);
	int binz = fieldmapX->GetZaxis()->FindBin(z);	
	if(binx==0 || biny==0 || binz==0) {
	//	cout << "Cannot interpolate at x = " << x << "  y = " << y << " z = " << z << endl;
	//	cout << " bin x = " << binx << " biny = " << biny << " binz = " << binz << endl;
	}

}
*/


/////////////////////////////////////////////////////////
//////////////////////////V3/////////////////////////////
/////////////////////////////////////////////////////////

TH3D* TBField::fieldmapX = new TH3D("field map Bx" , "Bx [T]; x[mm]; y[mm]; z[mm]", 402 , -201, 201, 402, -201, 201, 402, -201, 201);
TH3D* TBField::fieldmapY = new TH3D("field map By" , "By [T]; x[mm]; y[mm]; z[mm]", 402 , -201, 201, 402, -201, 201, 402, -201, 201);
TH3D* TBField::fieldmapZ = new TH3D("field map Bz" , "Bz [T]; x[mm]; y[mm]; z[mm]", 402, -201, 201, 402, -201, 201, 402, -201, 201);
using namespace std;

//_________________________________________________________________________
//  ----------------------------------
//   Ctors and Dtor
//  ----------------------------------

TBField::TBField()
{
}

TBField::~TBField()
{
}

//_________________________________________________________________________
//  ----------------------------------
//   Utility Methods
//  ----------------------------------
//_________________________________________________________________________



//function to be called once to set the vector field map
bool TBField::SetMagField() {
	

    TFile*file_map=new TFile("/data/smechbal/Fluka/NonUniB/V3/fieldmap1mm.root","READ");
	TNtuple*ntuple=(TNtuple*)file_map->Get("Field");

	// create 3D histograms	

    int nentries = ntuple->GetEntries();

	Float_t x;
	Float_t y;
	Float_t z;
	Float_t bx;
	Float_t by;
	Float_t bz;
	ntuple->SetBranchAddress("x",&x);
	ntuple->SetBranchAddress("y",&y);
	ntuple->SetBranchAddress("z",&z);
	ntuple->SetBranchAddress("bx",&bx);
	ntuple->SetBranchAddress("by",&by);
	ntuple->SetBranchAddress("bz",&bz);

	for(int ientry=0; ientry < nentries ;ientry++) {
		ntuple->GetEntry(ientry);

			    fieldmapX->SetBinContent(
				fieldmapX->GetXaxis()->FindBin(x),
				fieldmapX->GetYaxis()->FindBin(y),
				fieldmapX->GetZaxis()->FindBin(z),
				bx);
			    fieldmapY->SetBinContent(
				fieldmapY->GetXaxis()->FindBin(x),
				fieldmapY->GetYaxis()->FindBin(y),
				fieldmapY->GetZaxis()->FindBin(z),
				by);
			    fieldmapZ->SetBinContent(
				fieldmapZ->GetXaxis()->FindBin(x),
				fieldmapZ->GetYaxis()->FindBin(y),
				fieldmapZ->GetZaxis()->FindBin(z),
				bz);

	}
	file_map->Close();
	return true;
}


void TBField::Get(const double&x, const double&y, const double&z, double& BxFLUKA, double& ByFLUKA, double& BzFLUKA) {

//interpolate
	//cout << "Interpolate at x = " << x << "  y = " << y << " z = " << z << endl;
	//BxFLUKA = (fieldmapX->Interpolate(x,y,z));
	//ByFLUKA = (fieldmapY->Interpolate(x,y,z));
	//BzFLUKA = (fieldmapZ->Interpolate(x,y,z));
	int binx = fieldmapX->GetXaxis()->FindBin(x);
	int biny = fieldmapX->GetYaxis()->FindBin(y);
	int binz = fieldmapX->GetZaxis()->FindBin(z);	
	//cout << " bin x = " << binx << " biny = " << biny << " binz = " << binz << endl;
	//cout << "Interpolate at x = " << x << "  y = " << y << " z = " << z << endl;
	//Interpolate field map same way as in magfld.f
	BxFLUKA = fieldmapX->GetBinContent(binx,biny,binz);
	ByFLUKA = fieldmapY->GetBinContent(binx,biny,binz);
	BzFLUKA = fieldmapZ->GetBinContent(binx,biny,binz);
	if(binx==0 || biny==0 || binz==0) {
	//	cout << "Cannot interpolate at x = " << x << "  y = " << y << " z = " << z << endl;
		//cout << " bin x = " << binx << " biny = " << biny << " binz = " << binz << endl;
	}

}

TVector3 TBField::Get(const TVector3& v) {
  //cout << "In Get function KALTEST coord v.x=" << v.X() << "  , v.y=" << v.Y() << " , v.z=" << v.Z() << endl;
	
//convert coordinates from KalTest back to FLUKA, shift z by 10.63cm such that the center of magnet is at z = 0 mm (in FLUKA coord)

  double x = v.z();
  double y = v.x();
  double z = v.y() + OffsetMag; 

//  cout << "In Get function FLUKA coord x=" << x << "  y=" << y << "  z=" << z -OffsetMag << endl;
	
  double BxFLUKA;
  double ByFLUKA;
  double BzFLUKA; 
  Get(x,y,z,BxFLUKA,ByFLUKA,BzFLUKA);
//Convert from FLUKA to KalTest, magnetic field vector must point in the +Z direction
	
 double Bx = ByFLUKA;			// in [T]
 double By = BzFLUKA;			// in [T]
 double Bz = BxFLUKA;			// in [T]

	
/////////////////SANITY CHECK: SET ALL ENTRIES TO 0.3 T, TEST THAT RESULTS ARE SAME AS RECO IN UNI B////////
	/*
	double Bx = 0;
	double By = 0;
	double Bz = 0.3;
	*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
   double Bmag = TMath::Sqrt((Bx*Bx) + (By*By) + (Bz*Bz));
// cout<<"FLUKA coord Bx = " <<BxFLUKA<<", By = "<<ByFLUKA<< " , Bz = " << BzFLUKA << "  , Bmag = " << Bmag <<"\n";
 //cout<<"KALTEST coord GetGlobalBField:  Bx = " <<Bx<<", By = "<<By<< " , Bz = " << Bz << "  , Bmag = " << Bmag <<"\n";
 return TVector3(Bx, By, Bz);
}


TVector3 TBField::GetGlobalBfield(const TVector3& globalPosition)
{
	if(!fUseUniformBfield) {
		//cout << "non-uniform magnetic field! " << endl;
		TVector3 bfield = Get(globalPosition);

		return bfield;
	}

// AESOPLITE average magnetic field
	Double_t B0   = 0.30;  //in T

	if(fUseUniformBfield)
	{
		return TVector3(0,0,B0);
	}


}
