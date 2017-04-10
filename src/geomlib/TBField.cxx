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

TEveMagField* TBField::fField = 0;

Bool_t   TBField::fUseUniformBfield = kFALSE;
//Bool_t   TBField::fUseUniformBfield = kTRUE;
Double_t TBField::fFieldCoeff       = 1.;

TH3D* TBField::fieldmapX = new TH3D("field map Bx" , "Bx [kG]; x[cm]; y[cm]; z[cm]", 80 , -200, 200, 80, -200, 200, 10, -25, 25);
TH3D* TBField::fieldmapY = new TH3D("field map By" , "By [kG]; x[cm]; y[cm]; z[cm]", 80 , -200, 200, 80, -200, 200, 10, -25, 25);
TH3D* TBField::fieldmapZ = new TH3D("field map Bz" , "Bz [kG]; x[cm]; y[cm]; z[cm]", 80 , -200, 200, 80, -200, 200, 10, -25, 25);

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
	
	//do coordinate transform here when reading text file (x'=y , y'=z, z'=x, same for vector field)
	
     TNtuple *T= new TNtuple("fieldmap","fieldmap","z:x:y:bz:bx:by");
     T->ReadFile("/home/sarah/AESOPLITE/FLUKA/fieldmap.txt");
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
	//	cout<< x <<" ,"<<y <<" ,"<< z <<" ,"<<bx << " ," << by << " ," << bz << endl;
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
	cout << "Field Map INTERFACED WOOO" << endl;
	return true;
}


void TBField::Get(const double&x, const double&y, const double&z, double& Bx, double& By, double& Bz) {

	
	int binx = fieldmapX->GetXaxis()->FindBin(x);
	int biny = fieldmapX->GetYaxis()->FindBin(y);
	int binz = fieldmapX->GetZaxis()->FindBin(z);
	
	 Bx = (fieldmapX->GetBinContent(binx, biny, binz))/(10000);					//convert from [G] to [T]
	 By = (fieldmapY->GetBinContent(binx, biny, binz))/(10000);
     Bz = (fieldmapZ->GetBinContent(binx, biny, binz))/(10000);
	
//cout<<"DEBUG: "<<__LINE__<<": "<<x <<","<< y <<","<<binx<<" ,"<<biny<<": "<<Bx<<","<<By<< " ," << Bz <<"\n";

}


TVector3 TBField::Get(const TVector3& v) {
  double x = v.x();
  double y = v.y();
  double z = v.z();
  double Bx;
  double By;
  double Bz;
  Get(x,y,z,Bx,By,Bz);
 // cout << "Get function x=" << x << "  , y=" << y << " , z=" << z << endl;
//  cout << "Get function Bx=" << Bx << " , By=" << By << " , Bz=" << Bz << endl;
  return TVector3(Bx, By, Bz);
}


TVector3 TBField::GetGlobalBfield(const TVector3& globalPosition)
{
	if(!fUseUniformBfield) {
		//cout << "non-uniform magnetic field! " << endl;
		TVector3 bfield = Get(globalPosition);

		return bfield;
	}

	// ILD uniform magnetic field
	Double_t B0   = 0.35; 

	if(fUseUniformBfield)
	{
		return TVector3(0,0,B0);
	}

	// an artificial non-uniform magnetic field
	//cout << "using that part of function GetGlobalBField..." << endl;
    const Double_t zmax    = 3000.;
    const Double_t rmax    = 3000.;
    const Double_t coeffxy = fFieldCoeff / (zmax * rmax);
    const Double_t bmax    = 0.35; // B at the origin

    if(fUseUniformBfield)
    {
        return TVector3(0,0,bmax);
    }

    Double_t xg = globalPosition.X();
    Double_t yg = globalPosition.Y();
    Double_t zg = globalPosition.Z();

    Double_t bx = bmax * coeffxy * zg * xg; 
    Double_t by = bmax * coeffxy * zg * yg; 
    Double_t bz = bmax * (1. - coeffxy * zg * zg);

	TVector3 bfield(bx, by, bz);

    return bfield;
	
}