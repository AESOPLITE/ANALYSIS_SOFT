//*************************************************************************
//* ===================
//*  ALKalDetector Class
//* ===================
//*
//* (Description)
//*  Sample detector layers class for the AESOPLITE instrument 
//*  Describes configuration of seven plane silicon detector
//* (Requires)
//* (Provides)
//*     class ALKalDetector
//* (Update Recorded)
//*
//*   22/12/2016  S. Mechbal		Modified for AESOPLITE experiment
//***************************************************************************

#include "ALKalDetector.h"
#include "ALMeasLayer.h"
#include "ALHit.h"
#include "TRandom.h"
#include "TMath.h"
#include "TVector3.h"
#include <sstream>

const Int_t DetDummyLayer = 3; 		//let's not really worry about this right now
const Int_t DetLayer = 7;
const Double_t DetLayerList[]= {-1.5, -3.5, -9.5, -15.5, -17.5, -19.5, -21.5};	            //list of z postions of detectors in cm
const Double_t xwidth = 18.0;																//in cm
const Double_t zwidth = 00.4;																//in cm
const Double_t xoffset[] = {}; 
const Double_t yoffset[] = {};
Double_t ALKalDetector::fgBfield = 0.35;													//in T
ClassImp(ALKalDetector)


ALKalDetector::ALKalDetector(Int_t m): TVKalDetector(m)
{
	//Silicon material properties from PDG 
	
	Double_t A, Z, density, radlen;
	A = 28.0855;							// mass number
	Z = 14.0;								//atomic number
	density = 2.33;							//in gm/cm3
	radlen = 9.37;							//in cm
	TMaterial &Si = *new TMaterial("Si", "", A, Z, density, radlen, 0.);
	
	//define air as material 
	
	  Double_t A_prime       =  14.00674 * 0.7 + 15.9994 * 0.3;
      Double_t Z_prime       =  7.3;                               // atomic number
      Double_t density_prime = 1.205e-3;                        // [g/cm^3]
  	  Double_t radlen_prime  = 3.42e4;                          // [cm]
  	  TMaterial &air = *new TMaterial("Air", "", A_prime, Z_prime, density_prime, radlen_prime, 0.);
	
	//we now add the 7 layers 
	for (Int_t layer = 0; layer < DetLayer; layer++) {
	//reference vector on layer pointing towards the center (unclear whether this is correct) 
		if ((layer == 0) || (layer == 4) || (layer == 6)) {
			TVector3 xc(0., DetLayerList[layer], 0.);	
			TVector3 normal(0., 1.0, 0.);						//normal vector to the layer
			Add(new ALMeasLayer(Si, air, xc, normal, ALMeasLayer::kActive, ALMeasLayer::kNonBending));
		}
		else {
	    TVector3 xc(0., DetLayerList[layer], 0.);	
		TVector3 normal(0.,1.0, 0);						//normal vector to the layer
		Add(new ALMeasLayer(Si, air, xc, normal, ALMeasLayer::kActive, ALMeasLayer::kBending));
	 }
	}	
	SetOwner();
}

ALKalDetector::~ALKalDetector()
{
}