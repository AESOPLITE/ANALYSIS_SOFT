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

const Int_t DetDummyLayer = 7; 				
const Int_t DetLayer = 7;
const Double_t DetLayerList[]= {-1.46, -3.46, -9.46, -15.46, -17.46, -19.46, -21.46};	            //list of z postiions of top part of layer detectors in cm
const Double_t DetLayerListDummy[]= {-1.50, -3.50, -9.50, -15.50, -17.50, -19.50, -21.50};	            //list of z positions of bottom part of detectors in cm
const Double_t width = 9.0;					     //width from center of layer 											//in cm
const Double_t thick = -0.04;					 //thickness of layer																//in cm
const Double_t xoffset[] = {}; 
const Double_t yoffset[] = {};
Double_t ALKalDetector::fgBfield = 3;	
	
Bool_t active = ALMeasLayer::kActive;
Bool_t dummy  = ALMeasLayer::kDummy;										
ClassImp(ALKalDetector)


ALKalDetector::ALKalDetector(Int_t m): TVKalDetector(m)
{
	
	//Silicon material properties from PDG 

	Double_t A, Z, density, radlen;
	A = 28.0855;							// mass number
	Z = 14.0;								//atomic number
	density = 2.33;							//in g/cm3
	radlen = 9.37;							//in cm
	TMaterial &Si = *new TMaterial("Si", "", A, Z, density, radlen, 0.);
	
	//define air as material 
	
	  Double_t A_prime       =  14.00674 * 0.7 + 15.9994 * 0.3;
      Double_t Z_prime       =  7.3;                               // atomic number
      Double_t density_prime = 1.205e-3;                        // [g/cm^3]
  	  Double_t radlen_prime  = 3.42e4;                          // [cm]
  	  TMaterial &air = *new TMaterial("Air", "", A_prime, Z_prime, density_prime, radlen_prime, 0.);
	
    //define nitrogen as material 
	
	  Double_t A_second       =  14.00674;
      Double_t Z_second      =   7.0;                               // atomic number
      Double_t density_second = 1.17e-3;                        // [g/cm^3]
  	  Double_t radlen_second  = 3.26e4;                          // [cm]
  	  TMaterial &Ni = *new TMaterial("Ni", "", A_second, Z_second, density_second, radlen_second, 0.);
	
	//add measurement layers in bending/non bending plane

	/////////////////////////////////OLD CONFIGURATION/////////////////////////////////////
		for (Int_t layer = 0; layer < DetLayer; layer++) {
		if ((layer == 0) || (layer == 4) || (layer == 6)) {
			TVector3 xc(0., DetLayerList[layer], 0.);	
			TVector3 normal(0., 1.0, 0.);						
			Add(new ALMeasLayer(Si, Ni, xc, normal, ALMeasLayer::kActive, ALMeasLayer::kNonBending));
			xc.SetXYZ(0, (DetLayerList[layer] + thick), 0);
		    Add(new ALMeasLayer(Ni, Si, xc, normal, ALMeasLayer::kDummy, ALMeasLayer::kNonBending));

			
		} 
		else {
	    TVector3 xc(0., DetLayerList[layer], 0.);	
		TVector3 normal(0.,1.0, 0);						
		Add(new ALMeasLayer(Si, Ni, xc, normal, ALMeasLayer::kActive, ALMeasLayer::kBending));
	    xc.SetXYZ(0, (DetLayerList[layer] + thick), 0);
	    Add(new ALMeasLayer(Ni, Si, xc, normal, ALMeasLayer::kDummy, ALMeasLayer::kBending));

			
		}
	}		
	SetOwner();
}

ALKalDetector::~ALKalDetector()
{
}