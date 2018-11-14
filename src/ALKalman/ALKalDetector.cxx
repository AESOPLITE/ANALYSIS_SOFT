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
#include "LoadMCparameters.h"
#include <sstream>
#include <cmath>

using namespace std;


  
const int DetLayer = 7;
//const Double_t TckZPosTop[]= {-14.6, -34.6, -94.6, -154.6, -174.6, -194.6, -214.6};	            //in mm
//const Double_t TckZPosBottom[]= {-15.0, -35.0, -95.0, -155.0, -175.0, -195.0, -215.0};	        //in mm
const Double_t width = 90;					  //width from center of layer in mm
const Double_t thick = -0.4;					 //thickness of layer in mm
const Double_t step = -1;					//step between dummy layers in mm
const Double_t xoffset[] = {}; 
const Double_t yoffset[] = {};
Double_t ALKalDetector::fgBfield = 0.3;			//in T	
Bool_t active = ALMeasLayer::kActive;
Bool_t dummy  = ALMeasLayer::kDummy;		
Bool_t bending = ALMeasLayer::kBending;
Bool_t nonbending  = ALMeasLayer::kNonBending;	
Bool_t inuse = ALMeasLayer::kInUse;
Bool_t notinuse  = ALMeasLayer::kNotInUse;	
ClassImp(ALKalDetector)


ALKalDetector::ALKalDetector(Int_t m): TVKalDetector(m)
{


 //Load region numbers used in the MC geometry 
 int*TckReg=new int[7];
 int*TrigReg=new int[4];
 int*GReg=new int[1];
 int*ShellReg=new int[2];
 float*TckZPosTop=new float[7];
 float*TckZPosBottom=new float[7];
 float*TrigThresh=new float[4];
 float*GuardThresh=new float[1];
 for(int i=0;i<7;i++)TckReg[i]=0;
 for(int i=0;i<4;i++)TrigReg[i]=0;
 for(int i=0;i<1;i++)GReg[i]=0;
 for(int i=0;i<7;i++)TckZPosTop[i]=0;
 for(int i=0;i<7;i++)TckZPosBottom[i]=0;
 for(int i=0;i<4;i++)TrigThresh[i]=0;
 for(int i=0;i<1;i++)GuardThresh[i]=0;
  string MCparamfile="../src/ALSim/MCparameters.dat"; 
 LoadMCparameters(MCparamfile,TckReg,TrigReg,GReg,TckZPosTop,TrigThresh,GuardThresh,ShellReg);
 for(int i=0;i<7;i++){
  TckZPosTop[i]= TckZPosTop[i]*10;   //convert in mm
  TckZPosBottom[i]=TckZPosTop[i]-0.4;
  }
	
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

		for (Int_t layer = 0; layer < DetLayer; layer++) {
		if ((layer == 0) || (layer == 4) || (layer == 6)) {
			TVector3 xc(0., TckZPosTop[layer], 0.);	
			TVector3 normal(0., 1.0, 0.);		
		    Add(new ALMeasLayer(Si, Ni, xc, normal, active, nonbending, notinuse));	
	
			
		//add dummy layers, start from bottom of previous active layer
			if(layer==6) break;
		//how many dummies can we add?
			int ndummies = floor((TckZPosBottom[layer]-TckZPosTop[layer+1])/TMath::Abs(step));
			TVector3 xdummy;
			double temp = TckZPosBottom[layer];
		for(int i=0; i<ndummies;i++) {
			xdummy.SetXYZ(0, temp, 0);
			Add(new ALMeasLayer(Ni, Ni, xdummy, normal,dummy, nonbending, notinuse));
            temp+=step; 
		} 

	}
		else {
	    TVector3 xc(0., TckZPosTop[layer], 0.);	
		TVector3 normal(0.,1.0, 0);			
		Add(new ALMeasLayer(Si, Ni, xc, normal, active, bending, notinuse));
		
	
			
			
//add dummy layers, start from bottom of previous active layer
		int ndummies = floor((TckZPosBottom[layer]-TckZPosTop[layer+1])/TMath::Abs(step));
		TVector3 xdummy;
		double temp = TckZPosBottom[layer];
	for(int i=0; i<ndummies;i++) {
		xdummy.SetXYZ(0, temp, 0);
		Add(new ALMeasLayer(Ni, Ni, xdummy, normal,dummy, bending, notinuse));
		temp+=step; 
		} 

	}
	SetOwner();
 }
}

ALKalDetector::~ALKalDetector()
{
}
