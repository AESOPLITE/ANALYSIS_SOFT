#ifndef __ALKALDETECTOR__
#define __ALKALDETECTOR__

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

#include "TVector3.h"			//from ROOT
#include "TVKalDetector.h"		//from KalTrackLib
#include "ALMeasLayer.h"
#include "TAttDrawable.h"
#include "TBField.h"


class TVMeasLayer;
class TNode;

extern const int kNDetDummyLayer;
extern const int knDetLayer;	// = 7
extern const double kDetLayerList[];

class ALKalDetector : public TVKalDetector {
public:
   // Ctor and Dtor
    //ALKalDetector() {};			//default constructor
   ALKalDetector(Int_t m = 100); //m is the maximum number of measurement layers
   ~ALKalDetector();

   // Utility methods
   static Double_t GetBfield (const TVector3 &xx)
                               { return  TBField::GetGlobalBfield(xx).Mag(); }


private:
 static Double_t fgBfield;   // magnetic field [kG]

   ClassDef(ALKalDetector,1)   // Sample hit class
};

#endif

