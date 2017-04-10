#ifndef TBFIELD_H
#define TBFIELD_H
//*************************************************************************
//* =====================
//*  TBField Class
//* =====================
//*
//* (Description)
//*   A sigleton to hold information of detector system
//*   used in Kalman filter classes.
//* (Requires)
//* 	TObject
//* (Provides)
//* 	class TBField
//* (Update Recored)
//*   2013/01/31  Bo Li	 Original Version.
//*
//*************************************************************************

#include "TObject.h"     // from ROOT
#include "TVector3.h"
#include "TEveTrackPropagator.h"     // from ROOT
#include <TFile.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH3D.h>
//_____________________________________________________________________
//  ------------------------------
//   Detector system class
//  ------------------------------

class TBField : public TObject {
public:
   TBField();
   virtual ~TBField();

   // Utility methods
   static TVector3 GetGlobalBfield(const TVector3& globalPosition);
   static bool SetMagField();
   static void SetBfieldPtr(TEveMagField* b){ fField = b; } 
	static void Get(const double& posX, const double& posY, const double& posZ, double& Bx, double& By, double& Bz);
	static TVector3 Get(const TVector3& pos);
   //static void SetBfield   (TVector3      b){ fBfield   = b; }
   static void SetUseUniformBfield(Bool_t b) { fUseUniformBfield = b;    }
   static void SetBfieldCoeff(Double_t k)    { fFieldCoeff       = k;    }

   static bool IsUsingUniformBfield()        { return fUseUniformBfield; }

private:
   //static    TVector3      fBfield;     //a constant b field

   static TEveMagField* fField;   //a field map
   static Bool_t   fUseUniformBfield;
   static Double_t fFieldCoeff;
   static TH3D* fieldmapX;
   static TH3D* fieldmapY;
   static TH3D* fieldmapZ;

   ClassDef(TBField,1)  // Base class for detector system
};

#endif