#ifndef __ALHIT__
#define	__ALHIT__
#include "KalTrackDim.h"
#include "TVTrackHit.h"
#include "ALMeasLayer.h"

//*************************************************************************
//* ===================
//*  ALHit Class
//* ===================
//*
//* (Description)
//*  Sample measurement vector class as defined by ALMeasLayer for the AESOPLITE instrument 
//*  Describes the plane silicon detector
//* (Requires)
//* (Provides)
//*     class ALHit
//* (Update Recorded)
//*
//*   22/12/2016  S. Mechbal		Modified for AESOPLITE experiment
//***************************************************************************

class ALHit : public TVTrackHit {
public :
	
	//ALHit () {}		//default constructor
	ALHit (Int_t m = kMdim);
	
	ALHit(const ALMeasLayer &ms,
		  Double_t 			*x,
		  Double_t			*dx,
		 const TVector3		&xx,
		  Double_t 			  b,
		  Int_t 	          m= kMdim);
	
	virtual ~ALHit();
	
	virtual TKalMatrix XvToMv (const TVector3 &xv, Bool_t isbending) const;
	//virtual void 		DebugPrint(Option_t *opt = "")			const;
	
   // return raw position vector
   void      SetRawXv(const TVector3 &v3) {fRawxv = v3;}
   TVector3  GetRawXv() const {return fRawxv;}
private:
	TVector3 fRawxv;			//exact hit position
	
	ClassDef(ALHit, 1)		//Sample hit class
};

#endif