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

#include "ALHit.h"
#include "ALMeasLayer.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>

using namespace std;


ClassImp(ALHit)

//Constructors and destructors 

ALHit::ALHit(Int_t m)
	: TVTrackHit(m) {}

ALHit::ALHit(const ALMeasLayer &ms,			//measurement layer
			 Double_t 			*x,			//coordinate array
			 Double_t			 *dx,		//error on coordinate array
	   const TVector3 			 &xx,		//measurement vector
			 Double_t 			b,			//magnetic field
			 Int_t				m)			//dimension of measurement vector (m=2 in planar detector)
	: TVTrackHit(ms, x, dx, b, m)
{
}



ALHit::~ALHit()
{
}

//Implementation of public methods 

TKalMatrix ALHit::XvToMv(const TVector3 &xv, Bool_t isbending) const
{		
	
		TKalMatrix mv(kMdim, 1);

	  const ALMeasLayer &ml = dynamic_cast<const ALMeasLayer &>(GetMeasLayer());	
      Bool_t bendingplane = ml.IsBending();
	if (bendingplane) {
		mv(0,0) = xv.X();
		mv(1,0) = xv.Y();
	}
	else {
		mv(0,0) = xv.Y();
		mv(1,0) = xv.Z();
	}

	return mv;
}

