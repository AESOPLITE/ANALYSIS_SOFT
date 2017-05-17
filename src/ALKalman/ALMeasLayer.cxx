//*************************************************************************
//* ===================
//*  ALMeasLayer Class
//* ===================
//*
//* (Description)
//*  Sample measurement layer class for the AESOPLITE instrument used by ALHIT
//*  Describes the plane silicon detector
//* (Requires)
//* (Provides)
//*     class ALMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*   2011/06/17  D.Kamai           Modified to handle ladder structure.
//*
//*   22/12/2016  S. Mechbal		Modified for AESOPLITE experiment
//***************************************************************************


#include <iostream>
#include "ALMeasLayer.h"
#include "ALHit.h"
#include "ALKalDetector.h"
#include "TVTrack.h"
#include "TRandom.h"
#include "TMath.h"
#include "TRotMatrix.h"
#include "TVirtualPad.h"
#include "TBRIK.h"
#include "TNode.h"
#include "TString.h"

using namespace std;

#define _ALMeasLayerDebug_ 1

Bool_t   ALMeasLayer::kActive = kTRUE;
Bool_t   ALMeasLayer::kDummy = kFALSE;
Bool_t	 ALMeasLayer::kBending = kTRUE;
Bool_t   ALMeasLayer::kNonBending = kFALSE;
Double_t ALMeasLayer::fSigmaX = 0.0066;				//in cm	
Double_t ALMeasLayer::fSigmaY = 0.0115;				//in cm 
Double_t ALMeasLayer::fSigmaZ = 0.0066; 			   //in cm

ClassImp(ALMeasLayer)
	
ALMeasLayer::ALMeasLayer(TMaterial &min,			//material inside the layer
                         TMaterial &mout,		    //material outside the layer
			       const TVector3  &center,			//reference point on the layer (xc, yc, zc)
			       const TVector3  &normal,			//outward normal vector 
                   Bool_t     isactive,				//flag to tell if the layer is active
				   Bool_t	  isbending)				//flag to tell if the layer is in bending plane or non-bending
              : TVMeasLayer(min, mout, isactive),
                TPlane(center, normal), fIsBending(isbending)
		
{
}

ALMeasLayer::~ALMeasLayer()
{}

TKalMatrix ALMeasLayer::XvToMv(const TVector3 &xv, Bool_t isbending) const
{
	
	//Change of coordinates: bending plane (XY), non bending-plane (YZ) //   x'=y ; y'=z ; z'=x 
	
	TKalMatrix mv(kMdim, 1);
	Float_t X = xv.Y();
    Float_t Y = xv.Z();
	Float_t Z = xv.X();

    if (isbending) {
	mv(0,0) = X;
	mv(1,0) = Y;
	return mv;
}
	else {
		mv(0,0) = Y;
		mv(1,0) = Z;
    	return mv;
  
	}
	
}

TKalMatrix ALMeasLayer::XvToMv(const TVTrackHit &, const TVector3 &xv, Bool_t isbending) const
{return XvToMv(xv, isbending);}

TVector3 ALMeasLayer::HitToXv(const TVTrackHit &vht) const
{
	const ALHit &ht = dynamic_cast<const ALHit &>(vht);
	Bool_t bending = IsBending();
	TVector3 xraw= ht.GetRawXv();
	
	if (bending) {
	Double_t x = ht(0,0);
	Double_t y = ht(1,0);
	//Double_t z = GetXc().Z();								//default value of z-coordinate center of the plane for bending-plane, random but will be given huge error
	Double_t z = xraw.Z();
	return TVector3(x,y,z);
	
	}
	else {
	Double_t y = ht(0,0);
	Double_t z = ht(1,0);
    //Double_t x = GetXc().X();									//default value of x-coordinate center of the plane for non-bending plane, random but will be given huge error
    Double_t x = xraw.X();	
		return TVector3(x,y,z);
	}
}

//HitToXv for MC raw hits, no suppression of third coordinate 


TVector3 ALMeasLayer::RawHitToXv(const TVTrackHit &vht) const
{
	const ALHit &ht = dynamic_cast<const ALHit &>(vht);  
	Bool_t isbending = IsBending();   
	TVector3 xraw= ht.GetRawXv();
	
	if (isbending) {
	Double_t x = ht(0,0);	
	Double_t y = ht(1,0);
	Double_t z = xraw.Z();
   //cout << "Raw hit x=" << x << "   y=" << y << "  z=" << z << endl;
	return TVector3(x,y,z);

	}
	
    else {
	Double_t y = ht(0,0);
	Double_t z = ht(1,0);
    Double_t x = xraw.X();	
  // cout << "Raw hit x=" << x << "   y=" << y << "  z=" << z << endl;
   return TVector3(x,y,z);
	}
}





void ALMeasLayer::CalcDhDa(const TVTrackHit &vht, const TVector3 &xxv,
						   const TKalMatrix &dxphiada, TKalMatrix &H) const
{
   // Calculate
   //    H = (@h/@a) = (@phi/@a, @z/@a)^t
   // where
   //        h(a) = (phi, z)^t: expected meas vector
   //        a = (drho, phi0, kappa, dz, tanl, t0)
   //
	
	Int_t sdim = H.GetNcols();
	Int_t hdim = TMath::Max(5, sdim-1);
	Bool_t isbending = IsBending();
	
	//cout << "CalcDhDa isbending? " << isbending << endl;
	
   for (Int_t i=0; i<hdim; i++) {
	   if (isbending) {
      H(0,i) = dxphiada(0,i);
      H(1,i) = dxphiada(1,i);
   }
	
	  else {
	  H(0,i) = dxphiada(1,i);
      H(1,i) = dxphiada(2,i);
   	}
  }
  
   if (sdim == 6) {
      H(0,sdim-1) = 0.;
      H(1,sdim-1) = 0.;
   }
	//cout << "projector matrix H ..." << endl;
	//H.DebugPrint();
}





void ALMeasLayer::ProcessHit(const TVector3 &xx, TObjArray &hits, Bool_t isbending, Int_t index)
{
	//we proceed to a change of coordinate here to fit the helix parametric equation. The magnetic field points now in the +z direction.
	//bending plane is the XY plane

	

	if (isbending) {
	TKalMatrix h = XvToMv(xx, isbending); 
	Int_t m = 2;
	Double_t xraw = xx.Y();						//change of coordinates
    Double_t yraw = xx.Z();						//change of coordinates 	
	Double_t zraw = xx.X();						//change of coordinates 
	TVector3 xv;
	xv.SetXYZ(xraw,yraw,zraw);
	Double_t x = h(0,0);
	Double_t y = h(1,0);
	Double_t dx = GetSigmaX();				//rms of x in cm
	Double_t dy = GetSigmaY();				//rms of y in cm (400um/sqrt(12))
    x  += gRandom->Gaus(0., dx);   // smearing rphi
    y += gRandom->Gaus(0., dy);   // smearing z	
	Double_t meas[2];
	Double_t dmeas[2];
	meas[0] = x;
	meas[1] = y;
	dmeas[0] = dx;
	dmeas[1] = dy;
	Double_t b = ALKalDetector::GetBfield(xv);
	ALHit *aHit = new ALHit(*this, meas, dmeas, xx, b, m);
	aHit->SetRawXv(xv);
	hits.AddAt(aHit, index);

	}
	//non-bending plane is YZ plane
	else {
	TKalMatrix h = XvToMv(xx, isbending);
	Int_t m = 2;
	Double_t xraw = xx.Y();						//change of coordinates
    Double_t yraw = xx.Z();						//change of coordinates 	
	Double_t zraw = xx.X();						//change of coordinates 
	TVector3 xv;
	xv.SetXYZ(xraw,yraw,zraw);
	Double_t y = h(0,0);
	Double_t z = h(1,0);
	Double_t dy = GetSigmaY();				
	Double_t dz = GetSigmaZ();				
    y  += gRandom->Gaus(0., dy);   // smearing rphi
    z += gRandom->Gaus(0., dz);   // smearing z	
	Double_t meas[2];
	Double_t dmeas[2];
	meas[0] = y;
	meas[1] = z;
	dmeas[0] = dy;
	dmeas[1] = dz;	
	Double_t b = ALKalDetector::GetBfield(xv);		//won't be used for non-uniform mag field 
	ALHit *aHit = new ALHit(*this, meas, dmeas, xx, b,m);
	aHit->SetRawXv(xv);
	hits.AddAt(aHit, index);
	}
}


void ALMeasLayer::ProcessRawHit(const TVector3 &xx, TObjArray &hits, Bool_t isbending) {
	
	//cout << "Using function ProcessRawHit" << endl;
	TKalMatrix h(3,3);
	Int_t m = 2;
	Double_t x = xx.Y();						//change of coordinates
    Double_t y = xx.Z();						//change of coordinates 	
	Double_t z = xx.X();						//change of coordinates 
	TVector3 xv;
	xv.SetXYZ(x,y,z);
	Double_t dx = GetSigmaX();
	Double_t dy = GetSigmaY();
	Double_t dz = GetSigmaZ();
	h(0,0) = x;
	h(1,0) = y;
	h(2,0) = z;
	
	//cout << "ProcessRawHit  x=" << x << "  y=" << y << " z=" << z << endl;
	
	if (isbending) {
	Double_t meas[2];
	Double_t dmeas[2];
    meas[0] = x;
	meas[1] = y;
	dmeas[0] = dx;
	dmeas[1] = dy;
	Double_t b = ALKalDetector::GetBfield(xv);
	ALHit *RawHit = new ALHit(*this, meas, dmeas, xx, b, m);
    RawHit->SetRawXv(xv);
    hits.Add(RawHit);
	}
	
	else {
    Double_t meas[2];
	Double_t dmeas[2];
	meas[0] = y;
	meas[1] = z;
	dmeas[0] = dy;
	dmeas[1] = dz;
	Double_t b =ALKalDetector::GetBfield(xv);
	ALHit *RawHit = new ALHit(*this, meas, dmeas, xx, b, m);
    RawHit->SetRawXv(xv);
	hits.Add(RawHit);
	}
}
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
