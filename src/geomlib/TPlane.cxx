//*************************************************************************
//* ====================
//*  TPlane Class
//* ====================
//*
//* (Description)
//*   A class to implement a plane object.
//* (Requires)
//*     TVSurface
//* (Provides)
//*     class TPlane
//* (Update Recored)
//*   2004/10/30  A.Yamaguchi   Original version.  Currently fXc is 
//*                             supposed to be at the origin
//*
//*************************************************************************
//
#include "TPlane.h"
#include<iostream>

using namespace std;

//_____________________________________________________________________
//  -----------------------------------
//  Plane Class
//  -----------------------------------

ClassImp(TPlane)

TPlane::TPlane()
      : fXc(), fNormal()
{
} 

TPlane::TPlane(const TVector3 &xc)
      : fXc(xc), fNormal(xc.Unit())
{
} 

TPlane::TPlane(const TVector3 &xc, const TVector3 &n)
      : fXc(xc), fNormal(n.Unit())
{
} 

//_____________________________________________________________________
//  -----------------------------------
//  Calculate S
//  -----------------------------------
//
Double_t TPlane::CalcS(const TVector3 &xx) const
{
   return (xx - fXc) * fNormal;
}

//_____________________________________________________________________
//  -----------------------------------
//  Calculate (@S/@x)
//  -----------------------------------
//
TMatrixD TPlane::CalcDSDx(const TVector3 &xx) const
{
   TMatrixD dsdx(1,3);
   dsdx(0,0) = fNormal.X(); 
   dsdx(0,1) = fNormal.Y(); 
   dsdx(0,2) = fNormal.Z();
   return dsdx;
}

Bool_t TPlane::IsOnSurface(const TVector3 &xx) const
{
	//cout << "Calling TPlane IsOnSurface() " << endl;
//	cout << "center vector of surface fXc("<<fXc.X()<<","<<fXc.Y()<<","<<fXc.Z()<<")"<<endl;

#if 0
   return (xx - fXc) * fNormal == 0. ? kTRUE : kFALSE; 
#else
   return kTRUE;
#endif
} 


