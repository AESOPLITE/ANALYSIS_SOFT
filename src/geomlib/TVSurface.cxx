//*************************************************************************
//* ====================
//*  TVSurface Class
//* ====================
//*
//* (Description)
//*   This is the base class for various surfaces.
//* (Requires)
//*     TObject;
//* (Provides)
//*     class TVSurface
//* (Update Recored)
//*   2003/10/03  K.Fujii       Original version.
//*   2005/02/23  K.Fujii       Added new methods, Compare() and
//*                             GetSortingPolicy().
//*
//*************************************************************************
//
#include <iostream>
#include "TVSurface.h"
#include "TVTrack.h"

using namespace std;

//_____________________________________________________________________
//  -----------------------------------
//  Base Class for any surface
//  -----------------------------------

ClassImp(TVSurface)

//_____________________________________________________________________
//  -----------------------------------
//  Calculate crossing point with track
//  -----------------------------------
//
Int_t TVSurface::CalcXingPointWith(const TVTrack  &hel,
                                         TVector3 &xx,
                                         Double_t &phi,
                                         Double_t  eps) const
{
   return CalcXingPointWith(hel,xx,phi,0,eps);
}

Int_t TVSurface::CalcXingPointWith(const TVTrack  &hel,
				         TVector3 &xx,
				         Double_t &phi,
				         Int_t     mode,
				         Double_t  eps) const
{
//   cout << "Calling CalcXingPointWith() function, mode = " << mode << endl;
   static const Int_t       maxcount   = 100;
   static const Double_t    initlambda = 1.e-10;
   static const Double_t    lambdaincr = 10.;
   static const Double_t    lambdadecr = 0.1;
   
   //cout << "Initial xx at surface x = " << xx.X() << "  y = " << xx.Y() << "  z = " << xx.Z() << endl;
   xx = hel.CalcXAt(phi);
   Double_t  lastphi =  phi;
   Double_t  lasts   =  1.e10;
   Double_t  lambda  =  initlambda;

   TVector3  lastxx  =  xx;
   Int_t     count   =  0;

   Double_t s;

   while (1) {
      if (count > maxcount) {
         s       = lasts;
         phi     = lastphi;
         xx      = lastxx;
#if 0
         cout << "TVSurface::CalcXingPointWith:"
              << "   --- Loop count limit reached ---------- " << endl
              << "   phi    : " << phi    << endl
              << "   x      : " << xx.X() << " "
                                << xx.Y() << " "
                                << xx.Z() << endl 
              << "   s      : " << s      << endl
			  << "   eps    : " << eps    << endl
              << "   lambda : " << lambda << endl;
#endif
		//  cout << "TVSurface has reached the max number of count, return 0 " << endl;
         return 0;
      }
      count++;
      s  = CalcS(xx);
      if (TMath::Abs(s) < eps) break;
      if (TMath::Abs(s) < TMath::Abs(lasts)) {
         lasts   = s;
         lastphi = phi;
         lastxx  = xx;
         lambda *= lambdadecr;
      } else {
         s       = lasts;
         phi     = lastphi;
         xx      = lastxx;
         lambda *= lambdaincr;
      }
      TMatrixD dsdx   = CalcDSDx(xx);
      TMatrixD dxdphi = hel.CalcDxDphi(phi);
      TMatrixD dsdphi = dsdx * dxdphi;
      Double_t denom = (1 + lambda) * dsdphi(0,0);
      phi -= s / denom;
      xx   = hel.CalcXAt(phi);
	  
   }

   if( mode!=0 && hel.IsInB()){ // (+1,-1) = (fwd,bwd)
     const Int_t chg = (Int_t)TMath::Sign(1.1, hel.GetKappa());
     if( chg*phi*mode/hel.GetPtoR() > 0){
	//	 cout << "quit function at line 112 " << endl; 
       return 0;
     }
   }
   //cout << "TVSurface vector xx x = " << xx.X() << "  y = " << xx.Y() << " z = " << xx.Z() << endl;
   
  // cout << "TVSurface calling IsOnSurface() " << endl;
   return (IsOnSurface(xx) ? 1 : 0);
  
}

//_____________________________________________________________________
//  -----------------------------------
//  Compare to Surfaces
//  -----------------------------------
//
Int_t TVSurface::Compare(const TObject *obj) const
{
   Double_t me  = GetSortingPolicy();
   Double_t you = dynamic_cast<const TVSurface *>(obj)->GetSortingPolicy();
   return me < you ? -1 : (me > you ? +1 : 0);
}
