#ifndef __ALMEASLAYER__
#define __ALMEASLAYER__
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
//
#include "TVector3.h"
#include "TKalMatrix.h"
#include "TPlane.h"
#include "TVMeasLayer.h"
#include "KalTrackDim.h"
#include "TMath.h"
#include <sstream>
class TVTrackHit;

class ALMeasLayer : public TVMeasLayer, public TPlane {
	public:

	static Bool_t kActive;
	static Bool_t kDummy;
	static Bool_t kBending;
	static Bool_t kNonBending;
	static Bool_t kInUse;
	static Bool_t kNotInUse;
	//Construtors and destructor 


	ALMeasLayer(TMaterial &min,
				TMaterial &mout,
				const TVector3 &center,
				const TVector3 &normal,
				Bool_t isactive = ALMeasLayer::kActive,
				Bool_t isbending = ALMeasLayer::kBending,
         	    Bool_t isinuse = ALMeasLayer::kNotInUse);
	
	
	ALMeasLayer() {}			//default constructor
	
	virtual ~ALMeasLayer();
	
	//Pure virtual functions to be implemented
	
	virtual TKalMatrix XvToMv (const TVTrackHit &ht, const TVector3 &xv, Bool_t isbending) const;
	virtual TKalMatrix XvToMv (const TVector3 &xv, Bool_t isbending) const; 		
	virtual TVector3 HitToXv (const TVTrackHit &ht) const;
	virtual TVector3 RawHitToXv (const TVTrackHit &ht) const;
	virtual void CalcDhDa (const TVTrackHit &ht, const TVector3 &xv, const TKalMatrix &dxphiada, TKalMatrix &H) const;
	//virtual void ProcessHit(const TVector3 &xx, TObjArray &hits, Bool_t isbending);
	virtual void ProcessHit(const TVector3 &xx, TObjArray &hits, Bool_t isbending, Int_t index,Int_t nstrip);
	virtual void ProcessRawHit(const TVector3 &xx, TObjArray &hits, Bool_t isbending);


  	//function to set/get sigmas for unknown coordinate
	
	void SetSigmaX (Double_t a) {fSigmaX = a;}
	void SetSigmaY (Double_t a) {fSigmaY = a;}
	void SetSigmaZ (Double_t a) {fSigmaZ = a;}
	inline Bool_t  IsBending() const {return fIsBending;}
	inline Bool_t  IsInUse()   const {return fIsInUse;}
    inline virtual void Set_InUseFlag() {fIsInUse=kTRUE;}
	Double_t GetSigmaX() const {return fSigmaX;}
	Double_t GetSigmaY() const {return fSigmaY;}
	Double_t GetSigmaZ() const {return fSigmaZ;}
	

   
 private:
	
  static Double_t fSigmaX;  // sigma_x
  static Double_t fSigmaZ;  // sigma_z
  static Double_t fSigmaY;  // sigma_y
  Bool_t fIsBending;	//flag to tell whether the layer is in bending/non-bending plane
  Bool_t fIsInUse;		//flag to tell whether the layer has been chosen in the detector cradle
    
  ClassDef(ALMeasLayer, 1) 	//Sample measurement layer class

};


#endif
	
	
