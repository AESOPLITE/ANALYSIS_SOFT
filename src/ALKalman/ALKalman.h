
#ifndef __ALKALMAN__
#define __ALKALMAN__

#include "TROOT.h"
#include "TApplication.h"
#include "TRint.h"
#include "TRandom.h" 
#include "TFile.h"      
#include "TTree.h"  
#include "TKalDetCradle.h"    // from KalTrackLib
#include "TKalTrackState.h"   // from KalTrackLib
#include "TKalTrackSite.h"    // from KalTrackLib
#include "TKalTrack.h"        // from KalTrackLib
#include "ALMeasLayer.h"	  // from ALKalman
#include "ALKalDetector.h"	  // from ALKalman
#include "ALHit.h"			  // from ALKalman
#include "TFile.h"			  // from ROOT
#include "TBField.h"
#include "TTree.h"
#include "ALEvent.h"


#include <iostream>
#include "TNtupleD.h"		  // from ROOT
#include "TFile.h"			  // from ROOT
#include <iomanip>
#include <vector>
#include <algorithm>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <TROOT.h>
#include <Riostream.h>
#include "TChain.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLeaf.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLinearFitter.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TCut.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"
#include "TPDF.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TPostScript.h"
#include "TPaveText.h"
#include "TString.h"
#include "TH1F.h"
#include "TSystem.h"
#include "THStack.h"
#include "TNtuple.h"
#include "TFormula.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TH3F.h"
#include "TRandom3.h"
#include "TRandom2.h"
#include "TRandom1.h"
#include "TVector3.h"
#include "TSpline.h"
#include "Math/Interpolator.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include <iostream>
class ALKalman {

public:
    ALKalman(ALEvent *re) ; 
	virtual ~ALKalman();
	
	 int InitializeHelix(ALEvent *re, int InitType, bool secondIter, int DataType, TKalMatrix &state, TKalMatrix &covariant);	//InitType 0 = MCInit, 1 = PRInit, 2 = 3ptHelix	 
	 int InitialFit(TKalMatrix &svd_first,  TKalMatrix &C_first, TKalMatrix &svd_last, TKalMatrix &C_last, int type);
     int DoKF(ALEvent *re, int DataType, int InitType, bool secondIter);
	
	//class member

  TObjArray     *kalhits;    // hit buffer to hold original hits, include the backward hits
  TKalDetCradle *cradle;     // detctor system
  ALKalDetector *detector;   // detector
  ClassDef (ALKalman, 1)
		
};

#endif
