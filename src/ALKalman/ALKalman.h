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
#include "TNtupleD.h"		  // from ROOT
#include "TFile.h"			  // from ROOT
#include "TBField.h"
#include "TTree.h"
#include "ALEvent.h"

#include <iostream>
class ALKalman {

public:
ALKalman () {
   cout << "Object is being created" << endl;
}
	 virtual ~ALKalman() {};
	 void Reconstruct();
	
	ClassDef (ALKalman, 1)
};

#endif