
//*************************************************************************
//* ===================
//*  ALKalTest
//* ===================
//*
//* (Description)
//* 13/01/2017  S. Mechbal		Modified for AESOPLITE experiment
//*
//* This is the main program that executes the Kalman Filter. The first step is to load the ROOT 
//* output tree to TObjArray, create ALHit structure, initialize helix with three points and
//* do Kalman Filter. 
//* In KalTest, they have an additional class EXEventGen, which creates hit points from the intersection
//* of a generated helix with measurement surfaces. The points are added to a hit buffer, kalhits->

//* (Requires)
//* (Provides)
//*     class ALKalTest
//* (Update Recorded)
//*
//*   
//***************************************************************************

#include "ALKalman.h"		  // from ALKalman
#include "TTrackFrame.h"
#include <iomanip>
#include <iostream>
#include <cmath>


using namespace std;
static const double kMelectron = 0.5109989461e-3;		//electron mass in GeV
static const double kMuon = 0.105658371;				//muon mass in GeV
static const double kProton = 0.93827231;				//proton mass in GeV
static const double kAlpha = 3.727379378;				//alpha mass in GeV
static const double RestMassE   = 0.51099;				//mass electron in MeV
static const double RestMassMu = 105.6584;				//mass muon in MeV
static const double RestMassProton   = 938.2723;		//mass proton in MeV
static const double RestMassAlpha = 3727.37938;			//mass alpha in MeV
static const double c = TMath::C()*1000;				//speed of light in mm/s  
static const double average_B = 0.3;					//average magnetic field AESOPLITE	in [T]
static const double CovMElement=1e4;				//initial covariant matrix elements
static const double chi2lim= 1e10;					//limit on chi2 value of first iteration of filter
static const Bool_t kDir = kIterForward;				//first hit to last (does not work as well...)
static const Bool_t kDirFirst= kIterForward;
int uhitnid[7]={0,0,0,0,0,0,0};
bool testlayer[7]={false,false,false,false,false,false,false};
//Class constructor 
ALKalman::ALKalman(ALEvent *re) 
{
		   //cout <<endl;
   //cout << "Calling class constructor ALKalman() " << endl;

   // ===================================================================
   //  Prepare a detector
   // ===================================================================
   //ALEvent re;
   cradle = new TKalDetCradle();				
   detector = new ALKalDetector();
   cradle->Install(*detector); 	
   int nlayers = cradle->GetEntries();
//   cout << "Number of detectors in cradle = " << nlayers << endl;
     kalhits = new TObjArray(7);
	

   int nnhits = (int)re->get_Nhits();
   cout << "Event " << re->get_eventnumber() << "   number of hits = " << nnhits << endl;
   //Index of the hits used to reconstruct the track

   int ij = 0 ;
   int index;					//index of hits used
   int count = 0;
	
   //Use information of PatternRecognition to choose hits

   for(int j=0;j<nnhits;j++) 
  	 {
	      bool flagPR = re->get_hits().at(j)->get_flagPR();
	      bool fGhost = re->get_hits().at(j)->get_fGhost();
	      int k = re->get_hits().at(j)->get_k();
	    if(fGhost) continue;    
	    if(!flagPR) continue;	
	    count++;
	    int L = re->get_hits().at(j)->get_L();
	    int nstrip = re->get_hits().at(j)->get_nstrips();
		//coordinates fom MC or data
	    Float_t x=(re->get_hits().at(j)->get_x())*10;			//in mm
            Float_t y=(re->get_hits().at(j)->get_y())*10;           //in mm
            Float_t z=(re->get_hits().at(j)->get_z())*10;			//in mm
		//coordinate from PR fit
	    Float_t xPR=(re->get_hits().at(j)->get_xPR())*10;		//in mm
            Float_t yPR=(re->get_hits().at(j)->get_yPR())*10;       //in mm		
            Float_t zPR=(re->get_hits().at(j)->get_zPR())*10;		//in mm	
	   // cout << "L = " << L << ", x = " << x << ", y = " << y << " , z = " << z << endl;
	    TVector3 xx; 	    	                     							
	    if(L==0||L==4||L==6)//non-bending plane 
		{xx.SetXYZ(x,yPR,z);}
	    else 
		{xx.SetXYZ(xPR,y,z);}  	
	
		//loop over layer to choose active ones over dummy ones to record hit in the right layer
		for(int n=0; n<nlayers;n++) {
			ALMeasLayer &ms = *static_cast<ALMeasLayer *>(cradle->At(n));
			Bool_t isactive = ms.IsActive();
			Bool_t bending = ms.IsBending();
			Bool_t inuse	= ms.IsInUse();
			TVector3 Xc = ms.GetXc();
			if((isactive) && (!inuse)) {
				ms.Set_InUseFlag();
				ms.ProcessHit(xx, *kalhits, bending, L,nstrip);
		//		cout << "n= " << n <<", hit added to L = "<<  L << ", bending = " << bending << endl;
				testlayer[ij] = true;
			//We record the index of the hit so we will be able to record the reconstructed variables at the right place 
				uhitnid[ij]=L;
				ij++;
				break;
				}
	    	} //end loop over layers
	 
  	} //end j	
 	
     //if we have found one hit in each layer then we stop looking for more
   if(testlayer[0]&&testlayer[1]&&testlayer[2]&&testlayer[3]&&testlayer[4]&&testlayer[5]&&testlayer[6]) //cout << "retrieved the hits" <<endl;
     if(testlayer[0]&&testlayer[1]&&testlayer[2]&&testlayer[3]&&testlayer[4]&&testlayer[5]&&testlayer[6]==false) 
      {
       cout << "We did not find hits on all layer!!" <<endl;
      }
	if(count > 7) {
		cout << "Too many hits, PatternReco messed up! " << endl;
	}
   //cout << "Indexes of the hits used in the reconstruction" << endl;
   for(int i=0;i<7;i++) cout << uhitnid[i] << " "; //Check the layer
   //cout << endl;
//  cout << "number of entries in kalhit buffer is " << kalhits->GetEntries() << endl;
 
 for (int ind = 0; ind < kalhits->GetEntries(); ind++) 
	{
	       index = uhitnid[ind];	
	       ALHit hitd = *dynamic_cast<ALHit *>(kalhits->At(index));
	       TVector3 Pos = hitd.GetRawXv();
//	       cout << "kalhit index = " << ind << ", index = " <<  Pos.Y() << endl;
	}		



}

//class destructor
ALKalman::~ALKalman()
{
  delete kalhits;
  delete cradle;
  delete detector;

}

int ALKalman::InitializeHelix(ALEvent *re, int InitType, bool secondIter, int DataType, TKalMatrix &state, TKalMatrix &covariant)	
{
	//Choose 3 hits for initial 3-points helix
	cout << "InitializeHelix() " << endl;
	double kappa_init, phi0_init, tanL_init;
	kappa_init=phi0_init=tanL_init=0;
	double kappa_MC, phi0_MC, tanL_MC;
	kappa_MC = phi0_MC=tanL_MC=0;
	double kappa_PR, phi0_PR, tanL_PR;
	kappa_PR = phi0_PR=tanL_PR=0;
	double kappa_3, phi0_3,  tanL_3;
	kappa_3=phi0_3=tanL_3=0;
	int type;
	double deflec = re->get_deflecPR();
	double slope_PR = re->get_slopePR();
	double Q =  TMath::Sign(1,deflec);
	//cout << "slope from PR = " <<  slope_PR <<endl;
	   if (DataType==1) { 
	   if (re->get_T2()) {			//if CK fired, electrons or positrons
	      if (Q>0) type=4;
	      else type=3;
	      }
	    else   {			//if CK not fired
	      if (Q>0) type=10;
	      else type=11;	       			
	    }
	  } //if data
	  //if MC	
	   else type = re->get_typeMC();	//if MC				   
		   
        static TKalMatrix svd_first(kSdim,1);
         svd_first(0,0) = 0.;
         svd_first(3,0) = 0.;
         if (kSdim == 6) svd_first(5,0) = 0.;
	 int lTop = 99;
	 for(int k=0;k<re->get_Nhits();k++)
            {
             int Lindex=(int)re->get_hits().at(k)->get_L();
             bool flagPR = re->get_hits().at(k)->get_flagPR();
             if (!flagPR) continue;
             if (Lindex < lTop) lTop=Lindex;
			}	
           //cout << "check, the topmost layer is " << lTop << endl;  
	   for(int k=0;k<re->get_Nhits();k++) 
	    {
	     int Lindex=(int)re->get_hits().at(k)->get_L(); 		
             bool flagPR = re->get_hits().at(k)->get_flagPR();
	     if (!flagPR) continue;	   
	     if(Lindex==lTop) {
		cout  << "The topmost layer is " << lTop << endl;								
	       //Convert from FLUKA to KALTEST coordinates!!
		double cxL0 = re->get_hits().at(k)->get_cy();
		double cyL0 = re->get_hits().at(k)->get_cz();
		double czL0 = re->get_hits().at(k)->get_cx();
		double eMC =(re->get_hits().at(k)->get_eMC());       	
		TVector3 p0(eMC*cxL0, eMC*cyL0, eMC*czL0);                             		
		phi0_MC = TMath::ATan2(-p0.X(), p0.Y());
		kappa_MC = Q/TMath::Sqrt(p0.X()*p0.X() + p0.Y()*p0.Y());
		tanL_MC = p0.Z()/TMath::Sqrt(p0.X()*p0.X() + p0.Y()*p0.Y());
		
	//get PR variables							
	       //Convert from FLUKA to KALTEST coordinates!!
		double cxPR = re->get_hits().at(k)->get_cyPR();
		double cyPR = re->get_hits().at(k)->get_czPR();
		double czPR = re->get_hits().at(k)->get_cxPR();
		double EkPR = re->get_EkPR();   
      // -----------------------------------------------
      // Create initial helix using PR parameters at L0
      // -----------------------------------------------  	   
		
	   TVector3 p0PR(EkPR*cxPR, EkPR*cyPR, EkPR*czPR);                               
	   phi0_PR = TMath::ATan2(-p0PR.X(), p0PR.Y());
	   kappa_PR = Q/TMath::Sqrt(p0PR.X()*p0PR.X() + p0PR.Y()*p0PR.Y());
	   tanL_PR = p0PR.Z()/TMath::Sqrt(p0PR.X()*p0PR.X() + p0PR.Y()*p0PR.Y());
		//tanL_PR=fabs(kappa_PR)*p0PR.Z();

		 } //if L0
	}	//end for k

        int i1, i2, i3;
        if (kDir == kIterBackward) {
         i3 = 0;
         i1 = kalhits->GetEntries() - 1;
         i2 = i1 / 2;
		}
         else {
         i1 = 0;
         i3 = kalhits->GetEntries() - 1;
         i2 = i3 / 2;
      	 }
		 // ----------------------------------------
		 // Create initial helix using 3 hits
		 // ----------------------------------------
/*
		 ALHit   &h1 = *dynamic_cast<ALHit *>(kalhits->At(i1));   // first hit
		 ALHit   &h2 = *dynamic_cast<ALHit *>(kalhits->At(i2));   // last hit
		 ALHit   &h3 = *dynamic_cast<ALHit *>(kalhits->At(i3));   // middle hit
		 TVector3 x1 = h1.GetMeasLayer().HitToXv(h1);
		 TVector3 x2 = h2.GetMeasLayer().HitToXv(h2);
		 TVector3 x3 = h3.GetMeasLayer().HitToXv(h3);
		 Double_t init_bfield1 = h1.GetBfield();				  //in T
		 Double_t init_bfield2= h2.GetBfield();
		 Double_t init_bfield3 = h3.GetBfield();
		 THelicalTrack helstart(x1, x2, x3, init_bfield2 , kDir); // initial helix
	    	 phi0_3 = helstart.GetPhi0();
	     	 kappa_3 = helstart.GetKappa();
	     	 tanL_3 = helstart.GetTanLambda();
*/	     
	     if(InitType==0) {				//MC init
	     	kappa_init=kappa_MC;
	     	tanL_init=tanL_MC;
	     	phi0_init=phi0_MC;
	     	}
	     else if(InitType==1) {			//PR init
	     	kappa_init=kappa_PR;
	     	tanL_init=tanL_PR;
	     	phi0_init=phi0_PR;
		}

	//     else {
	//     phi0_init=phi0_3;
	//     kappa_init=kappa_3;
	//     tanL_init=tanL_3;
//		}
        

	  svd_first(1,0) = phi0_init;
	  svd_first(2,0) = kappa_init;
          svd_first(4,0) = tanL_init;
	  //cout << "phi_init = " << phi0_init << "  cpa_init = " << kappa_init << "  tanL = " << tanL_init << endl;
	  //cout << "phi_MC = " << phi0_MC << "  cpa_MC = " << kappa_MC << "  tanL_MC = " << tanL_MC << endl;
	  //cout << "phi_PR = " << phi0_PR << "  cpa_PR = " << kappa_PR << "  tanL_PR = " << tanL_PR << endl;
	  //cout << "phi_3 = " << phi0_3 << "  cpa_3 = " << kappa_3 << "  tanL_3 = " << tanL_3 << endl;
      static TKalMatrix C_first(kSdim,kSdim);
       for (Int_t k=0; k<kSdim; k++) {
		   C_first(k,k) = CovMElement;  //huge error matrix to start with
	   }   								
      	
    

//for two iterations
	if(secondIter) {
	
   //  Do Kalman Filter with initial helix, return state vector and covariant matrix at last site
  
   TKalMatrix svd_last(kSdim,1);
   TKalMatrix C_last(kSdim,kSdim);
   int initfit = InitialFit(svd_first, C_first, svd_last, C_last, type);
  
   state = svd_last;
   covariant = C_last;
   if(initfit == 0) {
        //cout << "initfit set to 0, exiting InitializeHelix" << endl;
	return 0;
}
	}  //end if 2nd iteration

//if only one iteration
 else {
	 
    state = svd_first;
    covariant =  C_first;
     //save initial helix parameters

 }
	 
//save initial state vector and covariant matrix
	re->set_phi0_init(state(1,0));
	re->set_cpa_init(state(2,0));
	re->set_tanl_init(state(4,0));
      //  re->set_Cov_init(covariant);
      //	state.DebugPrint("initial state vector");
    //covariant.DebugPrint("initial cov matrix");
     //cout << "end initialize helix function " << endl;
     return type;

}	
	
	 

 int ALKalman::InitialFit(TKalMatrix &svd_first, TKalMatrix &C_first, TKalMatrix &svd_last,TKalMatrix &C_last, int type)	 
{
	cout << "InitialFit called, doing first iteration KF " << endl;  
	
// //cout << " number of entries kalhit buffer = " << kalhits->GetEntries() << endl;
	kalhits->SetOwner(false);						

      
      // ---------------------------
      //  Create a dummy site: sited
      // ---------------------------
      int i1 = (kDirFirst == kIterForward) ? 0 : kalhits->GetEntries()-1;
      ALHit hitd = *dynamic_cast<ALHit *>(kalhits->At(i1));
      hitd(0,1) = 1.e6;   // give a huge error to x
      hitd(1,1) = 1.e6;   // give a huge error to y
      TKalTrackSite &sited = *new TKalTrackSite(hitd);
      sited.SetOwner();               //site own states   

      // ---------------------------
      //  Set dummy state to sited
      // ---------------------------
      

      // svd_first.DebugPrint("initial state vector");
      // C_first.DebugPrint("initial covariant matrix");
      sited.Add(new TKalTrackState(svd_first,C_first,sited,TVKalSite::kPredicted));
      sited.Add(new TKalTrackState(svd_first,C_first,sited,TVKalSite::kFiltered));

      // ---------------------------
      //  Add sited to the kaltrack
      // ---------------------------

      TKalTrack kalfirst;    // a track is a kal system
      if(type==3 || type==4)  kalfirst.SetMass(kMelectron);
      else if (type == 10 || type == 11) kalfirst.SetMass(kMuon);
      else if (type == 1) kalfirst.SetMass(kProton);
	  else if (type == -6) kalfirst.SetMass(kAlpha);
      kalfirst.SetOwner();   // kaltrack owns sites
      kalfirst.Add(&sited);  // add the dummy site to the track

      
      TIter first(kalhits, kDirFirst);
      ALHit *hitp =  0;

      while ((hitp = dynamic_cast<ALHit *>(first())))
        {
          const ALMeasLayer &ml = dynamic_cast<const ALMeasLayer &>(hitp->GetMeasLayer());	
          TVector3 xv = ml.HitToXv(*hitp);
          TKalTrackSite  &site = *new TKalTrackSite(*hitp); // create a site for this hit
          TVector3 bfield = TBField::GetGlobalBfield(xv);
          if (!kalfirst.AddAndFilter(site))
           { // add and filter this site
            delete &site;
           }   
         }
   
   kalfirst.SmoothBackTo(1);                          // smooth back to first site

  // ============================================================
  //  Get the last site then return it back
  // ============================================================
    
   TKalTrackState *theLastState = dynamic_cast<TKalTrackState*> (&(kalfirst.GetCurSite().GetCurState()));
   THelicalTrack Hel_last = theLastState->GetHelix();
   svd_last(0,0) = 0.0;
   svd_last(1,0) = Hel_last.GetPhi0();
   svd_last(2,0) = Hel_last.GetKappa();
   svd_last(4,0) = Hel_last.GetTanLambda();
   svd_last(5,0) = 0.0;
   if (kSdim == 6) svd_last(5,0) = 0.;
   svd_last.DebugPrint("state vector from 1st KF iteration");
   C_last = theLastState->GetCovMat();
   Double_t mom = Hel_last.GetMomentum();
  cout << "Momentum from fist filter iteration p = " << mom * 1000 << " MeV" << endl; 
   double chi2first = kalfirst.GetChi2();
   if(fabs(chi2first) > chi2lim) {
    cout << "Initial fit chi2 too high, throwing event away" << endl;	
    return 0;
	}
   else return 1;
}




int ALKalman::DoKF(ALEvent *re, int DataType, int InitType, bool secondIter)
{	
	
  cout << "DoKF() function called" << endl; 
  static TKalMatrix state(kSdim,1);
  static TKalMatrix covariant(kSdim,kSdim);
  int type = InitializeHelix(re, InitType, secondIter, DataType, state, covariant);
  if(type==0) {
	//cout << " exiting DoKF function" << endl;
	return 0;
	} 	

		
      // ---------------------------
      //  Create a dummy site: sited
      // ---------------------------
      int i1 = (kDir == kIterForward) ? uhitnid[0] : uhitnid[6];
     
      ALHit hitd = *dynamic_cast<ALHit *>(kalhits->At(i1));
      hitd(0,1) = 1.e6;   // give a huge error to x
      hitd(1,1) = 1.e6;   // give a huge error to y
      TKalTrackSite &sited = *new TKalTrackSite(hitd);
      sited.SetOwner();               //site own states
		
	  // ---------------------------
      //  Set dummy state to sited
      // ---------------------------
         sited.Add(new TKalTrackState(state,covariant,sited,TVKalSite::kPredicted));
         sited.Add(new TKalTrackState(state,covariant,sited,TVKalSite::kFiltered));

      // ---------------------------
      //  Add sited to the kaltrack
      // ---------------------------

      TKalTrack kaltrack;    // a track is a kal system
      if(type==3 || type==4)  kaltrack.SetMass(kMelectron);
      else if (type == 10 || type == 11) kaltrack.SetMass(kMuon);
      else if (type == 1) kaltrack.SetMass(kProton);
	  else if (type == -6) kaltrack.SetMass(kAlpha);
      kaltrack.SetOwner();   // kaltrack owns sites
      kaltrack.Add(&sited);  // add the dummy site to the track


TIter next(kalhits, kDir);
ALHit *hitp = 0;
int usedSites=0;
int hit_index = (kDir == kIterForward) ? uhitnid[0] : uhitnid[6];
      
 while ((hitp = dynamic_cast<ALHit *>(next()))) {
      //cout << "-----------------------Next Hit--------------------------- " << endl;
      const ALMeasLayer &ml = dynamic_cast<const ALMeasLayer &>(hitp->GetMeasLayer());	
      Bool_t isactive = ml.IsActive();
      Bool_t bend = ml.IsBending();
      TVector3 xv = ml.HitToXv(*hitp);
      TVector3 xraw = hitp->GetRawXv();
      TKalTrackSite  &site = *new TKalTrackSite(*hitp); // create a site for this hit
      TVector3 TVbfield = TBField::GetGlobalBfield(xv);
      double B = TVbfield.Mag();
     //cout << "x = "  << xv.X() << " y ="  << xv.Y() << " z ="  << xv.Z() << endl;
     //cout << "B = (" << TVbfield.X() << "," << TVbfield.Y() << "," << TVbfield.Z() << ")"  << endl;	
//		  //cout << "Phi = "  <<  xv.Phi() << " Perp = " << xv.Perp() << endl;
	  		 	 
      if (!kaltrack.AddAndFilter(site))
       { // add and filter this site
        //site.DebugPrint();
      // kaltrack.GetState(TVKalSite::kFiltered).DebugPrint();
//        cerr << " site discarded!" << endl;
        delete &site;
       } 
/////////////////////////////////IF SITE ACCEPTED/////////////////////////////////
      else 
       {
        
	TVKalState *state_fil = (TVKalState*) &(site.GetCurState());
        THelicalTrack hel_fil = (dynamic_cast<TKalTrackState *>(state_fil))->GetHelix();
        TVector3 x_fil=hel_fil.CalcXAt(0.0);
        Double_t mom = hel_fil.GetMomentum();	
        TVector3 x_glob = site.GetGlobalX();
	//TVector3 loctoloc = 
 
		   
		  //smooth back

	
        //   TVector3 x_smoothed = kaltrack.SmoothBackLastSite(hit_index+1);                  
	 	//  TVKalState &state_smoothed = site.GetState(TVKalSite::kSmoothed);
                //  THelicalTrack hel_smoothed = (dynamic_cast<TKalTrackState &>(state_smoothed)).GetHelix();
                  //TVector3 x_smoothed=hel_smoothed.CalcXAt(0.0);  
		  //cout << "\t xraw   =("<< xraw.X()<<",  "<<xraw.Y()<<", "<<xraw.Z()<<"): \n";
                  //cout << "\t x_fil=("<< x_fil.X()<<",  "<<x_fil.Y()<<", "<<x_fil.Z()<<"): \n";
                //  //cout << "\t x_smoothed=("<< x_smoothed.X()<<",  "<<x_smoothed.Y()<<", "<<x_smoothed.Z()<<"): \n";
	          //cout << "\t x_glob=("<< x_glob.X()<<",  "<<x_glob.Y()<<", "<<x_glob.Z()<<"): \n";	
		   
    
   // -------------------------------------------------------------------------------
   //  Add correspondance between ALTckhit structure and ALHit, only if active layer
   // -------------------------------------------------------------------------------
 if(isactive) {
	 float xreco, yreco, zreco, ereco;
	 if(TBField::IsUsingUniformBfield()) {
//(filtered) reconstructed position
	//!!!!!!  CONVERT BACK TO FLUKA COORDINATES !!!!!!!
	 xreco = x_fil.Z();
     yreco = x_fil.X();
	 zreco = x_fil.Y();
	 }
	 else {
	 xreco = x_glob.Z();
     	 yreco = x_glob.X();
	 zreco = x_glob.Y();	
	 }
		
        //directional cosines, find tangent to momentum vector
	double deflec = re->get_deflecPR();
	double Q =  TMath::Sign(1,deflec);
        TVector3 pivot = hel_fil.CalcXAt(0.0);
        TMatrixD dxdphi = hel_fil.CalcDxDphi(0.0);
	float alpha = 1/(c*B);				// in s/mm/T					
        TVector3 momentum((-Q/alpha)*dxdphi(0,0),(-Q/alpha)*dxdphi(1,0),(-Q/alpha)*dxdphi(2,0));	        // convert matrix diagonal to momentum vector
	TVector3 vtan = momentum.Unit();
	float mag = vtan.Mag();
	float theta = vtan.Theta();
        float phi   = vtan.Phi();
        //!!!!!!  CONVERT BACK TO FLUKA COORDINATES !!!!!!!
        float cxreco = (vtan.Z())/mag;
	float cyreco = (vtan.X())/mag;
	float czreco = (vtan.Y())/mag;

	float tot=(cxreco*cxreco) + (cyreco*cyreco) + (czreco*czreco);
//	//cout << "tot = " << tot << "  cxreco = " << cxreco << "  cyreco = " << cyreco << "  czreco = " << czreco << endl;
//	//cout << "cxMC = " << re->get_hits().at(hit_index)->get_cx() << " cyMC = " << re->get_hits().at(hit_index)->get_cy() << "  czMC = " << re->get_hits().at(hit_index)->get_cz() <<endl;
//		//cout << "cxPR = " << re->get_hits().at(hit_index)->get_cxPR() << " cyPR = " << re->get_hits().at(hit_index)->get_cyPR() << "  czPR = " << re->get_hits().at(hit_index)->get_czPR() <<endl;
		   
	//kinetic energy of particle   
	    if(type==3 || type==4)  ereco = sqrt(mom*mom + (RestMassE)*(RestMassE))  - RestMassE;   //for electrons
	  else if (type ==10 || type == 11) ereco = sqrt(mom*mom + (RestMassMu)*(RestMassMu))  - RestMassMu;   //for muons
	  else if (type ==1) ereco = sqrt(mom*mom + (RestMassProton)*(RestMassProton))  - RestMassProton;   //for protons
	  else if (type ==-6) ereco = sqrt(mom*mom + (RestMassAlpha)*(RestMassAlpha))  - RestMassAlpha;   //for alphas
	//fill variables 
    
	re->set_hxreco(uhitnid[hit_index],xreco);
	re->set_hyreco(uhitnid[hit_index],yreco);
	re->set_hzreco(uhitnid[hit_index],zreco);
	re->set_hcxreco(uhitnid[hit_index],cxreco);
	re->set_hcyreco(uhitnid[hit_index],cyreco);
	re->set_hczreco(uhitnid[hit_index],czreco);
	re->set_hereco(uhitnid[hit_index], ereco);
	re->get_hits().at(uhitnid[hit_index])->set_fUsed(true);
	usedSites++;	
	//cout << "reconstructed info added hit_index = " << hit_index << " at index " << uhitnid[hit_index] << endl;
	//cout << "zreco = " << zreco << "  fUsed bool set to true" << endl;

 		} //end if active
 }//end else  

	    if (kDir == kIterBackward) hit_index--;
	else hit_index++;	   
      } //end while
	
     //Smooth Back
//int isite = 1 ;
 kaltrack.SmoothBackTo(1);
     TVKalSite &cursite = static_cast<TVKalSite &>(*kaltrack[1]);

   // =======================================================================
   //  Find tangent to helix at first site and extrapolate injection point
   // =======================================================================   
   //cout << "Find tangent to helix at site and extrapolate injection point" <<endl;
   TVKalState *state_first =(TVKalState*) &(cursite.GetCurState());
   THelicalTrack hel_first = (dynamic_cast<TKalTrackState *>(state_first))->GetHelix();
   TVector3 pivot = hel_first.CalcXAt(0.0);
   TMatrixD dxdphi = hel_first.CalcDxDphi(0.0);                                 // tangent vector at destination surface
   TVector3 vtan(dxdphi(0,0),dxdphi(1,0),dxdphi(2,0));    
/*
   for(int isite=0;isite<7;isite++) {
     //cout << "site " << isite << endl;       
     TVKalSite &cursite = static_cast<TVKalSite &>(*kaltrack[isite]);
 	  
   // =======================================================================
   //  Find tangent to helix at first site and extrapolate injection point
   // =======================================================================   
   //cout << "Find tangent to helix at site and extrapolate injection point" <<endl;
   TVKalState *state_first =(TVKalState*) &(cursite.GetCurState());
   THelicalTrack hel_first = (dynamic_cast<TKalTrackState *>(state_first))->GetHelix();
   TVector3 pivot = hel_first.CalcXAt(0.0);
   TMatrixD dxdphi = hel_first.CalcDxDphi(0.0);					// tangent vector at destination surface
   TVector3 vtan(dxdphi(0,0),dxdphi(1,0),dxdphi(2,0));				// convert matrix diagonal to vector
  // //cout << "Tangent vector vtan=("<< vtan.X()<<",  "<<vtan.Y()<<", "<<vtan.Z()<<"): \n";
  
      double deflec = re->get_deflecPR();
        double Q =  TMath::Sign(1,deflec);
        double B = 0.3;
        float alpha = 1/(c*B);                          // in s/mm/T
        TVector3 momentum((-Q/alpha)*dxdphi(0,0),(-Q/alpha)*dxdphi(1,0),(-Q/alpha)*dxdphi(2,0));                // convert matrix diagonal to $
        TVector3 vtan = momentum.Unit();
        float mag = vtan.Mag();
        float theta = vtan.Theta();
        float phi   = vtan.Phi();
        //!!!!!!  CONVERT BACK TO FLUKA COORDINATES !!!!!!!
        float cxreco = (vtan.Z())/mag;
        float cyreco = (vtan.X())/mag;
        float czreco = (vtan.Y())/mag;
//cout << "site " << isite << "cx = " << cxreco << "   cy = " << cyreco << " cz = " << czreco << endl;

}

*/
   //Parametric equation of tangent line
   Double_t yinjection = 35; 
   Double_t t = (yinjection - pivot.Y())/(vtan.Y());
   //!!!!!!  CONVERT BACK TO FLUKA COORDINATES !!!!!!!
   double X0reco = pivot.Z() + vtan.Z() * t;
   double Y0reco = pivot.X() + vtan.X() * t;
   double Z0reco = pivot.Y() + vtan.Y() * t;
  
   //Get directional cosines of tangent line
   Double_t theta = vtan.Theta();
   Double_t phi   = vtan.Phi();
   Double_t CX0reco = sin(theta) * cos(phi);
   Double_t CY0reco = sin(theta) * sin(phi);
   Double_t CZ0reco = cos(theta);

   // ============================================================
   //  Fill reconstructed variables
   // ============================================================
  
    
   re->set_X0reco(X0reco); 
   re->set_Y0reco(Y0reco); 
   re->set_Z0reco(Z0reco); 
   re->set_CX0reco(CX0reco);
   re->set_CY0reco(CY0reco);
   re->set_CZ0reco(CZ0reco);
   int    ndf  = kaltrack.GetNDF();
   re->set_ndf(ndf);
   double chi2 = kaltrack.GetChi2();
   re->set_chi2(chi2);
   double cl   = TMath::Prob(chi2, ndf);
   re->set_cl(cl);
   double cpa  = cursite.GetCurState()(2, 0);
   re->set_cpa(cpa);
   double tanl  = cursite.GetCurState()(4, 0);
   re->set_tanl(tanl);
   double phi0  = cursite.GetCurState()(1, 0);
   re->set_phi0(phi0);
   double rho =  cursite.GetCurState()(0, 0);
   re->set_d0(rho);
   double dz =  cursite.GetCurState()(5, 0);
   re->set_dz(dz);
   //Cov_final = kaltrack.GetCurSite().GetCurState().GetCovMat();
   TVKalState *theLastState = (TVKalState*) &(cursite.GetCurState());
   const TKalMatrix Cov_final = theLastState->GetCovMat();
   double rhoerr2 = Cov_final(0,0);
   re->set_d0err2(rhoerr2);
   double phi0err2 = Cov_final(1,1);
   re->set_phi0err2(phi0err2);
   double cpaerr2 = Cov_final(2,2);
   re->set_cpaerr2(cpaerr2);
   double dzerr2 = Cov_final(5,5);
   re->set_dzerr2(dzerr2);
   double tanlerr2 = Cov_final(4,4);
   re->set_tanlerr2(tanlerr2);
   //cout << "phi0last = " << phi0 << "  cpalast = " << cpa << "  tanl_last " << tanl << endl;
   //cout << "phi0err2 = " << phi0err2 << "  cpaerr2 = " << cpaerr2 << "  tanlerr2 = " << tanlerr2 << endl;
  // re->set_Cov_last(Cov_final);
   double pt = 0;
   if(cpa!=0)pt=fabs(1.0/cpa);
   double pz = pt * tanl;
   double p0  = 1000 * pt * sqrt(1+tanl*tanl);   //p = pt / sinTheta in MeV	
   re->set_p0reco(p0);
   double Ekreco;
   if(type==3 || type ==4 ) Ekreco =  sqrt(p0*p0 + (RestMassE)*(RestMassE))  - RestMassE;			//for electrons
   else if(type==10 || type ==11) Ekreco =  sqrt(p0*p0 + (RestMassMu)*(RestMassMu))  - RestMassMu;	//for muons
   else if(type==1) Ekreco =  sqrt(p0*p0 + (RestMassProton)*(RestMassProton))  - RestMassProton;	//for protons
   else if(type==-6) Ekreco =  sqrt(p0*p0 + (RestMassAlpha)*(RestMassAlpha))  - RestMassAlpha;		//for alphas


   re->set_Ekreco(Ekreco);
   
   cout << "preco = " << p0 << " MeV" << endl; 
   cout << "chi2   =  " << chi2 << endl;
   cout <<  "ndf   =  " << ndf << endl;
   cout << "used sites = " <<         usedSites << endl;
   //cout << "End of MakeRecoEventMC" <<endl;
   next.Reset(); 
	return 1;
 }







	



