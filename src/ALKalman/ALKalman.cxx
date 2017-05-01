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
//* of a generated helix with measurement surfaces. The points are added to a hit buffer, kalhits.

//* (Requires)
//* (Provides)
//*     class LKalTest
//* (Update Recorded)
//*
//*   
//***************************************************************************

#include "ALKalman.h"		  // from ALKalman
#include <iomanip>
#include <iostream>


using namespace std;
static const double kMelectron = 0.5109989461e-3;		//mass electron in GeV
static const Bool_t kDir = kIterBackward;				//direction for the filter (starts from last hit)
//static const Bool_t kDir = kIterForward;				//first hit to last (does not work as well...)


void ALKalman::Reconstruct()
{


   // ===================================================================
   //  Prepare a detector
   // ===================================================================
	TObjArray kalhits; 
	TKalDetCradle	cradle;				
    ALKalDetector detector;
	cradle.Install(detector); 		// install detector into its cradle
	//cradle.SwitchOffMS();     	//debug: switch off multiple scattering
    //cradle.SwitchOffDEDX();		//debug: switch off dE/dX losses 
	TBField *bfieldmap = new TBField;
	//Bool_t magfieldmap = bfieldmap->SetMagField();
    Bool_t bApply2Iter = true;									//if initialize with first fit
    Int_t counter = 0;
	
	
  // ===================================================================
  //  Prepare a Root tree output
  // ===================================================================	 
	
	TFile hfile("KalRecon_30MeV_Smoothed.root","RECREATE","KalTest");
	TNtupleD *hTrackMonitor = new TNtupleD("Reconstruction", "", "ndf:chi2:cl:phi0:tanl:cpa:p0");
		
  // ===================================================================
  //  Read R00T MC input file
  // ===================================================================

	Bool_t flagMC = true;
    TFile*filein=new TFile("/home/sarah/AESOPLITE/ALanalysis-master/ALanalysis-master/RawEvents/RawEventMC_30MeV_uniform.root","READ");
    TTree *tree;
  	if(flagMC)tree = (TTree*)filein->Get("MC");
  
   ALEvent *e = new ALEvent;      			
   tree->SetBranchAddress("event",&e); 
   int nentries=tree->GetEntries();
   for (int i=0;i<nentries;i++) {				//loop over all events
   tree->GetEntry(i);
    int nnhits = (int)e->get_Nhits();
	       
		  if (nnhits == 7) {	       //until PatternRecognition ready
			  
		  Double_t E0 = e->get_EkMC() *1000 ;
		  Double_t X0 = e->get_X0MC();
	      Double_t Y0 = e->get_Y0MC();
		  Double_t Z0 = e->get_Z0MC();
		 // cout << "Event " << i << "E0 =" << E0 << " MeV,		X0=" << X0 << "  Y0=" << Y0 << "  Z0=" << Z0 << endl;
          for(int j=0;j<nnhits;j++) {
          Float_t X=((e->get_hits().at(j))->get_xin()+(e->get_hits().at(j))->get_xout())/2;
          Float_t Z=((e->get_hits().at(j))->get_zin()+(e->get_hits().at(j))->get_zout())/2;
          Float_t Y=((e->get_hits().at(j))->get_yin()+(e->get_hits().at(j))->get_youtn())/2;
          TVector3 xx;                      						
          xx.SetXYZ(X,Y,Z);
		  TVector3 xv;
		  xv.SetXYZ(Y,Z,X);
		  TVector3 bfield = TBField::GetGlobalBfield(xv);	
		  ALMeasLayer &ms = *static_cast<ALMeasLayer *>(cradle.At((j*2)));					//point to kActive layer
          Bool_t bending = ms.IsBending();
		  Bool_t active = ms.IsActive();
	      ms.ProcessHit(xx, kalhits, bending);
       			 }
      		 }
      
	 
    if ((kalhits.GetEntries() == 7)) {
		
    // ============================================================
    //  Do Kalman Filter in backward direction for the initial helix
    // ============================================================   
      
    THelicalTrack helstart;
    TKalMatrix C_start(kSdim,kSdim);							
    InitialBackwardFit(kalhits, helstart, C_start);
		
      // ---------------------------
      //  Create a dummy site: sited
      // ---------------------------
	  Int_t i1 = (kDir == kIterBackward) ? kalhits.GetEntries()-1 : 0;
      ALHit hitd = *dynamic_cast<ALHit *>(kalhits.At(i1));
      hitd(0,1) = 1.e6;   // give a huge error to x
      hitd(1,1) = 1.e6;   // give a huge error to y
      TKalTrackSite &sited = *new TKalTrackSite(hitd);
      sited.SetOwner();               //site own states

      
      // ---------------------------
      //  Set dummy state to sited
      // ---------------------------
      
      static TKalMatrix svd(kSdim,1);
      svd(0,0) = 0;
      svd(1,0) = helstart.GetPhi0();
      svd(2,0) = helstart.GetKappa();
      svd(3,0) = 0.;
      svd(4,0) = helstart.GetTanLambda();
      if (kSdim == 6) svd(5,0) = 0.;

    static TKalMatrix C(kSdim,kSdim);
    if(bApply2Iter){
       C = C_start;				 
    }
    else
    {
      for (Int_t i=0; i<kSdim; i++) {
		C(i,i) = 0.05;         // dummy error matrix
      }
	}
      

      sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kPredicted));
      sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kFiltered));

      // ---------------------------
      //  Add sited to the kaltrack
      // ---------------------------

      TKalTrack kaltrack;    // a track is a kal system
	  kaltrack.SetMass(kMelectron);	 //electrons are the reconstructed particles
      kaltrack.SetOwner();   // kaltrack owns sites
      kaltrack.Add(&sited);  // add the dummy site to the track
	
      
	TIter next(&kalhits, kDir);
    ALHit *hitp =  0;

       while ((hitp = dynamic_cast<ALHit *>(next())))  {
		   const ALMeasLayer &ml = dynamic_cast<const ALMeasLayer &>(hitp->GetMeasLayer());	
           TVector3 xv = ml.HitToXv(*hitp);
		   TVector3 xraw = hitp->GetRawXv();
		  TKalTrackSite  &site = *new TKalTrackSite(*hitp);   // create a site for this hit
		   		TVector3 bfield = TBField::GetGlobalBfield(xv);
         if (!kaltrack.AddAndFilter(site)) {  		          // add and filter this site
		//	cerr << " site discarded!" << endl;
         	delete &site;
      	} 
		   else {
			 TVKalState *state_fil = (TVKalState*) &(site.GetCurState());
	 		 THelicalTrack hel_fil = (dynamic_cast<TKalTrackState *>(state_fil))->GetHelix();
	  		 TVector3 x_fil=hel_fil.CalcXAt(0.0);

	  		TVKalState *state_exp = &site.GetState(TVKalSite::kPredicted);
	  		THelicalTrack hel_exp = (dynamic_cast<TKalTrackState *>(state_exp))->GetHelix();
	  		TVector3 x_exp=hel_exp.CalcXAt(0.0);
			Double_t mom = hel_fil.GetMomentum();
	 
			   
	//debug: comparing the expected and filtered points	   
/*
	   cout << "\t xraw =("<< xraw.X()<<",  "<<xraw.Y()<<", "<<xraw.Z()<<"): \n";
	   cout << "\t xv   =("<< xv.X()<<",  "<<xv.Y()<<", "<<xv.Z()<<"): \n";
	   cout << "\t x_exp=("<< x_exp.X()<<",  "<<x_exp.Y()<<", "<<x_exp.Z()<<"): \n";
	   cout << "\t x_fil=("<< x_fil.X()<<",  "<<x_fil.Y()<<", "<<x_fil.Z()<<"): \n";
       cout << "\t momentum = " << mom * 1000 << " MeV \n";
*/
			}   
       } 
      
	 kaltrack.SmoothBackTo(1);                          // smooth back to first site
	
      // ============================================================
      //  Monitor Fit Result
      // ============================================================

      Int_t    ndf  = kaltrack.GetNDF();
      Double_t chi2 = kaltrack.GetChi2();
      Double_t cl   = TMath::Prob(chi2, ndf);
      Double_t cpa  = kaltrack.GetCurSite().GetCurState()(2, 0);
      Double_t tanl  = kaltrack.GetCurSite().GetCurState()(4, 0);
      Double_t phi0  = kaltrack.GetCurSite().GetCurState()(1, 0);
	  Double_t rho =  kaltrack.GetCurSite().GetCurState()(0, 0);
	  Double_t pt = fabs(1.0/cpa);
      Double_t pz = pt * tanl;
  	  Double_t p0  = 1000 * pt * sqrt(1+tanl*tanl);   
      hTrackMonitor->Fill(ndf, chi2, cl, phi0, tanl, cpa, p0);
	  counter++;
    }
        kalhits.Delete();                       //clear buffer at the end of event
     }   									   //end of loop over events
	
  hfile.Write();
  filein->Close();
  Plot();
}

	
void ALKalman::InitialBackwardFit(TObjArray &kalhits, THelicalTrack &Hel_1st,TKalMatrix &C_1st) 
{
	      
   //loop over the entries in kalhit buffer
	       
      Int_t i1, i2, i3;
      if (kDir == kIterBackward) {
         i3 = 0;
         i1 = kalhits.GetEntries() - 1;
         i2 = i1 / 2;
      } else {
         i1 = 0;
         i3 = kalhits.GetEntries() - 1;
         i2 = i3 / 2;
      }                       
     
      // ---------------------------
      //  Create a dummy site: sited
      // ---------------------------
      ALHit hitd = *dynamic_cast<ALHit *>(kalhits.At(i1));
      hitd(0,1) = 1.e6;   // give a huge error to x
      hitd(1,1) = 1.e6;   // give a huge error to y
      TKalTrackSite &sited = *new TKalTrackSite(hitd);
      sited.SetOwner();               //site own states
      
      // ----------------------------------------
      // Create initial helix with MC full 3D hits 
      // ----------------------------------------
      ALHit   &h1 = *dynamic_cast<ALHit *>(kalhits.At(i1));   // first hit
      ALHit   &h2 = *dynamic_cast<ALHit *>(kalhits.At(i2));   // last hit
      ALHit   &h3 = *dynamic_cast<ALHit *>(kalhits.At(i3));   // middle hit
      TVector3 x1 = h1.GetMeasLayer().HitToXv(h1);
      TVector3 x2 = h2.GetMeasLayer().HitToXv(h2);
      TVector3 x3 = h3.GetMeasLayer().HitToXv(h3);
	  Double_t init_bfield1 = h1.GetBfield();				  //in kG
	  Double_t init_bfield2= h2.GetBfield();
	  Double_t init_bfield3 = h3.GetBfield();
     
      THelicalTrack helstart(x1, x2, x3, init_bfield1 , kDir); // initial helix 
      
      // ---------------------------
      //  Set dummy state to sited
      // ---------------------------
      
      static TKalMatrix svd(kSdim,1);
      svd(0,0) = 0.;
      svd(1,0) = helstart.GetPhi0();
      svd(2,0) = helstart.GetKappa();
      svd(3,0) = 0.;
      svd(4,0) = helstart.GetTanLambda();
      if (kSdim == 6) svd(5,0) = 0.;


     static TKalMatrix C(kSdim,kSdim);
      for (Int_t k=0; k<kSdim; k++) {
         C(k,k) = 0.005;   								// dummy error matrix
      }

      sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kPredicted));
      sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kFiltered));

      // ---------------------------
      //  Add sited to the kaltrack
      // ---------------------------

      TKalTrack kaltrack;    // a track is a kal system
	  kaltrack.SetMass(kMelectron);
      kaltrack.SetOwner();   // kaltrack owns sites
      kaltrack.Add(&sited);  // add the dummy site to the track

      
	TIter next(&kalhits, kDir);
    ALHit *hitp =  0;

       while ((hitp = dynamic_cast<ALHit *>(next())))  {
		   const ALMeasLayer &ml = dynamic_cast<const ALMeasLayer &>(hitp->GetMeasLayer());	
           TVector3 xv = ml.HitToXv(*hitp);
		   TVector3 xraw = hitp->GetRawXv();
		   TKalTrackSite  &site = *new TKalTrackSite(*hitp); // create a site for this hit
		   TVector3 bfield = TBField::GetGlobalBfield(xv);
         if (!kaltrack.AddAndFilter(site)) {  		         // add and filter this site
         	 
			 delete &site;
			}   
	   }
		 
	// kaltrack.SmoothBackTo(1);                          // smooth back to first site

  // ============================================================
  //  Get the last site then return it back
  // ============================================================
    
  TKalTrackState *theLastState = dynamic_cast<TKalTrackState*> (&(kaltrack.GetCurSite().GetCurState()));
   Hel_1st = theLastState->GetHelix();
   C_1st = theLastState->GetCovMat();
   Double_t mom = Hel_1st.GetMomentum();
 // cout << "Momentum from fist filter iteration p = " << mom * 1000 << " MeV" << endl; 
}


void ALKalman::Plot()
{
	
	  // ============================================================
      //  Plot Fit Result
      // ============================================================ 
	
   TFile*file = new TFile("KalRecon_30MeV_Smoothed.root","READ");
   TNtupleD *ntuple = (TNtupleD*)file->Get("Reconstruction");
  TH1D *p0hist = new TH1D("p0hist", "Reconstructed momentum 30MeV electrons", 100, 0.0, 100);

   TCanvas* c1 = new TCanvas("c1", "reconstruction", 1);
   c1->Divide(2,2);
   c1->cd(1);
   p0hist->GetXaxis()->SetTitle("p0 (in Mev/c)");
   p0hist->GetXaxis()->SetRangeUser(0, 100);
 ntuple->Draw("p0>>p0hist");
    

  c1->cd(2);
  ntuple->Draw("phi0");
  c1->cd(3);
  ntuple->Draw("tanl");
  c1->cd(4);
  ntuple->Draw("cpa");

  c1->Update();
  c1->Print("Recontruction_30MeV.pdf");
	
}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	