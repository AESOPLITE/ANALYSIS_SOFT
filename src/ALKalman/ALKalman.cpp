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
//*     class ALKalTest
//* (Update Recorded)
//*
//*   
//***************************************************************************

#include "ALKalman.h"		  // from ALKalman
#include <iomanip>
#include <iostream>


using namespace std;

static const Bool_t kDir = kIterBackward;	//direction for the filter (starts from last hit)
//static const Bool_t kDir = kIterForward;


//main program

void ALKalman::Reconstruct()
{
	
	//fFile = new TFile("h.root","RECREATE","Kalman Filter for AESOPLITE MC track");
 // fTree = new TTree("t", "Kalman Filter for AESOPLITE MC track");

   // ===================================================================
   //  Prepare a detector
   // ===================================================================
	
	
    TObjArray kalhits; 
	TKalDetCradle	cradle;				
    ALKalDetector detector;
	cradle.Install(detector); // install detector into its cradle
   TBField *bfield = new TBField;
	Bool_t magfield = bfield->SetMagField();		 
  
  //Read ROOT MC output file
  
  bool flagMC = true;
  TFile*filein=new TFile("/home/sarah/AESOPLITE/ALanalysis-master/ALanalysis-master/RawEventMC.root","READ");
  TTree *tree;
  if(flagMC)tree = (TTree*)filein->Get("MC");
  
   //Define variables to read event
   ALEvent *e = new ALEvent;      			
 //Set address to access event data
  tree->SetBranchAddress("event",&e); 
 // Get number of event in Tree
 int nentries=tree->GetEntries();

 for (int i=0;i<6;i++) {
   tree->GetEntry(i);
   //loop over hits in event
      int nnhits = (int)e->get_Nhits();
	 // TObjArray rawhits; 						//store original MC coordinates 
        //FOR NOW choose only events with 7 detector hits, process hits
		    if (nnhits == 7) {
		  Double_t E0 = e->get_EkMC() *1000 ;
		  cout << "Injection energy is " << E0 << " MeV " << endl;
          for(int j=0;j<nnhits;j++) {
          Float_t X=((e->get_hits().at(j))->get_xin()+(e->get_hits().at(j))->get_xout())/2;//Mean of xin and xout
          Float_t Z=((e->get_hits().at(j))->get_zin()+(e->get_hits().at(j))->get_zout())/2;
          Float_t Y=((e->get_hits().at(j))->get_yin()+(e->get_hits().at(j))->get_youtn())/2;
		  cout << "Hit " << j << " x = "<< Y << "   y= "<< Z << "   z= "<< X << endl;
          TVector3 xx;                      						//measurement vector
          xx.SetXYZ(X,Y,Z);	
		  ALMeasLayer &ms = *static_cast<ALMeasLayer *>(cradle.At(j));
          Bool_t bending = ms.IsBending();
	      ms.ProcessHit(xx, kalhits, bending);
		 // ms.ProcessRawHit(xx, rawhits, bending);
        }
       }
      
   //loop over the entries in kalhit buffer
	 
    if ((kalhits.GetEntries() == 7)) {
      
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
	  Double_t bfield = h1.GetBfield();
//	  cout << "The magnetic field magnitude is B = " << bfield << " T " << endl;
//	  cout << "Raw hit vector x1 =" << x1.X() << "  y1=" << x1.Y() << " z1=" << x1.Z() << endl;
//	  cout << "Raw hit vector x2 =" << x2.X() << "  y2=" << x2.Y() << " z2=" << x2.Z() << endl;
//	  cout << "Raw hit vector x3 =" << x3.X() << "  y3=" << x3.Y() << " z3=" << x3.Z() << endl;
      THelicalTrack helstart(x1, x2, x3, h1.GetBfield(), kDir); // initial helix 
      
      // ---------------------------
      //  Set dummy state to sited
      // ---------------------------
      
      static TKalMatrix svd(kSdim,1);
      svd(0,0) = 0.;
      svd(1,0) = 3.1654;
      svd(2,0) = -89.2282;
      svd(3,0) = 0.;
      svd(4,0) = 0.0960282;
      if (kSdim == 6) svd(5,0) = 0.;
	  cout << "rho=" << helstart.GetDrho() << " Phi0=" << helstart.GetPhi0() << "  , kappa=" << helstart.GetKappa() << "  , dz=" << helstart.GetDz() << "  TanL=" << helstart.GetTanLambda() << endl;


     static TKalMatrix C(kSdim,kSdim);
      for (Int_t i=0; i<kSdim; i++) {
         C(i,i) = 0.005;   // dummy error matrix
      }

      sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kPredicted));
      sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kFiltered));

      // ---------------------------
      //  Add sited to the kaltrack
      // ---------------------------

      TKalTrack kaltrack;    // a track is a kal system
      kaltrack.SetOwner();   // kaltrack owns sites
      kaltrack.Add(&sited);  // add the dummy site to the track

      
	TIter next(&kalhits, kDir);
    ALHit *hitp =  0;

       while ((hitp = dynamic_cast<ALHit *>(next())))  {
		   cout << "-----------------------Next Hit--------------------------- " << endl;
		   const ALMeasLayer &ml = dynamic_cast<const ALMeasLayer &>(hitp->GetMeasLayer());	
           TVector3 xv = ml.HitToXv(*hitp);
		   TVector3 xraw = hitp->GetRawXv();
		  TKalTrackSite  &site = *new TKalTrackSite(*hitp); // create a site for this hit
		   		TVector3 bfield = TBField::GetGlobalBfield(xv);
cout		 << "x = "  << xv.X()
			 << " y ="  << xv.Y()
			 << " z ="  << xv.Z()
			 << ", B = (" << bfield.X() << "," << bfield.Y() << "," << bfield.Z() << ")" 
			 << endl;
         if (!kaltrack.AddAndFilter(site)) {  		         // add and filter this site
        	site.DebugPrint();
         	kaltrack.GetState(TVKalSite::kFiltered).DebugPrint();
			cerr << " site discarded!" << endl;
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
	  cout << "\t xraw =("<< xraw.X()<<",  "<<xraw.Y()<<", "<<xraw.Z()<<"): \n";
	  cout << "\t xv   =("<< xv.X()<<",  "<<xv.Y()<<", "<<xv.Z()<<"): \n";
	  cout << "\t x_exp=("<< x_exp.X()<<",  "<<x_exp.Y()<<", "<<x_exp.Z()<<"): \n";
	  cout << "\t x_fil=("<< x_fil.X()<<",  "<<x_fil.Y()<<", "<<x_fil.Z()<<"): \n";
     cout << "\t momentum = " << mom * 1000 << " MeV \n";

			}   
       } 
      
	kaltrack.SmoothBackTo(6);                          // smooth back.
	
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
  	  Double_t p  = 1000 * pt * sqrt(1+tanl*tanl);   //p = pt / sinTheta 	
 	 cout << "rho=" << rho << endl;
	 cout << "phi0=" << phi0 << endl;
	 cout << "tanl=" << tanl << endl;
	 cout << "cpa" << cpa << endl;	
	 cout << "p is " << p << " MeV " << endl; 
    }
 //   kalhits.Clear();                         //clear buffer at the end of event
   // rawhits.Clear();      
     }   //end of loop over events
  
  filein->Close();
	
 }

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	