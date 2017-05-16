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
#include <cmath>


using namespace std;
static const double kMelectron = 0.5109989461e-3;		//mass electron in GeV
static const double RestMass   = 0.51099;				//mass electron in MeV
static const double sigmaXZ = 0.0066;
static const double sigmaY =  0.0115;
static const Bool_t kDir = kIterBackward;				//direction for the filter (starts from last hit)
static const double CovMElement=5.0e-2;			//initial covariant matrix elements
//static const Bool_t kDir = kIterForward;				//first hit to last (does not work as well...)

//Class constructor



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
	Int_t ndf;
	Double_t chi2, cl, d0, phi0, cpa, dz, d0err2, phi0err2, cpaerr2, dzerr2, tanlerr2;
	Double_t tanl, p0, kE_reco, pullP, pullX, pullY, pullZ, X0_reco, Y0_reco, Z0_reco;
  // ===================================================================
  //  Prepare a Root tree output
  // ===================================================================	 
	 
	
	TFile hfile("KalRecon_100MeV_Smoothed.root","RECREATE","KalTest");
  	TTree *t = new TTree("t", "ALKalman");

 	
	
  // ===================================================================
  // Track parameters and errors
  // ===================================================================	
	t->Branch("ndf", &ndf,"ndf/I");
    t->Branch("chi2", &chi2,"chi2/D");
	t->Branch("cl", &cl,"cl/D");
	t->Branch("d0", &d0,"d0/D");
	t->Branch("phi0", &phi0,"phi0/D");
	t->Branch("cpa", &cpa,"ndf/D");
	t->Branch("dz", &dz,"dz/D");
	t->Branch("d0err2", &d0err2,"d0err2/D");
	t->Branch("phi0err2", &phi0err2,"phi0err2/D");
	t->Branch("cpaerr2", &cpaerr2,"cpaerr2/D");
	t->Branch("dzerr2", &dzerr2,"dzerr2/D");
	t->Branch("tanlerr2", &tanlerr2,"tanlerr2/D");
	t->Branch("tanl", &tanl,"tanl/D");

  // ===================================================================
  //  Pull distributions
  // ===================================================================	
	
	
    t->Branch("p0", &p0,"p0/D");
    t->Branch("kE_reco", &kE_reco,"kE_reco/D");
	t->Branch("X0_reco", &X0_reco,"X0_reco/D");
	t->Branch("Y0_reco", &Y0_reco,"Y0_reco/D");
	t->Branch("Z0_reco", &Z0_reco,"Z0_reco/D");
    t->Branch("pullP", &pullP,"pullP/D");
    t->Branch("pullX", &pullX,"pullX/D");
    t->Branch("pullY", &pullY,"pullY/D");
    t->Branch("pullZ", &pullZ,"pullZ/D");

	

	TNtupleD *hTrackMonitor = new TNtupleD("Reconstruction", "", "ndf:chi2:cl:d0:phi0:cpa:dz:tanl:d0err2:phi0err2:cpaerr2:dzerr2:tanlerr2:p0:pullP");
	TNtupleD *pullMonitor = new TNtupleD("Pulls", "", "pullX:pullY:pullZ:kE_reco:X0_reco:Y0_reco:Z0_reco");

  // ===================================================================
  //  Read R00T MC input file
  // ===================================================================

	Bool_t flagMC = true;
    TFile*filein=new TFile("/home/sarah/AESOPLITE/ALanalysis-master/ALanalysis-master/RawEvents/RawEventMC_100MeV_uniform.root","READ");
    TTree *tree;
  	if(flagMC)tree = (TTree*)filein->Get("MC");
  
   ALEvent *e = new ALEvent;      			
   tree->SetBranchAddress("event",&e); 
   int nentries=tree->GetEntries();
   for (int i=0;i<nentries;i++) 
     {//loop over all events
      tree->GetEntry(i);
      int nnhits = (int)e->get_Nhits();
      Double_t E0 = e->get_EkMC() *1000 ;
      Double_t pMC = 100;
      if (nnhits == 7)
       {       //until PatternRecognition ready
        Double_t E0 = e->get_EkMC() *1000 ;
        Double_t X0 = e->get_X0MC();			//injection point
        Double_t Y0 = e->get_Y0MC();			//injeciton point
        Double_t Z0 = e->get_Z0MC();			//injection point 

 // cout << "Event " << i <<	"	X0=" << Y0 << "  Y0=" << Z0 << "  Z0=" << X0 << endl;
        for(int j=0;j<nnhits;j++) 
          {
           Float_t X=((e->get_hits().at(j))->get_xin()+(e->get_hits().at(j))->get_xout())/2;
           Float_t Z=((e->get_hits().at(j))->get_zin()+(e->get_hits().at(j))->get_zout())/2;
           Float_t Y=((e->get_hits().at(j))->get_yin()+(e->get_hits().at(j))->get_yout())/2;
           TVector3 xx;                      						
           xx.SetXYZ(X,Y,Z);
           TVector3 xv;
           xv.SetXYZ(Y,Z,X);
           TVector3 bfield = TBField::GetGlobalBfield(xv);	
           ALMeasLayer &ms = *static_cast<ALMeasLayer *>(cradle.At((j*2)));					//point to kActive layer
           Bool_t bending = ms.IsBending();
           Bool_t active = ms.IsActive();
           ms.ProcessHit(xx, kalhits, bending);
          }//j
       }//if
      
 
    if ((kalhits.GetEntries() == 7))
     {	
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
      if(bApply2Iter)
       {
        C = C_start;				 
       }
      else
       {
        for (Int_t i=0; i<kSdim; i++)
          {
           C(i,i) = CovMElement;         // dummy error matrix
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

      while ((hitp = dynamic_cast<ALHit *>(next()))) 
        {
         const ALMeasLayer &ml = dynamic_cast<const ALMeasLayer &>(hitp->GetMeasLayer());	
         TVector3 xv = ml.HitToXv(*hitp);
         TVector3 xraw = hitp->GetRawXv();
         TKalTrackSite  &site = *new TKalTrackSite(*hitp);   // create a site for this hit
         TVector3 bfield = TBField::GetGlobalBfield(xv);
         if (!kaltrack.AddAndFilter(site))
          {// add and filter this site
           //cerr << " site discarded!" << endl;
           delete &site;
           } //if
         else
          {
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
           pullX = xraw.X()-x_fil.X();
           pullY = xraw.Y()-x_fil.Y();
           pullZ = xraw.Z()-x_fil.Z();
  
           t->Fill();	   
           pullMonitor->Fill(pullX/sigmaXZ, pullY/sigmaY, pullZ/sigmaXZ);
          }   //else
        }//while 

        // =======================================================================
        //  Find tangent to helix at first site and extrapolate injection point
        // =======================================================================   

        TVKalState *state_first =(TVKalState*) &(kaltrack.GetCurSite().GetCurState());
        THelicalTrack hel_first = (dynamic_cast<TKalTrackState *>(state_first))->GetHelix();
        TVector3 pivot = hel_first.CalcXAt(0.0);
        TMatrixD dxdphi = hel_first.CalcDxDphi(0.0);					// tangent vector at destination surface
        TVector3 vtan(dxdphi(0,0),dxdphi(1,0),dxdphi(2,0));				// convert matrix diagonal to vector
        //cout << "Tangent vector vtan=("<< vtan.X()<<",  "<<vtan.Y()<<", "<<vtan.Z()<<"): \n"; 

        //Parametric equation of tangent line
        Double_t yinjection = 35;
        Double_t t = (yinjection - pivot.Y())/(vtan.Y());
        X0_reco = pivot.X() + vtan.X() * t;
        Y0_reco = pivot.Y() + vtan.Y() * t;
        Z0_reco = pivot.Z() + vtan.Z() * t;
        //cout << "Reconstructed injection point	X0= " << X0_reco << "  Y0=" << yinjection << "  Z0=" << Z0_reco << endl;

        //Smooth back to first site	
        Int_t isite = 1;
        kaltrack.SmoothBackTo(isite);									 
        TVKalSite &cursite = static_cast<TVKalSite &>(*kaltrack[isite]);

        // ============================================================
        //  Monitor Fit Result
        // ============================================================ 

        Int_t    ndf  = kaltrack.GetNDF();									//degrees of freedom	
        Double_t chi2 = kaltrack.GetChi2();									//chi2
        Double_t cl   = TMath::Prob(chi2, ndf);								//confidence level
        Double_t d0  =  cursite.GetCurState()(0, 0 ); 
        Double_t phi0  = cursite.GetCurState()(1, 0 ); 
        Double_t cpa  = cursite.GetCurState()(2, 0 ); 
        Double_t dz   = cursite.GetCurState()(3, 0 ); 
        Double_t tanl  = cursite.GetCurState()(4, 0 ); 
        //  cout << "Event " << i << "  1/tanl = " << 1/(tanl) << endl;
        const TKalMatrix& covK = cursite.GetCurState().GetCovMat() ; 

        // errors^2 of track parameters
        Double_t d0err2  = covK( 0 , 0 )   ;
        Double_t phi0err2 = covK( 1 , 1 )   ;
        Double_t cpaerr2 = covK( 2 , 2 )   ;
        Double_t dzerr2  = covK( 3 , 3 )   ;
        Double_t tanlerr2 = covK( 4 , 4 )   ;

        Double_t pt = fabs(1.0/cpa);
        Double_t pz = pt * tanl;
        Double_t p0  = 1000 * pt * sqrt(1+tanl*tanl);
        cout << "Final reconstructed momentum p0=" << p0 << endl;
        Double_t kE_reco =  sqrt(p0*p0 + (RestMass)*(RestMass))  - RestMass;	//kinetic energy of incoming particle in MeV
        //Double_t pullP = (pMC - p0)/34.65;	//replace numerical value with var. sigma_p
        pullMonitor->Fill(kE_reco, X0_reco, Y0_reco, Z0_reco);
        hTrackMonitor->Fill(ndf, chi2, cl, d0, phi0, cpa, dz, tanl, d0err2, phi0err2, cpaerr2, dzerr2, tanlerr2 , p0, pullP);
        //t->Fill();
        counter++;
        kaltrack.Delete();
       }//while
      kalhits.Delete();                       //clear buffer at the end of event
     }//end of loop over events

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
         C(k,k) = CovMElement;   								// dummy error matrix
      }

      sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kPredicted));
      sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kFiltered));

      // ---------------------------
      //  Add sited to the kaltrack
      // ---------------------------

      TKalTrack kalfirst;    // a track is a kal system
      kalfirst.SetMass(kMelectron);
      kalfirst.SetOwner();   // kaltrack owns sites
      kalfirst.Add(&sited);  // add the dummy site to the track

      
      TIter next(&kalhits, kDir);
      ALHit *hitp =  0;

       while ((hitp = dynamic_cast<ALHit *>(next())))
         {
          const ALMeasLayer &ml = dynamic_cast<const ALMeasLayer &>(hitp->GetMeasLayer());	
          TVector3 xv = ml.HitToXv(*hitp);
          TVector3 xraw = hitp->GetRawXv();
          TKalTrackSite  &site = *new TKalTrackSite(*hitp); // create a site for this hit
          TVector3 bfield = TBField::GetGlobalBfield(xv);
          if (!kalfirst.AddAndFilter(site))
           { // add and filter this site
            delete &site;
           }   
         }
  
  // kalfirst.SmoothBackTo(1);                          // smooth back to first site

  // ============================================================
  //  Get the last site then return it back
  // ============================================================
    
   TKalTrackState *theLastState = dynamic_cast<TKalTrackState*> (&(kalfirst.GetCurSite().GetCurState()));
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
	
   TFile*file = new TFile("KalRecon_100MeV_Smoothed.root","READ");
   TNtupleD *ntuple = (TNtupleD*)file->Get("Reconstruction");
   TNtupleD *pulls = (TNtupleD*)file->Get("Pulls");
   TH1D *p0hist = new TH1D("p0hist", "Reconstructed momentum 100MeV electrons", 300, 0.0, 300);
   TH1D *confl = new TH1D("confl", "P(chi2,ndf) distribution", 100, 0.0, 1);
   TH1D *pullP = new TH1D("pullp", "Momentum pull", 100, -10, 10);
   TH1D *pull_X = new TH1D("pullx", "X pull", 100, -10, 10);
   TH1D *pull_Y = new TH1D("pully", "Y pull", 100, -10, 10);
   TH1D *pull_Z = new TH1D("pullz", "Z pull", 100, -10, 10);
   Double_t pMC = 100;
  
   gStyle->SetOptFit(1111);
   TCanvas* c1 = new TCanvas("c1", "reconstruction", 1);
   c1->Divide(3,2);
   c1->cd(1);
   p0hist->GetXaxis()->SetTitle("p0 (in Mev/c)");
   p0hist->GetXaxis()->SetRangeUser(0, 300);
   ntuple->Draw("p0>>p0hist");			
   Double_t sigma_p = p0hist->GetStdDev();

   c1->cd(2);
   
   pullP->GetXaxis()->SetTitle("(pMC - pFit)");
   pullP->GetXaxis()->SetRangeUser(-10, 10);
   ntuple->Draw("pullP>>pullp");
  // pullP->Sumw2();
 //  pullP->Divide(sigma);
   //pullP->Draw();
	
	c1->cd(3);
  //ntuple->Draw("phi0");
   pull_X->GetXaxis()->SetTitle("(xMC-xReco)/sigmaXZ");
   pull_X->GetXaxis()->SetRangeUser(-10,10);
   pulls->Draw("pullX>>pullx");
	
	
  c1->cd(4);
 // ntuple->Draw("tanl");
   pull_Y->GetXaxis()->SetTitle("(yMC - yReco)/sigmaY");
   pull_Y->GetXaxis()->SetRangeUser(-10, 10);
   pulls->Draw("pullY>>pully");
   c1->cd(5);
  //ntuple->Draw("cpa");
   pull_Z->GetXaxis()->SetTitle("(zMC - zReco)/sigmaXZ");
   pull_Z->GetXaxis()->SetRangeUser(-10, 10);
   pulls->Draw("pullZ>>pullz");


  c1->cd(6);
  ntuple->Draw("cl>>confl");
 
   
  
	
  c1->Update();
  c1->Print("Recontruction_100MeV.pdf");
	
}

//Function to reconstruct one event at the time 

void ALKalman::MakeRecoEvent(TBField *bfield, ALEvent *re)
{
   // ===================================================================
   //  Prepare a detector
   // ===================================================================
   TObjArray kalhits; 
   TKalDetCradle	cradle;				
   ALKalDetector detector;
   cradle.Install(detector); 	
   
   Bool_t bApply2Iter = true;									//if initialize with first fit
   int nnhits = (int)re->get_Nhits();
   
   for(int j=0;j<nnhits;j++) 
     {
      Float_t X=((re->get_hits().at(j))->get_xin()+(re->get_hits().at(j))->get_xout())/2;//Mean of xin and xout
      Float_t Z=((re->get_hits().at(j))->get_zin()+(re->get_hits().at(j))->get_zout())/2;
      Float_t Y=((re->get_hits().at(j))->get_yin()+(re->get_hits().at(j))->get_yout())/2;
      TVector3 xx;                      						
      xx.SetXYZ(X,Y,Z);	
      //cout << "Hitxx " << j << " x = "<< xx.Y() << "   y= "<< xx.Z()<< "   z= "<< xx.X() << endl;
      ALMeasLayer &ms = *static_cast<ALMeasLayer *>(cradle.At(j));
      Bool_t bending = ms.IsBending();
      ms.ProcessHit(xx, kalhits, bending);
     }//end j
      
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
   if(bApply2Iter)
    {
     C = C_start;
    }
   else
    {
     for (int ij=0; ij<kSdim; ij++)
       {
	C(ij,ij) = CovMElement;         // dummy error matrix
       }//end ij
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

   while ((hitp = dynamic_cast<ALHit *>(next()))) 
     {
      //cout << "-----------------------Next Hit--------------------------- " << endl;
      const ALMeasLayer &ml = dynamic_cast<const ALMeasLayer &>(hitp->GetMeasLayer());	
      TVector3 xv = ml.HitToXv(*hitp);
      TVector3 xraw = hitp->GetRawXv();
      TKalTrackSite  &site = *new TKalTrackSite(*hitp); // create a site for this hit
      TVector3 TVbfield = TBField::GetGlobalBfield(xv);
      //cout << "x = "  << xv.X() << " y ="  << xv.Y() << " z ="  << xv.Z();
      //cout << ", B = (" << TVbfield.X() << "," << TVbfield.Y() << "," << TVbfield.Z() << ")"  << endl;
      if (!kaltrack.AddAndFilter(site))
       { // add and filter this site
        site.DebugPrint();
        kaltrack.GetState(TVKalSite::kFiltered).DebugPrint();
        cerr << " site discarded!" << endl;
        delete &site;
       } 
      else 
       {
        TVKalState *state_fil = (TVKalState*) &(site.GetCurState());
        THelicalTrack hel_fil = (dynamic_cast<TKalTrackState *>(state_fil))->GetHelix();
        TVector3 x_fil=hel_fil.CalcXAt(0.0);
        TVKalState *state_exp = &site.GetState(TVKalSite::kPredicted);
        THelicalTrack hel_exp = (dynamic_cast<TKalTrackState *>(state_exp))->GetHelix();
        TVector3 x_exp=hel_exp.CalcXAt(0.0);
        Double_t mom = hel_fil.GetMomentum();
        //debug: comparing the expected and filtered points
        //cout << "\t xraw =("<< xraw.X()<<",  "<<xraw.Y()<<", "<<xraw.Z()<<"): \n";
        //cout << "\t xv   =("<< xv.X()<<",  "<<xv.Y()<<", "<<xv.Z()<<"): \n";
        //cout << "\t x_exp=("<< x_exp.X()<<",  "<<x_exp.Y()<<", "<<x_exp.Z()<<"): \n";
        //cout << "\t x_fil=("<< x_fil.X()<<",  "<<x_fil.Y()<<", "<<x_fil.Z()<<"): \n";
        //cout << "\t momentum = " << mom * 1000 << " MeV \n";
       }//end else   
      } //end while
	
	  
   // =======================================================================
   //  Find tangent to helix at first site and extrapolate injection point
   // =======================================================================   

   TVKalState *state_first =(TVKalState*) &(kaltrack.GetCurSite().GetCurState());
   THelicalTrack hel_first = (dynamic_cast<TKalTrackState *>(state_first))->GetHelix();
   TVector3 pivot = hel_first.CalcXAt(0.0);
   TMatrixD dxdphi = hel_first.CalcDxDphi(0.0);					// tangent vector at destination surface
   TVector3 vtan(dxdphi(0,0),dxdphi(1,0),dxdphi(2,0));				// convert matrix diagonal to vector
   //cout << "Tangent vector vtan=("<< vtan.X()<<",  "<<vtan.Y()<<", "<<vtan.Z()<<"): \n";

   //Parametric equation of tangent line
   Double_t yinjection = 35; 
   Double_t t = (yinjection - pivot.Y())/(vtan.Y());
   double X0reco = pivot.X() + vtan.X() * t;
   double Y0reco = pivot.Y() + vtan.Y() * t;
   double Z0reco = pivot.Z() + vtan.Z() * t;
   
    //Get directional cosines of tangent line
   Double_t theta = vtan.Theta();
   Double_t phi   = vtan.Phi();
   Double_t CX0reco = sin(theta) * cos(phi);
   Double_t CY0reco = sin(theta) * sin(phi);
   Double_t CZ0reco = cos(theta);
   kaltrack.SmoothBackTo(1);                          // smooth back to first site
 
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
   double cpa  = kaltrack.GetCurSite().GetCurState()(2, 0);
   re->set_cpa(cpa);
   double tanl  = kaltrack.GetCurSite().GetCurState()(4, 0);
   re->set_tanl(tanl);
   double phi0  = kaltrack.GetCurSite().GetCurState()(1, 0);
   re->set_phi0(phi0);
   double rho =  kaltrack.GetCurSite().GetCurState()(0, 0);
   double pt =0;
   if(cpa!=0)pt=fabs(1.0/cpa);
   double pz = pt * tanl;
   double p0  = 1000 * pt * sqrt(1+tanl*tanl);   //p = pt / sinTheta 	
   re->set_p0reco(p0);
   double Ekreco =  sqrt(p0*p0 + (RestMass)*(RestMass))  - RestMass;
   re->set_Ekreco(Ekreco);
   
   //cout << "Reconstructed injection point	X0= " << X0reco << "  Y0=" << yinjection << "  Z0=" << Z0reco << endl;
   //cout << "p = " << p0 << " MeV" << endl; 

 }
 
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
