
#include "MainRecoEventMC.h"

using namespace std;

int main(int argc, char*argv[]) 
{

 if(argc!=8)
  {
   cout << "Wrong number of parameters!!!!" << endl;
   cout << "The program needs 5 input parameters:"<< endl;
   cout << "First is Fluka type of simulated particles" <<endl;
   cout << "Second is energy in MeV or GeV" <<endl;
   cout << "Third is random seed configuration" << endl;
   cout << "Fourth is first cycle to reconstruct (Starts at 1)" <<endl;
   cout << "Fifth is the number of cycle to reconstruct" <<endl;
   cout << "Six is tag of reconstruction" <<endl;
   cout << "Seven is MC source tag " << endl;
   return -1;
  }
 //Fluka type of particle
 int type=(int) atoi(argv[1]);           //3 for electrons
 int Ene=(int) atoi(argv[2]);           //Energy 
 int seed=(int) atoi(argv[3]);           //random seed configuration 
 int Ncycles=(int) atoi(argv[4]);      //first cycle to reconstruct
 int Ncycles2=(int) atoi(argv[5]);    //last cycle
 int DataType= 0;                    //datatype, 0=MC, 1 = data
 string RecoInd=argv[6];            //string Reco Index: allows to distinct between types of reconstruction
 string source=argv[7];
 int InitType=1;		   //0=MC, 1=PR, 2=3pt helix
 bool secondIter=true;
 
 //Load region numbers used in the MC geometry 
 int*TckReg=new int[7];
 int*TrigReg=new int[4];
 int*GReg=new int[1];
 int*ShellReg=new int[2];
 float*TckZPos=new float[7];
 float*TrigThresh=new float[4];
 float*GuardThresh=new float[1];
 for(int i=0;i<7;i++)TckReg[i]=0;
 for(int i=0;i<4;i++)TrigReg[i]=0;
 for(int i=0;i<2;i++)ShellReg[i]=0;
 for(int i=0;i<1;i++)GReg[i]=0;
 for(int i=0;i<7;i++)TckZPos[i]=0;
 for(int i=0;i<4;i++)TrigThresh[i]=0;
 for(int i=0;i<1;i++)GuardThresh[i]=0;
	
 string MCparamfile="../src/ALSim/MCparameters.dat"; 
 
 LoadMCparameters(MCparamfile,TckReg,TrigReg,GReg,TckZPos,TrigThresh,GuardThresh,ShellReg);	
 float ladderOffsetLeft[7] = { 8.889, 8.866, 8.904, 8.871, 8.886, 8.873, 8.884 };
 float ladderOffsetRight[7] = { 98.472, 98.447, 98.478, 98.465, 98.460, 98.466 };

 //Set Magnetic field map
 TBField *bfield = new TBField();
 bfield->SetUseUniformBfield(false);

 bool FlagMagF=false;
 FlagMagF=bfield->SetMagField();
 if(FlagMagF)cout << "Field Map loaded" << endl;
 else 
  {
   cout << "There is an issue when loading the Field Map" << endl;
   return 1;
  }  

 //Set Magnetic field map for Runge-Kutta method
string fN = "/home/smechbal/ANALYSIS_SOFT/src/RKFitter/fieldmap5mm.bin";
FieldMap *fM = new FieldMap(fN, "binary", 81);
 //filename structure

 string startfile="aesopliteNonUniB_V4";
 string endfile="_fort.99";
 
 //Input files 


 string Inppath="/home/smechbal/MCproduction/AESOPLITE/rootfiles/NonUniB/V4";
 //Output files
 string Outpath="/home/smechbal/MCproduction/AESOPLITE/rootfiles/NonUniB/V4";

 
 for(int j=Ncycles;j<Ncycles+Ncycles2;j++)//Number of cycles
      {
   TFile *file;
   TFile *fileout;		  
  // for electrons, energies in MeV
       //input file
       if(type==3 || type==4) {
       file=new TFile(Form("%s/%d/%s/RawEvent_%s_%d_%dMeV%d%03d%s.root",Inppath.c_str(),type,source.c_str(),startfile.c_str(),type,Ene,seed,j,endfile.c_str()),"READ");
       cout << "Input file is open" <<endl;
       //output file
      fileout=new TFile(Form("%s/%d/%s/RecoEvent_%s_%d_%dMeV%d%03d%s_%s.root",Outpath.c_str(),type,source.c_str(),startfile.c_str(),type,Ene,seed,j,endfile.c_str(),RecoInd.c_str()),"RECREATE");
  }
		  //for muons, energies in GeV
       else {
		  
	     cout << Form("%s/%d/RawEvent_%s_%d_%dGeV%03d%s.root",Inppath.c_str(),type,startfile.c_str(),type,Ene,j,endfile.c_str()) <<endl;
       //input file
       file=new TFile(Form("%s/%d/RawEvent_%s_%d_%dGeV%03d%s.root",Inppath.c_str(),type,startfile.c_str(),type,Ene,j,endfile.c_str()),"READ");
       cout << "Input file is open" <<endl;
       //output file
       fileout=new TFile(Form("%s/%d/RecoEvent_%s_%d_%dGeV%03d%s_%s.root",Outpath.c_str(),type,startfile.c_str(),type,Ene,j,endfile.c_str(),RecoInd.c_str()),"RECREATE");
 }
       //Get Tree from the input file
       TTree *tree = (TTree*)file->Get("MC");

       //Define variables to read event
       ALEvent *e = new ALEvent();
       //Set address to access event data
       tree->SetBranchAddress("event",&e);  

       // Create a TTree
       TTree *REtree = new TTree("MC"," Reco event MC");
       //Define variables to make the Reco event
       ALEvent *re=new ALEvent();
       // Create a branch with event
       REtree->Branch("Revent",&re); 
       // Get number of events in Tree
       int nentries=tree->GetEntries();
       cout << "Number  of events: " << nentries << endl;
       //Loop over the events
       // This is the main loop
       for (int k=0;k<nentries;k++)
//     	for(int k=0;k<1000;k++) 
    	{
	  
          tree->GetEntry(k); //Load the entry k in the variable e  
          //Copy the raw event into a reco event with same structure 
          re=new ALEvent();
          re->Copy(e);  
          int nnhits = re->get_Nhits();	  
          //Number of layers with hit(s)
          int NL= re->get_NLayers();
          //Number of layers with hit(s) in bending/non-bending plane
          int NLB = re->get_Layer(1) + re->get_Layer(2) + re->get_Layer(3) + re->get_Layer(5);
          int NLNB =  re->get_Layer(0) + re->get_Layer(4) + re->get_Layer(6);	

//	  if (nhit!=7) continue;
//	  cout << "Event " << k << ", nhits = " << nhit << endl;		 
          /////////////////////
          //Pattern Recognition
          /////////////////////
			  		   
          ALPatternRecognition* TestPattern = new ALPatternRecognition();
	  int PR = TestPattern->FindPattern(re,DataType,TckZPos,ladderOffsetLeft,ladderOffsetLeft,TrigThresh);	 	  
	  if(PR==0) {
	 	REtree->Fill();
		delete re;
		delete TestPattern;
		continue;
		}
	
	 //////////////////////////////
	 //Kalman Filter reconstruction
	 //////////////////////////////

	 else {
	  ALKalman* KalmanReco = new ALKalman(re);
	  int KF = KalmanReco->DoKF(re, DataType,InitType, secondIter);
          if(KF==0) {
          REtree->Fill();
          delete re;
	  delete KalmanReco;
          continue;
            }
	  delete KalmanReco;
	   }	
/*
	 //////////////////////////////
	 //Runge-Kutta reconstruction//
	 //////////////////////////////

 else 
	 	{
           if(NL ==7 && nnhits > 6) { 					//number of hit condition	                       		 
	   double event = re->get_eventnumber();
	   double*zLayers=new double[7];
	   double e0MC=0;
	   double p0MC=0;
	   double x0MC=0;
	   double y0MC=0;
	   double z0MC=0;
	   double cx0MC=0;
	   double cy0MC=0;
	   double cz0MC=0;
	   double kappa_MC=0;
	   double p0PR=0;
	   double x0PR=0;
	   double y0PR=0;
	   double z0PR=0;
	   double cx0PR=0;
	   double cy0PR=0;
	   double cz0PR=0;
	   double kappa_PR=0;
	   int typeMC =re->get_typeMC();
	   int QMC=0;
	   if(typeMC==3 || typeMC==11) QMC = -1;
	   if(typeMC==4 || typeMC==10 || typeMC==1) QMC = 1;
	   if(typeMC==-6) QMC = 2;
	   float mass=0.000511;//electon mass in GeV
	   if(typeMC==11 || typeMC==10) mass=0.10566;	//muon mass in MeV
	   if(typeMC==1) mass=0.93827;					//proton mass in GeV
	   if(typeMC==-6) mass=3.72739;					//alpha mass in GeV 


	   for(int i=0;i<7;i++)zLayers[i]=0;
		for(int j=0;j<nnhits;j++) 
		{ 	
			double zz=10*(re->get_hits().at(j))->get_z();			//in mm
			int L = int(re->get_hits().at(j)->get_mregMC())%11;
	                bool flagPR = re->get_hits().at(j)->get_flagPR();
			if(!flagPR) continue;
			zLayers[L] = zz;
			if(L==0)       
			{
				
				//MC variables on L0
				e0MC = (re->get_hits().at(j)->get_eMC());               //energy at L0 in GeV
				p0MC = TMath::Sqrt((e0MC*e0MC) - (mass*mass));
                x0MC=(re->get_hits().at(j)->get_x())*10;  // Randomly move the starting point and direction
                y0MC=(re->get_hits().at(j)->get_y())*10;  // Randomly move the starting point and direction           
                z0MC=(re->get_hits().at(j)->get_z())*10;
				cx0MC= re->get_hits().at(j)->get_cx();
				cy0MC= re->get_hits().at(j)->get_cy();
				cz0MC= re->get_hits().at(j)->get_cz();
				TVector3 p0(e0MC*cx0MC, e0MC*cy0MC, e0MC*cz0MC);                             		
				kappa_MC = QMC/p0MC;
				
				//PR variables at L0
			    p0PR = re->get_p0PR();
				double deflecPR = re->get_deflecPR();
				double QPR = TMath::Sign(1,deflecPR);
                x0PR=(re->get_hits().at(j)->get_xPR())*10;  // Randomly move the starting point and direction
                y0PR=(re->get_hits().at(j)->get_yPR())*10;  // Randomly move the starting point and direction           
                z0PR=(re->get_hits().at(j)->get_zPR())*10;
                cx0PR= re->get_hits().at(j)->get_cxPR();
                cy0PR= re->get_hits().at(j)->get_cyPR();
                cz0PR= re->get_hits().at(j)->get_czPR();
                kappa_PR = QPR/p0PR;						// TVector3 p0(e0MC*cx0MC, e0MC*cy0MC, e0MC*cz0MC)
				
				//cout << "Event " << event <<endl;
		        //cout << "x0MC = " << x0MC << "  y0MC = " << y0MC << "  cx0MC = " << cx0MC << "  cy0MC = " << cy0MC << " kappaMC = " << kappa_MC << endl;
                //cout << "x0PR = " << x0PR << "  y0PR = " << y0PR << "  cx0PR = " << cx0PR << "  cy0PR = " << cy0PR << " kappaPR = " << kappa_PR << endl;
	 }
		}			//end loop over hits
			 
	 // Create a data list to pass to the fitting routine
	char oLayer[7] = { 'n','b','b','b','n','b','n' };    // Orientation, nonbending or bending
	TkrData *Td = new TkrData();
	for (int lyr = 0; lyr < 7; lyr++) {Td->addLyr(oLayer[lyr], TckZPos[lyr], ladderOffsetLeft[lyr], ladderOffsetRight[lyr]);}
	for(int j=0;j<nnhits;j++)
		{ 	
		bool flagPR = re->get_hits().at(j)->get_flagPR();
		bool fGhost = re->get_hits().at(j)->get_fGhost();
		if(!flagPR) continue;
		if(fGhost) continue;
		int k = re->get_hits().at(j)->get_k();
		int L = int(re->get_hits().at(j)->get_mregMC())%11;
		int nstrip = re->get_hits().at(j)->get_nstrips();
		//coordinates fom MC or data
		double xx=(re->get_hits().at(j)->get_x())*10;			//in mm
		double yy=(re->get_hits().at(j)->get_y())*10;           //in mm
		double zz=(re->get_hits().at(j)->get_z())*10;			//in mm
		//coordinate from PR fit
		double xPR=(re->get_hits().at(j)->get_xPR())*10;		//in mm
		double yPR=(re->get_hits().at(j)->get_yPR())*10;        //in mm		
		double zPR=(re->get_hits().at(j)->get_zPR())*10;		//in mm	
		//Fill in hits informations onto TkrData class
		if(L==0||L==4||L==6) Td->addHit(L,xx);
		else Td->addHit(L,yy);	
		} // end loop on hits

	  //Do RK fit
	  	bool verbose = false;
		bool MCS = false;
	        double stepSize = 5;
		int alg = 0;
		RKfitter *rkf = new RKfitter(verbose, 0., fM, Td, MCS, stepSize, alg);
		vector<int> hits = {0,0,0,0,0,0,0};
   		double guess[5] = {x0PR,y0PR,cx0PR,cy0PR,kappa_PR};
		rkf->fitIt(false, guess, hits);
		 
	//Get fit parameters and fill event TTree	 
		 
		rkf->print("test fit");		
		double a[5];
		rkf->tkrParams(a);
	    double e[5];
		rkf->errors(e);
	  	double p0reco = abs(1./ a[4]);		//fitted momentum in GeV
		double chi2 = rkf->chiSqr();			//chi2 of fit
	    int typereco = TMath::Sign(1,a[4]);		//positive or negative particle?
	    double cxL0 = a[1];						//directional cosine at beginning of track
	  	double cyL0 = a[2];
		double err_cpa = e[4];
		double cpa = typereco/p0reco;
		for (int lyr = 0; lyr < Td->nLayers; lyr++)
		{
			double r[3];
			rkf->getIntercept(lyr, r); 
		    double xreco = r[0];
			double yreco = r[1];
			re->get_hits().at(lyr)->set_xreco(xreco);
			re->get_hits().at(lyr)->set_yreco(yreco);
			re->get_hits().at(lyr)->set_fUsed(true);
			if(lyr==0) 
					{
					re->get_hits().at(lyr)->set_cxreco(cxL0);
					re->get_hits().at(lyr)->set_cyreco(cyL0);
					} 
			}
            cout << " event " << event << "  p0reco = " << p0reco << "GeV" << "  p0MC = " << p0MC << " GeV " << endl;

 
	 //Fill in reconstruction variables in TTree 
	    re->set_p0reco(p0reco*1000);		//save p0reco in MeV
	    re->set_chi2(chi2);
	    re->set_typereco(typereco);
	    re->set_cpa(cpa);
	    re->set_cpaerr2(err_cpa);
	   delete rkf;
 
 		 } //end nnhit condition
  } //if reconstructed by PR		


*/
   
          /////////////////////    
          //Fill the output file with the reconstructed event
          /////////////////////
           REtree->Fill();
          //Free memory    
           delete re;
           delete TestPattern;

         }//k End loop on events
       delete e;
       //Write tree in output file
       fileout->cd();
       REtree->Write();
       //Close files
       fileout->Close();
       file->Close();
      // REtree->Delete();
      }//j

 //delete pointers
  delete[] TckReg;
 delete[] TrigReg;
 delete[] GReg;
 delete[] ShellReg;
 delete[] TckZPos;
 delete[] GuardThresh;
 delete[] TrigThresh;
 return 0;
}