
#include "MainRecoEventMC.h"

using namespace std;

int main(int argc, char*argv[]) 
{

 if(argc!=6)
  {
   cout << "Wrong number of parameters!!!!" << endl;
   cout << "The program needs 5 input parameters:"<< endl;
   cout << "First is Fluka type of simulated particles" <<endl;
   cout << "Second is energy in MeV or GeV" <<endl;
   cout << "Third is first cycle to reconstruct (Starts at 1)" <<endl;
   cout << "Fourth is the number of cycle to reconstruct" <<endl;
   cout << "Fifth is tag of reconstruction" <<endl;
   return -1;
  }
 //Fluka type of particle
 int type=(int) atoi(argv[1]);           //3 for electrons
 int Ene=(int) atoi(argv[2]);           //Energy 
 int Ncycles=(int) atoi(argv[3]);      //first cycle to reconstruct
 int Ncycles2=(int) atoi(argv[4]);    //last cycle
 int DataType= 0;                    //datatype, 0=MC, 1 = data
 string RecoInd=argv[5];            //string Reco Index: allows to distinct between types of reconstruction
 int InitType=1;		   //0=MC, 1=PR, 2=3pt helix
 bool secondIter=false;
 
 //Load region numbers used in the MC geometry 
 int*TckReg=new int[7];
 int*TrigReg=new int[4];
 int*GReg=new int[1];
 float*TckZPos=new float[7];
 float*TrigThresh=new float[4];
 float*GuardThresh=new float[1];
 for(int i=0;i<7;i++)TckReg[i]=0;
 for(int i=0;i<4;i++)TrigReg[i]=0;
 for(int i=0;i<1;i++)GReg[i]=0;
 for(int i=0;i<7;i++)TckZPos[i]=0;
 for(int i=0;i<4;i++)TrigThresh[i]=0;
 for(int i=0;i<1;i++)GuardThresh[i]=0;
 string MCparamfile="../src/ALSim/MCparameters.dat"; 
 LoadMCparameters(MCparamfile,TckReg,TrigReg,GReg,TckZPos,TrigThresh,GuardThresh);
 	

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

 //filename structure

 string startfile="aesopliteNonUniB_V3";
 string endfile="_fort.99";
 
 //Input files 


 string Inppath="/home/smechbal/MCproduction/AESOPLITE/rootfiles/NonUniB/V3/Disc";
//string Inppath="/home/smechbal/MCproduction/AESOPLITE/rootfiles/NonUniB/V2_NotDisc";
 //Output files
 string Outpath="/home/smechbal/MCproduction/AESOPLITE/rootfiles/NonUniB/V3/Disc";
//string Outpath="/home/smechbal/MCproduction/AESOPLITE/rootfiles/NonUniB/V2_NotDisc";

 
 for(int j=Ncycles;j<Ncycles+Ncycles2;j++)//Number of cycles
      {
   TFile *file;
   TFile *fileout;		  
  // for electrons, energies in MeV
       //input file
       if(type==3 || type==4) {
       file=new TFile(Form("%s/%d/RawEvent_%s_%d_%dMeV%03d%s.root",Inppath.c_str(),type,startfile.c_str(),type,Ene,j,endfile.c_str()),"READ");
       cout << "Input file is open" <<endl;
       //output file
      fileout=new TFile(Form("%s/%d/RecoEvent_%s_%d_%dMeV%03d%s_%s.root",Outpath.c_str(),type,startfile.c_str(),type,Ene,j,endfile.c_str(),RecoInd.c_str()),"RECREATE");
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
     	//for(int k=0;k<1000;k++) 
    	{
	  
          tree->GetEntry(k); //Load the entry k in the variable e  
          //Copy the raw event into a reco event with same structure 
          re=new ALEvent();
          re->Copy(e);  
	  
	 
			 
          /////////////////////
          //Pattern Recognition
          /////////////////////
			  		   
          ALPatternRecognition* TestPattern = new ALPatternRecognition();
          int PR = TestPattern->FindPattern(de,DataType,zL,OffsetLL,OffsetRL,TrigThresh);	 
	  if(PR==0) {
	 	REtree->Fill();
		delete re;
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
          continue;
            }

	   }	
		 
          /////////////////////    
          //Fill the output file with the reconstructed event
          /////////////////////
           REtree->Fill();
          //Free memory    
           delete re;
         }//k End loop on events
       delete e;
       //Write tree in output file
       fileout->cd();
       REtree->Write();
       //Close files
       fileout->Close();
       file->Close();
      }//j


 return 0;
}
