#include "MainPatternRecoEventMC.h"
using namespace std;

int main(int argc, char*argv[]) 
{

 if(argc!=6)
  {
   cout << "Wrong number of parameters!!!!" << endl;
   cout << "The program needs 5 input parameters:"<< endl;
   cout << "First is Fluka type of simulated particles" <<endl;
   cout << "Second is energy in MeV" <<endl;
   cout << "Third is first cycle to reconstruct (Starts at 1)" <<endl;
   cout << "Fourth is the number of cycle to reconstruct" <<endl;
   cout << "Fifth is tag of reconstruction" <<endl;
   return -1;
  }
 //Fluka type of particle
 int type=(int) atoi(argv[1]); //3 for electrons
 int Ene=(int) atoi(argv[2]);   //Energy in MeV
 int Ncycles=(int) atoi(argv[3]);   //first cycle to reconstruct
 int Ncycles2=(int) atoi(argv[4]);   //last cycle
  //string Reco Index: allows to distinct between types of reconstruction
 string RecoInd=argv[5];

 
 //Load region numbers used in the MC geometry 
 int*TckReg=new int[7];
 int*TrigReg=new int[4];
 int*GReg=new int[1];for(int i=0;i<1;i++)GReg[i]=0;
 for(int i=0;i<7;i++)TckReg[i]=0;
 for(int i=0;i<4;i++)TrigReg[i]=0;
 for(int i=0;i<1;i++)GReg[i]=0;
 string MCparamfile="../src/ALSim/MCparameters.dat"; 
 
 LoadMCparameters(MCparamfile,TckReg,TrigReg,GReg); 
  

 //filename structure
 string startfile="aesopliteUniB_V1";
 string endfile="_fort.99";
 
 //Input files 
 string Inppath="/home/smechbal/MCproduction/AESOPLITE/rootfiles/UniB/V1";

 //Output files
 string Outpath="/home/smechbal/MCproduction/AESOPLITE/rootfiles/UniB/V1";

   
 for(int j=Ncycles;j<Ncycles+Ncycles2;j++)//Number of cycles
      {
       cout << Form("%s/%d/RawEvent_%s_%d_%dMeV%03d%s.root",Inppath.c_str(),type,startfile.c_str(),type,Ene,j,endfile.c_str()) <<endl;

       //input file
       TFile*file=new TFile(Form("%s/%d/RawEvent_%s_%d_%dMeV%03d%s.root",Inppath.c_str(),type,startfile.c_str(),type,Ene,j,endfile.c_str()),"READ");
       cout << "Input file is open" <<endl;
       //output file
        TFile*fileout=new TFile(Form("%s/%d/PatternReco_%s_%d_%dMeV%03d%s.root",Outpath.c_str(),type,startfile.c_str(),type,Ene,j,endfile.c_str()),"RECREATE");

       //Get Tree from the input file
       TTree *tree = (TTree*)file->Get("MC");
       //Define variables to read event
       ALEvent *e = new ALEvent();
       //Set address to access event data
       tree->SetBranchAddress("event",&e);  
       // Create a TTree
       TTree *REtree = new TTree("MC"," Pattern Reco event MC");
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
       // for(int k=0;k<100; k++)
         {
		  cout << "----------------------EVENT " << k << " -----------------------------------" << endl;
          tree->GetEntry(k); //Load the entry k in the variable e  
          //Copy the raw event into a reco event with same structure 
          re=new ALEvent();
          re->Copy(e);  
			 
          /////////////////////
          //Pattern Recognition
          /////////////////////
			    
          //This is where selection is made on the event that can be reconstructed 
          
          //only select for pattern recognition event with at least one hit per layer inner trigger==127 for now
          if (re->get_Ti()==127)
           {
           ALPatternRecognition* TestPattern = new ALPatternRecognition();
		  TestPattern->FindPattern(re);

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
