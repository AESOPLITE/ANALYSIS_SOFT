#include "MainRecoEvent.h"
using namespace std;

int main() 
{
 //This function apply the reconstruction to only one file 
 //
 //These 3 variables will be in argument on the function
 //The executable can easily be call in a shell script to run over multiple files 
 int typeT=3; 
 int Ene=30;
 int cycle =1;
  
 //Load region numbers used in the MC geometry 
 int*TckReg=new int[7];
 int*TrigReg=new int[4];
 int*GReg=new int[1];for(int i=0;i<1;i++)GReg[i]=0;
 for(int i=0;i<7;i++)TckReg[i]=0;
 for(int i=0;i<4;i++)TrigReg[i]=0;
 for(int i=0;i<1;i++)GReg[i]=0;
 string MCparamfile="/src/ALSim/MCparameters.dat"; 
 
 LoadMCparameters(MCparamfile,TckReg,TrigReg,GReg); 
  
 //Set Magnetic field map
 TFile*file_map=new TFile("/prod/fieldmap.root","READ");
 TBField *bfield = new TBField();
 bool FlagMagF=false;
 FlagMagF=bfield->SetMagField();
 if(FlagMagF)cout << "Field Map INTERFACED WOOO" << endl;
 else 
  {
   cout << "There is an issue when loading the Field Map" << endl;
   return 1;
  }
/* 	
 //filename structure
 string startfile="aesoplite";
 string endfile="_fort.99"; 
  
 //Input file 
 string Inppath="/home/psm/Documents/UDEL/AESOPLITE/Github/ALanalysis/rootfiles/R20_Alpha11_54";
 TFile*file   =new TFile(Form("%s/RawEvent_%s_%d_%dMeV%03d%s.root",Inppath.c_str(),startfile.c_str(),typeT,Ene,cycle,endfile.c_str()),"READ");
 //Output file  
 string Outpath="/home/psm/Documents/UDEL/AESOPLITE/Github/ALanalysis/rootfiles/R20_Alpha11_54";
 TFile*fileout=new TFile(Form("%s/RecoEvent_%s_%d_%dMeV%03d%s.root",Outpath.c_str(),startfile.c_str(),typeT,Ene,cycle,endfile.c_str()),"RECREATE");
  */
 
TFile*file=new TFile("/home/sarah/AESOPLITE/ALanalysis-master/ALanalysis-master/RawEvents/RawEventMC_100MeV_uniform.root","READ");

 	
 //Get Tree from the input file
 TTree *tree = (TTree*)file->Get("MC");
 //Define variables to read event
 ALEvent *e = new ALEvent;
 //Set address to access event data
 tree->SetBranchAddress("event",&e);  
  
  
  // Create a TTree
 TTree *REtree = new TTree("MC"," Reco event MC");
 //Define variables to make the Reco event
 ALEvent *re=new ALEvent;
 // Create a branch with event
 REtree->Branch("Revent",&re);  
  
 // Get number of event in Tree
 int nentries=tree->GetEntries();
 cout << "Number  of events: " << nentries << endl;
 //Loop over the events
 // This is the main loop
	for (int i=0;i<100;i++)
   {
    tree->GetEntry(i); //Load the entry i in the variable e  
    
    //Copy the raw event into a reco event with same structure 
    re=new ALEvent;
    re->Copy(e);  
    /////////////////////
    //Pattern Recognition
    /////////////////////
    //TestPC.FindPattern();

    /////////////////////
    //Track Reconstruction
    /////////////////////
    
	int nnhits = (int)re->get_Nhits();
	   
	   //only reconstruct events with 7 hits
 if (nnhits == 7) {
    ALKalman* TestKalman= new ALKalman();
    TestKalman->MakeRecoEvent(bfield,re);
 }
	double Ekreco=0;
    double X0reco=0;
    double Y0reco=0;
    double Z0reco=0;
    double CX0reco=0;
    double CY0reco=0;
    double CZ0reco=0;
    
    /////////////////////
    //Fill the reconstructed information
    /////////////////////
    re->set_Ekreco(Ekreco);
    re->set_X0reco(X0reco);
    re->set_Y0reco(Y0reco);
    re->set_Z0reco(Z0reco);
    re->set_CX0reco(CX0reco);
    re->set_CY0reco(CY0reco);
    re->set_CZ0reco(CZ0reco);   
    
    //Fill the output file with the reconstructed event
    /////////////////////
    //REtree->Fill();
    //Free memory    
    delete re;
   }//End loop on event

 delete e;
 
 
 //Write tree in output file
// fileout->cd();
// REtree->Write();
 //Close files
// fileout->Close();
 file->Close();
 

 return 0;
}
