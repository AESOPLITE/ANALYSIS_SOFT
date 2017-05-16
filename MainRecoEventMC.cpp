#include "MainRecoEventMC.h"
using namespace std;

int main() 
{
 //This function apply the reconstruction

 //Fluka type of particle
 int type=3; //3 for electrons
 //Number of energies
 int Nene=12;
 //Energies
 int Ene[12]={10,20,30,40,50,60,70,80,90,100,200,300};
 //Number of cycles per energy
 int* Ncycles=new int[Nene]; 
 for(int i=0;i<Nene;i++)Ncycles[i]=100;  ;
  
 
 //Load region numbers used in the MC geometry 
 int*TckReg=new int[7];
 int*TrigReg=new int[4];
 int*GReg=new int[1];for(int i=0;i<1;i++)GReg[i]=0;
 for(int i=0;i<7;i++)TckReg[i]=0;
 for(int i=0;i<4;i++)TrigReg[i]=0;
 for(int i=0;i<1;i++)GReg[i]=0;
 string MCparamfile="../src/ALSim/MCparameters.dat"; 
 
 LoadMCparameters(MCparamfile,TckReg,TrigReg,GReg); 
  
 //Set Magnetic field map
 //TFile*file_map=new TFile("fieldmap.root","READ");
 TBField *bfield = new TBField();
 bool FlagMagF=false;
 FlagMagF=bfield->SetMagField();
 if(FlagMagF)cout << "Field Map INTERFACED WOOO" << endl;
 else 
  {
   cout << "There is an issue when loading the Field Map" << endl;
   return 1;
  }

 //filename structure
 string startfile="aesopliteUniB_V1";
 string endfile="_fort.99";
 
 //Input files 
 string Inppath="/home/psmangeard/MCproduction/AESOPLITE/rootfiles/UniB/V1";

 //Output files
 string Outpath="/home/psmangeard/MCproduction/AESOPLITE/rootfiles/UniB/V1";
  
 //string Reco Index: allows to distinct between types of reconstruction
 string RecoInd="Test1";
 
 for(int i=0;i<Nene;i++)//Energies
// for(int i=0;i<1;i++)//Energies
   {
    for(int j=0;j<Ncycles[i];j++)//Number of cycles
//    for(int j=0;j<1;j++)//Number of cycles
      {
       cout << Form("%s/%d/RawEvent_%s_%d_%dMeV%03d%s.root",Inppath.c_str(),type,startfile.c_str(),type,Ene[i],j+1,endfile.c_str()) <<endl;

       //input file
       TFile*file=new TFile(Form("%s/%d/RawEvent_%s_%d_%dMeV%03d%s.root",Inppath.c_str(),type,startfile.c_str(),type,Ene[i],j+1,endfile.c_str()),"READ");
       cout << "Input file is open" <<endl;
       //output file
       TFile*fileout=new TFile(Form("%s/%d/RecoEvent_%s_%d_%dMeV%03d%s_%s.root",Outpath.c_str(),type,startfile.c_str(),type,Ene[i],j+1,endfile.c_str(),RecoInd.c_str()),"RECREATE");
      
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
       // Get number of events in Tree
       int nentries=tree->GetEntries();
       cout << "Number  of events: " << nentries << endl;
       //Loop over the events
       // This is the main loop
       for (int k=0;k<nentries;k++)
         {
          tree->GetEntry(k); //Load the entry k in the variable e  
          //Copy the raw event into a reco event with same structure 
          re=new ALEvent;
          re->Copy(e);  
          /////////////////////
          //Pattern Recognition
          /////////////////////

          /////////////////////
          //Track Reconstruction
          ///////////////////// 
          int nnhits = (int)re->get_Nhits();
   
          //This is where selection is made on the event that can be reconstructed 
          
          //only reconstruct events with 7 hits, one hit per layer inner trigger==127
          if (nnhits == 7 && re->get_Ti()==127)
           {
            ALKalman* TestKalman= new ALKalman();
            TestKalman->MakeRecoEvent(bfield,re);
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
   }//i

 return 0;
}