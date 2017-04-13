////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, December 8, 2016
////////////////////////////////////////////////////////////////////////////////////////// 


#include "ALEvent.h"
ClassImp(ALEvent)

int ReadRawEvent(bool);
/////////////////////////////
//Example of How to read the Raw events
//This example show how to make plots from raw events
/////////////////////////////
int ReadRawEvent(bool flagMC=true)
{
////////////////////
// Input parameters
////////////////////
//flagMC= true if MC(default)  
//flagMC=false if data 
////////////////////
  
  //Output file 
 TFile*filein=new TFile("RawEventMC.root","READ");
 
 //Get Tree from the input file
 TTree *tree;
 if(flagMC)tree = (TTree*)filein->Get("MC");
 else tree = (TTree*)filein->Get("Data"); //To be defined
  
 //Define variables to read event
 ALEvent *e = new ALEvent;

 //Set address to access event data
 tree->SetBranchAddress("event",&e); 

 // Get number of event in Tree
 int nentries=tree->GetEntries();
 cout << "Number  of events: " << nentries << endl;

 //Declare histograms
 //Example: Energy deposited in T1 in MeV
 TH1F*HistEneT1=new TH1F("HistEneT1","HistEneT1",100,0.,10);
 HistEneT1->GetXaxis()->SetTitle("Deposited energy in T1 (MeV)");
 HistEneT1->GetYaxis()->SetTitle("Events");
 HistEneT1->GetYaxis()->SetTitleOffset(2.);
 HistEneT1->SetLineColor(kBlack);
 //Loop over the events
 // This is the main loop
 for (int i=0;i<nentries;i++)
   {
    tree->GetEntry(i); //Load the entry i in the variable e 
    
    if(i%10000==0) cout << "Event: " << i <<endl;
    
    //Get the number of energy deposits T1 (size of the vector)
    int nT1=e->get_EneT1().size();
    //Loop over the number of energy deposits 
    for(int j=0;j<nT1;j++)
      {
       //Get the value 	
       double temp=e->get_EneT1().at(j);
       //Make your selection or calculatoin here
       temp*=1000.; // Here we change energy unit from GeV to MeV
       //Fill the histograms
       HistEneT1->Fill(temp);	
      }
    
   }//End loop on event

////////////////////   
//Display  
////////////////////
 
 TCanvas*can=new TCanvas("can","can",200,150,500,500);
 can->cd();
 gPad->SetLeftMargin(0.15);
 
 //Draw the histogram
 HistEneT1->Draw("hist");
////////////////////
//Close files
////////////////////
 //filein->Close();
 return 1;
}