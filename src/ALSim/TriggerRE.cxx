////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, January 10, 2017
////////////////////////////////////////////////////////////////////////////////////////// 


#include "TriggerRE.h"

////////////////////////////
//Analyse triggers : coincidence, efficiciency etc...
/////////////////////////////
int TriggerRE(bool flagMC=true)
{
////////////////////
// Input parameters
////////////////////
//flagMC= true if MC(default)  
//flagMC=false if data 
////////////////////

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
 int L=7;//number of tracker layers
 int E=8;//number of energies used here
 float Ene[8]={10.,20.,30.,50.,70.,100.,200.,300.};//Energies in MeV
 int NTConfig=6;//Number of configuration of triggers
 string TConfig[6]={"T1","T3","T4","T1&T3","T1&T4","T1&T3&T4"};
 //histograms of the position X(cm) of the hits
 TH1F****HitPosX=new TH1F***[L];

 for(int i=0;i<L;i++)
   {
    HitPosX[i]=new TH1F**[NTConfig];
    for(int j=0;j<NTConfig;j++)
      {
       HitPosX[i][j]=new TH1F*[E];

       for (int k=0;k<E;k++)
         {
          HitPosX[i][j][k]=new TH1F(Form("X, L%d, %s, E=%3d MeV",i+1,TConfig[j].c_str(),(int)Ene[k]),Form("X, L%d, %s, E=%3d MeV",i+1,TConfig[j].c_str(),(int)Ene[k]),20,-10,+10);  
	  	  HitPosX[i][j][k]->GetXaxis()->SetTitle("X (cm)");
	  	  HitPosX[i][j][k]->GetYaxis()->SetTitle("#entries");
          HitPosX[i][j][k]->GetYaxis()->SetTitleOffset(2.);
          HitPosX[i][j][k]->SetLineColor(1+k);  
          HitPosX[i][j][k]->SetLineWidth(2);
	  HitPosX[i][j][k]->SetTitle(Form("L%d, %s",i+1,TConfig[j].c_str()));
	 }//k
      }//j
   }//i

 //Energy deposited in scintillators in MeV
 TH1F**TrigEne=new TH1F*[3];
 
 for(int i=0;i<3;i++)
   {
    TrigEne[i]=new TH1F(Form("Energy %s",TConfig[i].c_str()),Form("Energy %s",TConfig[i].c_str()),60,0.,6);
    TrigEne[i]->GetXaxis()->SetTitle("Deposited energy (MeV)");
    TrigEne[i]->GetYaxis()->SetTitle("Events");
    TrigEne[i]->GetYaxis()->SetTitleOffset(2.);
    TrigEne[i]->SetLineColor(kBlack);     
   } 
   
  //Energy deposited per hit in Trackers in MeV
 TH1F**TrckEne=new TH1F*[L];
 
 for(int i=0;i<L;i++)
   {
    TrckEne[i]=new TH1F(Form("Energy L%d",i+1),Form("Energy L%d",i+1),100,0.,0.5);
    TrckEne[i]->GetXaxis()->SetTitle("Deposited energy (MeV)");
    TrckEne[i]->GetYaxis()->SetTitle("Events");
    TrckEne[i]->GetYaxis()->SetTitleOffset(2.);
    TrckEne[i]->SetLineColor(kBlack);     
   }   
   
   
  //Number of simulated event per enegry
  float N=100000.;
   
  //trigger configuration efficiencies
  
  float**NTrig=new float*[NTConfig];

  for(int j=0;j<NTConfig;j++)
    {
     NTrig[j]=new float[E];
     for(int k=0;k<E;k++)
       {
	NTrig[j][k]=0; 
       }//k
    }//j
  
 //Loop over the events
 // This is the main loop
 for (int i=0;i<nentries;i++)
   {
    tree->GetEntry(i); //Load the entry i in the variable e 
    bool* t=new bool[NTConfig];
    for(int j=0;j<NTConfig;j++) t[j]=true;
    if(i%100000==0) cout << "Event: " << i <<endl;
    
    float ene= 1000.*e->get_EkMC();//Energy from Gev to MeV	 
    int ee=0;	 
    for(int k=0;k<E;k++)//Select Energy of the injected particle
      {
       if(ene>=0.9*Ene[k] &&ene<=1.1*Ene[k]) {ee=k;k=E;}
      }
    //No threshold is applied on the signal in the scintillators
    //Should be done here in the future  
    if(!e->get_T1())t[0]=t[3]=t[4]=t[5]=false;
    if(!e->get_T3())t[1]=t[3]=t[5]=false;
    if(!e->get_T4())t[2]=t[4]=t[5]=false;

    //Fill total energy deposited in the scintillators
    int nT=e->get_EneT1().size();
    double tmpE=0;
    for(int j=0;j<nT;j++) tmpE+=1000.*e->get_EneT1().at(j);
    if(tmpE<0.6)t[0]=t[3]=t[4]=t[5]=false;
    if(t[0])TrigEne[0]->Fill(tmpE);
    nT=e->get_EneT3().size();
    tmpE=0;
    for(int j=0;j<nT;j++) tmpE+=1000.*e->get_EneT3().at(j);
    if(tmpE<0.6)t[1]=t[3]=t[5]=false;
    if(t[1])TrigEne[1]->Fill(tmpE);
    nT=e->get_EneT4().size();
    tmpE=0;
    for(int j=0;j<nT;j++) tmpE+=1000.*e->get_EneT4().at(j);
    if(tmpE<1.2)t[2]=t[4]=t[5]=false;
    if(t[2])TrigEne[2]->Fill(tmpE);
   
    //Increment trigger efficiency variables
    for(int k=0;k<NTConfig;k++)
      {
	if(t[k]) NTrig[k][ee]+=100./N;
      }
    //loop over the number of hits
    int nnhits=(int)e->get_Nhits();
    for(int j=0;j<nnhits;j++)
      {
       float X=((e->get_hits().at(j))->get_xin()+(e->get_hits().at(j))->get_xout())/2;//Mean of xin and xout
       int reg=(e->get_hits().at(j))->get_mregMC();
       int l=0;
       for(int k=0;k<L;k++)//Select layer of the hit
         {
	  if(TckReg[k]==reg) {l=k;k=L;}
	 }
	//Fill deposited energy at the hit 
       TrckEne[l]->Fill(1000.*(e->get_hits().at(j))->get_DeltaE());	
       if(1000.*(e->get_hits().at(j))->get_DeltaE()<0.09)continue;
       for(int k=0;k<NTConfig;k++)//Fill the appropriate histograms
         {
          if(t[k])HitPosX[l][k][ee]->Fill(X);
	 }//k
      }//j
 
   }//End loop on event


 //Make Graph of trigger efficiencies  
  TGraph**gTrig=new TGraph*[NTConfig];
  TMultiGraph*multiTrig=new TMultiGraph();
  TLegend*legmulti=new TLegend(0.5,0.2,0.9,0.3);
  legmulti->SetFillStyle(0);legmulti->SetBorderSize(0);legmulti->SetNColumns(2);
  for(int j=0;j<NTConfig;j++)
    {
     gTrig[j]=new TGraph(E,Ene,&NTrig[j][0]);
     gTrig[j]->SetMarkerStyle(20+j);
     gTrig[j]->SetMarkerColor(40+j);
     gTrig[j]->SetMarkerSize(1.5);
     gTrig[j]->SetLineColor(40+j);
     gTrig[j]->SetLineWidth(2);
     gTrig[j]->GetXaxis()->SetTitle("E_{electron} (MeV)");
     gTrig[j]->GetYaxis()->SetTitle("Efficiency (%), No threshold applied in scintillators");  
     gTrig[j]->SetTitle(Form("%s",TConfig[j].c_str())); 
     multiTrig->Add(gTrig[j],"lp");
     legmulti->AddEntry(gTrig[j],Form("%s",TConfig[j].c_str()),"lp");
    }//j   
   
   
////////////////////   
//Display  
////////////////////
 
 //Trackers
 TCanvas**can=new TCanvas*[L];
 
 //One single legend energy dependent
 TLegend*leg=new TLegend(0.4,0.12,0.7,0.40);
 leg->SetFillStyle(0); leg->SetBorderSize(0);
 for(int k=0;k<E;k++)//Select Energy of the injected particle
   {
    leg->AddEntry(HitPosX[0][0][k],Form("E_{e#minus}= %d MeV",(int) Ene[k]),"l");
   }  
 
 gStyle->SetOptStat(0);//Doesn't show the statistical box
 for(int i=0;i<L;i++)
   {
    can[i]=new TCanvas(Form("Layer %d",i+1),Form("Layer %d",i+1),200,10,1500,1000);
    can[i]->Divide(3,2);
    for(int j=0;j<NTConfig;j++)
      {
       can[i]->cd(j+1); 
       gPad->SetLeftMargin(0.15);
       HitPosX[i][j][0]->Draw();
       for(int k=1;k<E;k++)//Select Energy of the injected particle
         {
	  HitPosX[i][j][k]->Draw("same");
	 }
       if(j==0)leg->Draw("same");
       float max=0;
       for(int k=0;k<E;k++)//Determine maximum of entries
         {
	  if(HitPosX[i][j][k]->GetBinContent(HitPosX[i][j][k]->GetMaximumBin())>max)max=HitPosX[i][j][k]->GetBinContent(HitPosX[i][j][k]->GetMaximumBin());
	 }
       for(int k=0;k<E;k++)//Determine maximum of entries
         {
	  HitPosX[i][j][k]->SetMaximum(1.05*max);
	 }  
      }//j
    
   }//k

 //Energy per hits in tracker  
 TCanvas*canTrackEne=new TCanvas("canTrackEne","canTrackEne",200,10,2000,1000); 
 canTrackEne->Divide(4,2);
 for(int j=0;j<L;j++)
   {
    canTrackEne->cd(j+1); 
    gPad->SetLeftMargin(0.15);
    TrckEne[j]->GetYaxis()->SetTitleOffset(2);  
    TrckEne[j]->Draw("hists");
   }     
   
   
//Triggers
 TCanvas*canTrigmulti=new TCanvas("canTrigmulti","canTrigmulti",200,10,800,800); 
 canTrigmulti->cd();
 gPad->SetLeftMargin(0.15);
 multiTrig->Draw("alp");   
 multiTrig->GetXaxis()->SetTitle("E_{electron} (MeV)");
 multiTrig->GetYaxis()->SetTitle("Efficiency (%) - No threshold applied in scintillators");  
 multiTrig->GetYaxis()->SetTitleOffset(1.5);  
 multiTrig->SetMinimum(0.1);
 gPad->SetLogy(1);
 legmulti->Draw("same");  
 gPad->Update();  
 
 TCanvas*canTrig=new TCanvas("canTrig","canTrig",200,10,1500,1000); 
 canTrig->Divide(3,2);
 for(int j=0;j<NTConfig;j++)
   {
    canTrig->cd(j+1); 
    gPad->SetLeftMargin(0.15);
    gTrig[j]->GetYaxis()->SetTitleOffset(2);  
    gTrig[j]->Draw("alp");
   }  
 //Energy in Triggers  
 TCanvas*canTrigEne=new TCanvas("canTrigEne","canTrigEne",200,10,1500,500); 
 canTrigEne->Divide(3,1);
 for(int j=0;j<3;j++)
   {
    canTrigEne->cd(j+1); 
    gPad->SetLeftMargin(0.15);
    TrigEne[j]->GetYaxis()->SetTitleOffset(2);  
    TrigEne[j]->Draw("hists");
   }  
////////////////////
//Save plots into a pdf file
//////////////////  
   
  canTrig->Print("ALTriggerHits_17012017.pdf["); // No actual print, just open the file
  canTrig->Print("ALTriggerHits_17012017.pdf");// Actually print canvas to the file
  canTrigmulti->Print("ALTriggerHits_17012017.pdf");// Actually print canvas to the file
  canTrigEne->Print("ALTriggerHits_17012017.pdf");// Actually print canvas to the file
  canTrackEne->Print("ALTriggerHits_17012017.pdf");// Actually print canvas to the file
  for(int i=0;i<L;i++) can[i]->Print("ALTriggerHits_17012017.pdf");// Actually print canvas to the file
  canTrig->Print("ALTriggerHits_17012017.pdf]");        // No actual print, just close the file   
 
////////////////////
//Close files
////////////////////
 //filein->Close();
 return 1;
}