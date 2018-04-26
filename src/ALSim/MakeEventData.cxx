////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, September 11, 2017
////////////////////////////////////////////////////////////////////////////////////////// 

#include "MakeEventData.h"
#include "MakeRawBPDEvent.h"
#include "ALPatternRecognition.h"
#include "ALKalman.h"
#include "LoadDataparameters.h"
#include "TBox.h"
#include "headers.h"

int MakeEventData(string filename,int geoconfig, int FieldConf, bool TwoIter,string RecoID)
{
 //Load configuration parameter
 float* zL=new float[7];
 float*OffsetLL=new float[7];
 float*OffsetRL=new float[7];
 float*TrigThresh=new float[5];
 for(int i=0;i<7;i++)zL[i]=OffsetLL[i]=OffsetRL[i]=0;
 for(int i=0;i<5;i++)TrigThresh[i]=0;
 string paramfile=Form("../src/ALSim/Dataparameters%d.dat",geoconfig); 

 LoadDataparameters(paramfile,zL,OffsetLL,OffsetRL,TrigThresh);

//Load magnetic field map 
//Set Magnetic field map
// TFile*file_map=new TFile("../prod/fieldmap.root","READ");
//load 1mm grid magnetic field map
 TFile*file_map=new TFile("/home/sarah/AESOPLITE/ANALYSIS_SOFT/prod/fieldmap1mm.root","READ");
TBField *bfield = new TBField();

 bool FlagMagF=false;  
 cout << "Field Conf bool set to " << FieldConf << endl;

 bfield->SetUseUniformBfield(FieldConf);
 

 FlagMagF=bfield->SetMagField();
 if(FlagMagF)cout << "Field Map loaded" << endl;
 else 
  {
   cout << "There is an issue when loading the Field Map" << endl;
   return 1;
  }
 
	
 for(int i=0;i<7;i++)
   {
    cout << "L"<<i <<", zL:" << zL[i] ;
    cout << ", OffsetLL:" << OffsetLL[i] ;
    cout << ", OffsetRL:" << OffsetRL[i] << endl;
   }  
 cout << "T1 threshold: " << TrigThresh[0] <<endl;
 cout << "T2 threshold: " << TrigThresh[1] <<endl;
 cout << "T3 threshold: " << TrigThresh[2] <<endl;
 cout << "T4 threshold: " << TrigThresh[3] <<endl;
 cout << "Guard threshold: " << TrigThresh[4] <<endl;


 //Maximum of configuration to try fitting
 
 int MAXB=1000;
 int MAXNB=1000;

 //Strip pitch in cm
 float  strippitch=0.0228;

 float zz0=0;//in cm
 
 TCanvas*can;
 TMultiGraph*multi;
 TMultiGraph*multiB;
 TF1*tmpf;
 TF1*tmpfB;
 TF1*tmpfB2;


 //Input root file 
 //Get input filename and ASIC Length and EVT Length
 vector<string> input;
 //Split the line 
 input=split(&filename,' '); 
  
 TFile*file=new TFile(Form("%s.root",input.at(0).c_str()),"READ");
 cout << "Input file is open" <<endl;

 //Input root file 
 TFile*fileout=new TFile(Form("%s.EVENT_%s.root",input.at(0).c_str(),RecoID.c_str()),"RECREATE");
 cout << "Output file is created" <<endl;
 
 
 //Get Tree from the input file
 TTree *tree = (TTree*)file->Get("BPD");
 //Define variables to read event
 ALEvent *e = new ALEvent();
 //Set address to access event data
 tree->SetBranchAddress("event",&e);  
 // Create a TTree
 TTree *DEtree = new TTree("Data"," Data Event");
 //Define variables to make the Reco event
 ALEvent *de=new ALEvent();
 // Create a branch with event
 DEtree->Branch("event",&de); 
 // Get number of events in Tree
 int nentries=tree->GetEntries();
 cout << "Number  of events: " << nentries << endl;
 //Loop over the events
 // This is the main loop
 
 for (int k=0;k<nentries;k++)
   {
    tree->GetEntry(k); //Load the entry k in the variable e  
    //cout << "Make Event Data: "<<k <<endl;
    //Copy the raw event into a Data event with same structure 
    //if(k>2){k==nentries;continue;}
    if(k%100==0) cout << "Event " << k << endl;
	de=new ALEvent();
    de->Copy(e);  
     
    //Set boolean flag for triggers
    
    if(de->get_EneT1().at(0)>TrigThresh[0]) de->set_T1(true);
    if(de->get_EneT2().at(0)>TrigThresh[1]) de->set_T2(true);
    if(de->get_EneT3().at(0)>TrigThresh[2]) de->set_T3(true);
    if(de->get_EneT4().at(0)>TrigThresh[3]) de->set_T4(true);
    if(de->get_Eneg().at(0)>TrigThresh[4]) de->set_guard(true);

     uint8_t Ti=(uint8_t)de->get_Ti();
    //Number of layers wih hit(s)
    int NL=de->get_NLayers();
    //if(NL>7) cout << "ERROR ON NUMBER OF LAYERS" <<endl;
    int* Lay=new int[7];
    de->get_Layers(Lay);

    int NLB=0;
    int NLNB=0;
    //B layer with hits
    NLB+=(int)((Ti >>1) & 0x01);
    NLB+=(int)((Ti >>2) & 0x01);
    NLB+=(int)((Ti >>3) & 0x01);
    NLB+=(int)((Ti >>5) & 0x01);
    //NB layer with hits
    NLNB+=(int)((Ti >>0) & 0x01);
    NLNB+=(int)((Ti >>4) & 0x01);
    NLNB+=(int)((Ti >>6) & 0x01);

    //loop over the number of clusters
    int nnhits = (int)de->get_Nhits();

    for(int i=0;i<nnhits;i++)
      {
       float x=-999.;
       float y=-999.;
       float z=-999.;
       float offsetLL=0;//Left ladder offset
       float offsetRL=0;//Right ladder offset
       float center=0;
       int L=((de->get_hits()).at(i))->get_L();//0 to 6
       de->get_hits().at(i)->set_k(i);
       int nstrips=((de->get_hits()).at(i))->get_nstrips(); //0 if only one strip
       int fstripID=((de->get_hits()).at(i))->get_fstripID();
 
       //Determine the coordinate x,y,z parameters
       //Get the values from configuration file read at the beginning of the function
       z=zL[L];
       offsetLL=OffsetLL[L];
       offsetRL=OffsetRL[L];
         
       //Determine the coordinates from the strips
       //Equations from Sarah's email of September 4 2017.
       //Just calculate the mean value of the strip positions
       //Assume that the center of the coodinnates is 0 is equidistant from the 
       //two ladders
       center=(offsetRL+(offsetLL+384.*strippitch))/2.;
       //In non-bending plane: Coordinate X in cm
       if(L==0||L==4||L==6)
        {
         //First 6 chips: 0 to 5; strip number 1 to 384
         if(fstripID>0 &&fstripID<=384)
          {
           x=(fstripID-0.5)*strippitch+offsetLL-center; 
          }
         //Last 6 chips: 6 to 11; strip number 385 to 768
         if(fstripID>384)
          {
           x=(fstripID-384.5)*strippitch+offsetRL-center;
          }
          
         //Calculate the center of the cluster depending on the number of strips 
    
         x+=nstrips*strippitch/2.;             
         y=-999.;          
        }
       
       //In bending plane: Coordinate Y in cm
       if(L==1||L==2||L==3||L==5)
        {
         //First 6 chips: 0 to 5; strip number 1 to 384
         if(fstripID>0 &&fstripID<=384)
          {
           y=(fstripID-0.5)*strippitch+offsetLL-center; 
       
          }
         //Last 6 chips: 6 to 11; strip number 385 to 768
         if(fstripID>384)
          {
           y=(fstripID-384.5)*strippitch+offsetRL-center;
       
          }
         y+=nstrips*strippitch/2.;    
 
         x=-999.; 
        }

        
       //Fill the coordinates    
       ((de->get_hits()).at(i))->set_x(x);
       ((de->get_hits()).at(i))->set_y(y);
       ((de->get_hits()).at(i))->set_z(z);          
      }  //i nnhits

    ////////////////////////////////////
    //TRIGGER
    ////////////////////////////////////  
      
    //if not less than 5 layers were touched then don't try pattern recognition
    //cout << "Internal trigger: " << de->get_Ti() <<endl;
    if(NL<5)
     {   
      /////////////////////    
      //Fill the output file 
      /////////////////////
      DEtree->Fill();
      //Free memory    
      delete de;
      continue;
    }
    //T1&T3&T4
   /* if(de->get_EneT1().at(0)<0 || de->get_EneT3().at(0)<0 ||de->get_EneT4().at(0)<0)
     {   
      /////////////////////    
      //Fill the output file 
      /////////////////////
      DEtree->Fill();
      //Free memory    
      delete de;
      continue;
    }
    */
      
    ////////////////////////////////////
    //Pattern Recognition     Sarah's Code
    ////////////////////////////////////
   int DataType = 1; 
   ALPatternRecognition* PR = new ALPatternRecognition();
   int EventPR = PR->FindPattern(de,DataType);
   if(EventPR==0) {
         DEtree->Fill();
      //Free memory    
      delete de;
      continue;
      }
   
//  cout << "Event " << k << " deflecPR = " << de->get_deflecPR() << endl;
    /////////////////////    
    //RECONSTRUCTION
    /////////////////////    
     
    ALKalman* KF = new ALKalman(de);
    int EventKF=  KF->DoKF(de, DataType, 1, TwoIter);
    if(EventKF==0) {		   // delete KF;
	DEtree->Fill();
	//Free memory
	delete de;
	continue;
	}
    
    /////////////////////    
    //Fill the output file 
    /////////////////////
    DEtree->Fill();
    //Free memory    
    delete de;
   }//k End loop on events
 delete e;
 //Write tree in output file
 fileout->cd();
 DEtree->Write(0,TObject::kOverwrite);
 //Close files
 fileout->Close();
 file->Close();
      
 return 0;
}



