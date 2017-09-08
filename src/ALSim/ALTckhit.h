////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, October 28, 2016
////////////////////////////////////////////////////////////////////////////////////////// 

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <TROOT.h>
#include <Riostream.h>
#include "TChain.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLeaf.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLinearFitter.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TCut.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"
#include "TPDF.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TPostScript.h"
#include "TPaveText.h"
#include "TString.h"
#include "TH1F.h"
#include "TSystem.h"
#include "THStack.h"
#include "TNtuple.h"
#include "TFormula.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TH3F.h"
#include "TRandom3.h"
#include "TRandom2.h"
#include "TRandom1.h"



#ifndef ALTckhit_H
#define ALTckhit_H

using namespace std;

class ALTckhit:public TObject
{
 private:
 
  float mregMC; // region if MC
  float mtrackMC; // mtrack of track
  int typeMC;   //type of particle if MC
  float eMC;     //total energy if MC (entrance of the tracker layer)
  float xin;      //x coordinate of hit MC (entrance of the tracker layer)
  float yin;      //y coordinate of hit MC (entrance of the tracker layer)
  float zin;      //z coordinate of hit MC (entrance of the tracker layer)
  float xout;     //x coordinate of hit MC (exit the tracker layer)
  float yout;     //y coordinate of hit MC (exit the tracker layer)
  float zout;     //z coordinate of hit MC (exit the tracker layer)
  float age;    //Age of particle if MC  (entrance of the tracker layer)
  float cx;     //cosineX of momentum MC (Average in the tracker layer)
  float cy;     //cosineY of momentum MC (Average in the tracker layer)
  float cz;     //cosineZ of momentum MC (Average in the tracker layer)
  int flag;     //flag value: Track=1, Hit=0
  float DeltaE;     //Energy along the track or at hit location Hit=0  
  
  
  //Raw Data from cluster information
  int L;         //Layer from 0 to 6 top to bottom
  int chip;      //Chip ID: 0 to 11
  int nstrips;   //Number of strips in the cluster
  int fstrip;    //First strip ID from 0 to 63 
  int fstripID;  //First strip on the layer from 0 767
  //Coordinates of the cluster in cm determined from the raw data
  float x;
  float y;
  float z;

  
  //Information from Pattern Recognition
 
  float xPR;       //x coordinate of "chosen" hit MC (selected by Pattern Recognition)
  float yPR;       //y coordinate of "chosen" hit MC (selected by Pattern Recognition)
  float zPR;       //z coordinate of "chosen" hit MC (selected by Pattern Recognition)
  bool flagPR;	   // flag  if hit chosen for reconstruction, flag = 0 if not 
 
  //Reconstructed information
  float xreco;      //x coordinate of hit 
  float yreco;      //y coordinate of hit 
  float zreco;      //z coordinate of hit 
  float agereco;    //Age of particle  
  float cxreco;     //cosineX of momentum 
  float cyreco;     //cosineY of momentum 
  float czreco;     //cosineZ of momentum    
  float ereco;     //kinetic energy 
  int k;     	   //kth hit in event  
 
 public: 
   //Constructors
   ALTckhit();// Default
   //Destructor
   ~ALTckhit(){};  
     ////////////////////////////////
   //"setting" member methods
   ////////////////////////////////
   void set_mregMC(float a){mregMC=a;}
   void set_mtrackMC(float a){mtrackMC=a;}
   void set_typeMC(int a){typeMC=a;}
   void set_eMC(float a){eMC=a;}
   void set_xin(float a){xin=a;}
   void set_yin(float a){yin=a;}
   void set_zin(float a){zin=a;}
   void set_xout(float a){xout=a;}
   void set_yout(float a){yout=a;}
   void set_zout(float a){zout=a;}
   void set_age(float a){age=a;}
   void set_cx(float a){cx=a;}
   void set_cy(float a){cy=a;}
   void set_cz(float a){cz=a;}
   void set_flag(int a){flag=a;}
   void set_DeltaE(float a){DeltaE=a;}
   ////////////////////////////////
   void set_xreco(float a){xreco=a;}
   void set_yreco(float a){yreco=a;}
   void set_zreco(float a){zreco=a;}
   void set_cxreco(float a){cxreco=a;}
   void set_cyreco(float a){cyreco=a;}
   void set_czreco(float a){czreco=a;}
   void set_ereco(float a){ereco=a;}
   void set_agereco(float a){agereco=a;}
   void set_k(int a){k=a;}
   ////////////////////////////////
   void set_L(int a){L=a;}         
   void set_chip(int a){chip=a;}      
   void set_nstrips(int a){nstrips=a;}   
   void set_fstrip(int a){fstrip=a;}    
   void set_fstripID(int a){fstripID=a;}  
   void set_x(float a){x=a;}
   void set_y(float a){y=a;}
   void set_z(float a){z=a;}    
   /////////////////////////////////
   
   void set_xPR(float a){xPR=a;}
   void set_yPR(float a){yPR=a;}
   void set_zPR(float a){zPR=a;}
   void set_flagPR(bool a){flagPR=a;}
  
   ////////////////////////////////
   //"Getting" member methods
   ////////////////////////////////
   float get_mregMC( ){return mregMC;}
   float get_mtrackMC( ){return mtrackMC;}
   int get_typeMC( ){return typeMC;}
   float get_eMC( ){return eMC;}
   float get_xin( ){return xin;}
   float get_yin( ){return yin;}
   float get_zin( ){return zin;}
   float get_xout( ){return xout;}
   float get_yout( ){return yout;}
   float get_zout( ){return zout;}
   float get_cx( ){return cx;}
   float get_cy( ){return cy;}
   float get_cz( ){return cz;}
   int get_flag( ){return flag;}
   float get_DeltaE( ){return DeltaE;}
     
   ////////////////////////////////
   int get_L(){return L;}         
   int get_chip(){return chip;}      
   int get_nstrips(){return nstrips;}   
   int get_fstrip(){return fstrip;}    
   int get_fstripID(){return fstripID;}  
   float get_x(){return x;}
   float get_y(){return y;}
   float get_z(){return z;}    
   ////////////////////////////////

   float get_xPR() {return xPR;}
   float get_yPR() {return yPR;}
   float get_zPR() {return zPR;}
   bool  get_flagPR() {return flagPR;}
   
   ////////////////////////////////
   float get_xreco( ){return xreco;}
   float get_yreco( ){return yreco;}
   float get_zreco( ){return zreco;}
   float get_cxreco( ){return cxreco;}
   float get_cyreco( ){return cyreco;}
   float get_czreco( ){return czreco;}
   float get_ereco( ){return ereco;}
   float get_agereco( ){return agereco;}
   int get_k( ){return k;}

   ////////////////////////////////
   ClassDef(ALTckhit,1)
};


#endif
