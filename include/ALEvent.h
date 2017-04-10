////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, October 28, 2016
////////////////////////////////////////////////////////////////////////////////////////// 

#include "ALTckhit.h"

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


using namespace std;

#ifndef ALEVENT_H
#define ALEVENT_H


class ALEvent:public TObject
{
 private:
  
   int eventnumber; //Event number
   //Monte Carlo information: Truth variable names finish with MC 
   int ncase; 
   int typeMC; //type of particle
   double EkMC;   //kinetic energy of the particle
   double X0MC,Y0MC,Z0MC;//Coordinates of the partcle at the injection point 
   double CX0MC,CY0MC,CZ0MC; //Incidence cosines of the partcle at the injection point 
   
   // Hits information
   int Nhits; //Number of hits in the event
   
   //Reconstruction information: variables finish with reco
   int typereco; //type of particle
   double Ekreco;   //kinetic energy of the particle
   double X0reco,Y0reco,Z0reco;//Coordinates of the partcle at the injection point 
   double CX0reco,CY0reco,CZ0reco; //Incidence cosines of the partcle at the injection point 
   
   //Hits information
   std::vector<ALTckhit*> hits;  
   
   //Triggers
   bool T1;
   bool T2;
   bool T3;
   bool T4;
   bool guard;
   std::vector<double> EneT1;  
   std::vector<double> EneT2;  
   std::vector<double> EneT3;
   std::vector<double> EneT4;
   std::vector<double> Eneg;   
   std::vector<double> timeT1;  
   std::vector<double> timeT2;  
   std::vector<double> timeT3;
   std::vector<double> timeT4;
   std::vector<double> timeg;     
 public: 
   //Constructors
   ALEvent();// Default
   //Destructor
   ~ALEvent(){};   
   ////////////////////////////////
   //"setting" member methods
   ////////////////////////////////
   void set_eventnumber(int a){eventnumber=a;}
   void set_ncase(int a){ncase=a;}
   void set_typeMC(int a){typeMC=a;}
   void set_EkMC(double a){EkMC=a;}
   void set_X0MC(double a){X0MC=a;}
   void set_Y0MC(double a){Y0MC=a;}
   void set_Z0MC(double a){Z0MC=a;}
   void set_CX0MC(double a){CX0MC=a;}
   void set_CY0MC(double a){CY0MC=a;}
   void set_CZ0MC(double a){CZ0MC=a;}
   ////////////////////////////////
   void set_Nhits(int a){Nhits=a;}
   void add_Nhits(){Nhits++;}
   ////////////////////////////////
   void set_typereco(int a){typereco=a;}
   void set_Ekreco(double a){Ekreco=a;}
   void set_X0reco(double a){X0reco=a;}
   void set_Y0reco(double a){Y0reco=a;}
   void set_Z0reco(double a){Z0reco=a;}
   void set_CX0reco(double a){CX0reco=a;}
   void set_CY0reco(double a){CY0reco=a;}
   void set_CZ0reco(double a){CZ0reco=a;}
   void add_hit(ALTckhit* h){hits.push_back(h);Nhits++;}
   void set_T1(bool a){T1=a;}
   void set_T2(bool a){T2=a;}
   void set_T3(bool a){T3=a;}
   void set_T4(bool a){T4=a;}
   void set_guard(bool a){guard=a;}
   void add_EneT1(double a){EneT1.push_back(a);}
   void add_EneT2(double a){EneT2.push_back(a);}
   void add_EneT3(double a){EneT3.push_back(a);}
   void add_EneT4(double a){EneT4.push_back(a);}
   void add_Eneg(double a){Eneg.push_back(a);}
   void add_timeT1(double a){timeT1.push_back(a);}
   void add_timeT3(double a){timeT3.push_back(a);}
   void add_timeT4(double a){timeT4.push_back(a);}
   void add_timeg(double a){timeg.push_back(a);}   
   
   ////////////////////////////////
   //"Getting" member methods
   ////////////////////////////////
   int get_eventnumber(){return eventnumber;}
   int get_ncase(){return ncase;}
   int get_typeMC(){return typeMC;}
   double get_EkMC(){return EkMC;}
   double get_X0MC(){return X0MC;}
   double get_Y0MC(){return Y0MC;}
   double get_Z0MC(){return Z0MC;}
   double get_CX0MC(){return CX0MC;}
   double get_CY0MC(){return CY0MC;}
   double get_CZ0MC(){return CZ0MC;}
   ////////////////////////////////
   int get_Nhits(){return Nhits;}
   ////////////////////////////////
   int get_typereco(){return typereco;}
   double get_Ekreco(){return Ekreco;}
   double get_X0reco(){return X0reco;}
   double get_Y0reco(){return Y0reco;}
   double get_Z0reco(){return Z0reco;}
   double get_CX0reco(){return CX0reco;}
   double get_CY0reco(){return CY0reco;}
   double get_CZ0reco(){return CZ0reco;}
   std::vector<ALTckhit*> get_hits(){return hits;}
   bool get_T1(){return T1;}
   bool get_T2(){return T2;}
   bool get_T3(){return T3;}
   bool get_T4(){return T4;}
   bool get_guard(){return guard;}
   std::vector<double>  get_EneT1(){return EneT1;}
   std::vector<double>  get_EneT2(){return EneT2;}
   std::vector<double>  get_EneT3(){return EneT3;}
   std::vector<double>  get_EneT4(){return EneT4;}
   std::vector<double>  get_Eneg(){return Eneg;}
   std::vector<double>  get_timeT1(){return timeT1;}
   std::vector<double>  get_timeT2(){return timeT2;}
   std::vector<double>  get_timeT3(){return timeT3;}
   std::vector<double>  get_timeT4(){return timeT4;}
   std::vector<double>  get_timeg(){return timeg;}

   ////////////////////////////////
   ClassDef(ALEvent,1)

};

#endif