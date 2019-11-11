//////////////////////////////////////////////////////////////////////////////////////////
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, September 8, 2017
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
#include "ALEvent.h"

//ClassImp(ALTckhit)
//ClassImp(ALEvent)



vector<string> split (string*,char);
double s2lf(string*);
float s2f(string*);
int s2i(string*);
int Partition(int ,int ,float []);
void Quick_sort(int ,int ,float []);
void extractdate(int *y,int*m,int*d,string*str);
void extracttime(int *h,int*m,int*s,string*str);
int DecodeASIShort(string data,vector<ALTckhit*>* Hh,int*);
int DecodeASILong(string data,vector<ALTckhit*>* Hh,int*,int*,int *);
int MakeRawBPDEvent(string);
int MakeRawBPDEventIT(string);
