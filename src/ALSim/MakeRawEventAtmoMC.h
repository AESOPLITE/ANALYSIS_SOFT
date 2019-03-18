////////////////////////////////////////////////////////////////////////////////////////////////////
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, October 28, 2016
///			   Sarah Mechbal, smechbal@ucsc.edu
///	   Santa Cruz Institute for Particle Physics, University of California, Santa Cruz, March 2019 
///////////////////////////////////////////////////////////////////////////////////////////////////


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
#include "LoadMCparameters.h"
#include "tools.h"



int MakeRawEventAtmoMC(int,int,string,string,string,string,string);
