////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, October 28, 2016
////////////////////////////////////////////////////////////////////////////////////////// 

void write99toroot()
{

 TFile*file=new TFile("aesoplite.root","RECREATE");
 //TNtuple*ntuple=new TNtuple("BoundX","BoundX","ncase:mreg:newreg:jtrack:ek:xsco:ysco:zsco:atrack:cxtrck:cytrck:cztrck");
 TNtuple*ntuple=new TNtuple("Track","Track","ncase:mreg:mtrack:type:age:e:x:y:z:cx:cy:cz:Edep:flag");
 ntuple->ReadFile("aesoplite001_fort.99");
 
 //Can read multiple files. Ouputs will be merge in an single file
 // ntuple->ReadFile("aesoplite002_fort.99");
 // ntuple->ReadFile("aesoplite003_fort.99");

 ntuple->Write();
 file->Close();
  
}