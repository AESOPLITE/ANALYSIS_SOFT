////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, October 28, 2016
////////////////////////////////////////////////////////////////////////////////////////// 

void write99toroot()
{

 char* Outpath="rootfiles";
 char* Inppath=".";

	//////////////for electrons/positrons///////// 
	
 int type=3; //3 for electrons, 4 for positrons
 int Nene=12;
 int Ene[12]={10,20,30,40,50,60,70,80,90,100,200,300};
 
/*
	////////////for mu+-////////////
	
 int type=11;					//10 for mu+
 int Nene=3;
 int Ene[3]={1,2,3};
*/	
///////////////////////////////////////////////////////	
 int* Ncycles=new int[Nene]; 
 for(int i=0;i<Nene;i++)Ncycles[i]=100;  
 
 TFile*file;
 TNtuple*ntuple;
 
 char* startfile="aesopliteNonUniB_V1";
 char* endfile="_fort.99";

 for(int i=3;i<4;i++)
   {
   // for(int j=0;j<Ncycles[i];j++)
	   for(int j=0;j<10;j++)
      {
		  
	////////for electrons/////////////////////////
		  
       cout << Form("%s/%d/%s_%d_%dMeV%03d%s",Inppath,type,startfile,type,Ene[i],j+1,endfile) <<endl;
       file=new TFile(Form("%s/%d/%s/%s_%d_%dMeV%03d%s.root",Inppath,type,Outpath,startfile,type,Ene[i],j+1,endfile),"RECREATE");
       ntuple =new TNtuple("Track","Track","ncase:mreg:mtrack:type:age:p:x:y:z:cx:cy:cz:Edep:flag");
       ntuple->ReadFile(Form("%s/%d/%s_%d_%dMeV%03d%s",Inppath,type,startfile,type,Ene[i],j+1,endfile));
		  
	////////////for muons////////////////////////////	  
		/*  
	   cout << Form("%s/%d/%s_%d_%dGeV%03d%s",Inppath,type,startfile,type,Ene[i],j+1,endfile) <<endl;
       file=new TFile(Form("%s/%d/%s/%s_%d_%dGeV%03d%s.root",Inppath,type,Outpath,startfile,type,Ene[i],j+1,endfile),"RECREATE");
       ntuple =new TNtuple("Track","Track","ncase:mreg:mtrack:type:age:e:x:y:z:cx:cy:cz:Edep:flag");
       ntuple->ReadFile(Form("%s/%d/%s_%d_%dGeV%03d%s",Inppath,type,startfile,type,Ene[i],j+1,endfile));
		  
		*/  
       ntuple->Write();
       file->Close();   
      }//j
    }//i
}
