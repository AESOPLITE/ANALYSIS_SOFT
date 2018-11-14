////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, October 28, 2016
////////////////////////////////////////////////////////////////////////////////////////// 

#include "MakeRawEventMC.h"
int MakeRawEventMC(int typeT,int Ene,int cycle,string Inppath,string Inppath2,string Outpath,string startfile,string endfile)
{

cout << "Calling MakeRawEventMC" << endl;
//Load region numbers used in the MC geometry 
 int*TckReg=new int[7];
 int*TrigReg=new int[4];
 int*GReg=new int[1];
 int*ShellReg=new int[2];
 float*TckZPos=new float[7];
 float*TrigThresh=new float[4];
 float*GuardThresh=new float[1];
 for(int i=0;i<7;i++)TckReg[i]=0;
 for(int i=0;i<4;i++)TrigReg[i]=0;
 for(int i=0;i<2;i++)ShellReg[i]=0;
 for(int i=0;i<1;i++)GReg[i]=0;
 for(int i=0;i<7;i++)TckZPos[i]=0;
 for(int i=0;i<4;i++)TrigThresh[i]=0;
 for(int i=0;i<1;i++)GuardThresh[i]=0;
	
 string MCparamfile="../src/ALSim/MCparameters.dat"; 
 
 LoadMCparameters(MCparamfile,TckReg,TrigReg,GReg,TckZPos,TrigThresh,GuardThresh,ShellReg);

 TFile *file;
 TFile *fileout;
 //Input file 
 if(typeT==3 || typeT==1)
  {
   file=new TFile(Form("%s/%d/%s/%s_%d_%dMeV%03d%s.root",Inppath.c_str(),typeT,Inppath2.c_str(),startfile.c_str(),typeT,Ene,cycle,endfile.c_str()),"READ");
   cout << "Input file is open" <<endl;
   //Output file 
   fileout=new TFile(Form("%s/%d/RawEvent_%s_%d_%dMeV%03d%s.root",Outpath.c_str(),typeT,startfile.c_str(),typeT,Ene,cycle,endfile.c_str()),"RECREATE");
   cout << "Output file is created" <<endl;
  }
 else
  {
   file=new TFile(Form("%s/%d/%s/%s_%d_%dGeV%03d%s.root",Inppath.c_str(),typeT,Inppath2.c_str(),startfile.c_str(),typeT,Ene,cycle,endfile.c_str()),"READ");
   cout << "Input file is open" << endl;
   fileout=new TFile(Form("%s/%d/RawEvent_%s_%d_%dGeV%03d%s.root",Outpath.c_str(),typeT,startfile.c_str(),typeT,Ene,cycle,endfile.c_str()),"RECREATE"); 
   cout << "Output file is created " << endl;
   } 

	
//Get ntuple from the input file
 TNtuple*ntuple=(TNtuple*)file->Get("Track");
 cout << "Got the ntuple from the input file" <<endl;

 //Define variables to read ntuple
 float ncase=0; 
 float mreg=0;
 float mtrack=0;
 float type=0;
 float EkMC=0;
 float pMC=0; 
 float x=0;
 float y=0;
 float z=0;
 float age=0; 
 float cx=0;
 float cy=0;
 float cz=0; 
 float Edep=0; 
 float flag=0; 

 //Set addresses to access ntuple data
 
 ntuple->SetBranchAddress("ncase",&ncase); 
 ntuple->SetBranchAddress("mreg",&mreg);
 ntuple->SetBranchAddress("mtrack",&mtrack);
 ntuple->SetBranchAddress("type",&type);
 ntuple->SetBranchAddress("p",&pMC); 
 ntuple->SetBranchAddress("x",&x);
 ntuple->SetBranchAddress("y",&y);
 ntuple->SetBranchAddress("z",&z);
 ntuple->SetBranchAddress("age",&age); 
 ntuple->SetBranchAddress("cx",&cx);
 ntuple->SetBranchAddress("cy",&cy);
 ntuple->SetBranchAddress("cz",&cz); 
 ntuple->SetBranchAddress("flag",&flag); 
 ntuple->SetBranchAddress("Edep",&Edep); 

 
 cout << "Set Branch Address done" <<endl;

 
 // Get number of entries in ntuple (number of lines in the output fluka file) 
 int nentries=ntuple->GetEntries();
 cout << "Number  of entries: " << nentries << endl;
 int ievt=-1; // Identify the event
 // Create a TTree
 TTree *tree = new TTree("MC","Raw event MC");
 ALEvent *e = new ALEvent();
 // Create a branch with event
 tree->Branch("event",&e);  

 int j=0;
 int nh=0;
 double EneFoam=0;
 double EneShell=0;
 double EneT1=0;
 double EneT2=0; 
 double EneT3=0; 
 double EneT4=0; 
 double Eneg=0; 
 double timeFoam=0;
 double timeShell=0;
 double timeT1=0;
 double timeT2=0; 
 double timeT3=0; 
 double timeT4=0;
 int prevreg=0; 
 double timeg=0; 
 double nFoam=0;				//number of segments of track in insulation foam
 double nShell=0;				//number of segments of track in shell
 double nT1=0; 						//number of segments of track in T1
 double nT2=0;
 double nT3=0; 						
 double nT4=0; 
 double ng=0; 
 int iFoam=0;
 int iShell=0;
 int iT1=0; 
 int iT2=0;
 int iT3=0; 
 int iT4=0; 
 int ig=0;
//Temporary variable for internal triggers 
 int*Titmp=new int[7];
 for (int i=0;i<7;i++)Titmp[i]=0;
 
 int nL=0;
 int iL=0;
 ALTckhit* h=new ALTckhit();

 
 //Loop over the entries and create the events
 // This is the main loop
 for (int i=0;i<nentries;i++)
   {
    ntuple->GetEntry(i); //Load the entry i in the variables  
    if(i%100000==0) cout << "entry: " << i <<endl;

    ////////////////////////////////// 
    // Select the types of particle that produce hits
    ////////////////////////////////// 
    if(prevreg==ShellReg[0] && prevreg!=mreg)//Foam insulation
     {
      if(nFoam>0)
       {
        e->add_EneIsofoam(EneFoam); e->add_timeIsofoam(timeFoam);
		EneFoam=0;timeFoam=0;nFoam=0;iFoam=0;
       }
     }	 
    if(prevreg==ShellReg[1] && prevreg!=mreg)//Aluminium Shell
     {
      if(nShell>0)
       {
        e->add_EneShell(EneShell); e->add_timeShell(timeShell);
		EneShell=0;timeShell=0;nShell=0;iShell=0;
       }
     }	
	 
    if(prevreg==TrigReg[0] && prevreg!=mreg)//T1
     {
      if(nT1>0)
       {
        e->add_EneT1(EneT1); e->add_timeT1(timeT1);
		if (EneT1 > TrigThresh[0]) e->set_T1(true);
		EneT1=0;timeT1=0;nT1=0;iT1=0;
       }
     }
    if(prevreg==TrigReg[1] && prevreg!=mreg)//T2
     {
      if(nT2>0)
       {
        e->add_EneT2(EneT2);
        // e->add_timeT2(timeT2);
        if (EneT2 > TrigThresh[1]) e->set_T2(true);
        EneT2=0;timeT2=0;nT2=0;iT2=0;
       }
     }

    if(prevreg==TrigReg[2] && prevreg!=mreg)//T3
     {
      if(nT3>0)
       {
        e->add_EneT3(EneT3); e->add_timeT3(timeT3);
		if(EneT3 > TrigThresh[2]) e->set_T3(true);
        EneT3=0;timeT3=0;nT3=0;iT3=0;
       }
     }
    if(prevreg==TrigReg[3] && prevreg!=mreg)//T4
     {
      if(nT4>0)
       {
        e->add_EneT4(EneT4); e->add_timeT4(timeT4);
		if (EneT4 > TrigThresh[3]) e->set_T4(true);
        EneT4=0;timeT4=0;nT4=0;iT4=0;
       }
     }
    if(prevreg==GReg[0] && prevreg!=mreg)//guard
     {
      if(ng>0)
       {
        e->add_Eneg(Eneg); e->add_timeg(timeg);
		if (Eneg > TrigThresh[4]) e->set_guard(true);
	    Eneg=0;timeg=0;ng=0;ig=0;
       }
     }
    if(prevreg>=TckReg[0] && prevreg<=TckReg[6]&& prevreg!=mreg)
     {       
      if(nL>0){h->set_cx(h->get_cx()/nL);h->set_cy(h->get_cy()/nL);h->set_cz(h->get_cz()/nL);}
      //End of the hit. Add the hit to the event
      h->set_k(nh);
      e->add_hit(h);
      nh++;
      //Check the layer for internal trigger
      for (int ij=0;ij<7;ij++)
        {
	 if (h->get_mregMC()== TckReg[ij]) {Titmp[ij]=1;ij=7;}
	}
      //Reset the hit structure
      h=new ALTckhit();
      nL=0;
      iL=0;
     }  
    if(ncase!=ievt)//Start a new event
     {
      j++;
      nh=0;

      if(i==0)//if first event
       {
	ievt=ncase; 
	delete e;
	e=new ALEvent();
        e->set_eventnumber(j);
        e->set_ncase(ncase);
	//First line  with a given ncase contains the source particle info
        e->set_typeMC(type);
        e->set_pMC(pMC);
   //Calculate kinetic energy from momentum at point of injection
		float mass=0.000511;								//electron mass in GeV
        if(type==11 || type ==10)	mass=0.10566;			//muon mass in GeV
	    EkMC = TMath::Sqrt(pMC*pMC+mass*mass);
		e->set_EkMC(EkMC);
                e->set_X0MC(x);
                e->set_Y0MC(y);
                e->set_Z0MC(z);
        e->set_CX0MC(cx);
        e->set_CY0MC(cy);
        e->set_CZ0MC(cz);
	//First attempt to fill T2
      //  e->set_T2(false);
       }	//if
      else
       {
	//Fill internal trigger information of the previous event
	int t=0;
	for (int ij=0;ij<7;ij++) t+=Titmp[ij]*(int)TMath::Power(2,ij);
	e->set_Ti(t);  
        //Fill the previous event in the output file  
	tree->Fill(); 
	delete e;
    e=new ALEvent();
	h=new ALTckhit();
	nh=0;
	ievt=ncase; 
        e->set_eventnumber(j);
        e->set_ncase(ncase);
        //First line with a given ncase contains the source particle info
        e->set_typeMC(type);
        e->set_pMC(pMC);
   //Calculate kinetic energy from momentum at point of injection
		float mass=0.000511;								//electron mass in GeV
        if(type==11 || type ==10)	mass=0.10566;			//muon mass in GeV
	    float EkMC = TMath::Sqrt(pMC*pMC+mass*mass);
		e->set_EkMC(EkMC);
        e->set_X0MC(x);
        e->set_Y0MC(y);
        e->set_Z0MC(z);
        e->set_CX0MC(cx);
        e->set_CY0MC(cy);
        e->set_CZ0MC(cz);
        //Test to set T2
       // e->set_T2(false);
        //Reset internal trigger
        for (int ij=0;ij<7;ij++)Titmp[ij]=0;
       }     //else
     }  //if
    
    ////////////////////////////////// 
    //Check to see if the particle has crossed one of the scintillator or shell
    //(r19=air, r1=T1, r6=T3, r7=guard, r11=Tracker1 ,..., r17=Tracker7, r18=T4)
    ////////////////////////////////// 
   
	//Shell & Insulation
    if (mreg == ShellReg[0])//Insulating foam
      {
       if(i==iFoam+1){EneFoam+=Edep;nFoam++;}
       else{nFoam=1;EneFoam=Edep;timeFoam=age;}
       iFoam=i;
      }
    if (mreg == ShellReg[1])//Aluminium Shell
      {
       if(i==iShell+1){EneShell+=Edep;nShell++;}
       else{nShell=1;EneShell=Edep;timeShell=age;}
       iShell=i;
      }	 
	 
    //Triggers
    if (mreg == TrigReg[0])//T1
      {
       if(i==iT1+1){EneT1+=Edep;nT1++;}
       else{nT1=1;EneT1=Edep;timeT1=age;}
       iT1=i;
      }
    if (mreg == TrigReg[1])//T2
      {
       if(i==iT2+1){EneT2+=Edep;nT2++;}
       else{nT2=1;EneT2=Edep;timeT2=age;}
       iT2=i;
      }

    if (mreg == TrigReg[2])//T3
      {
       if(i==iT3+1){EneT3+=Edep;nT3++;}
       else{nT3=1;EneT3=Edep;timeT3=age;}
       iT3=i;
      }
    if (mreg == TrigReg[3])//T4
      {
       if(i==iT4+1){EneT4+=Edep;nT4++;}
       else{nT4=1;EneT4=Edep;timeT4=age;}
       iT4=i;
      }
    if (mreg == GReg[0])//guard
      {
       if(i==ig+1){Eneg+=Edep;ng++;}
       else{ng=1;Eneg=Edep;timeg=age;}
       ig=i;
      }
    
    for (int ii=0;ii<7;ii++)
      {
       if(mreg == TckReg[ii])
        {
	 if(i==iL+1) //if the segment of track is not the first of the hit 
	  {
	       h->set_DeltaE(h->get_DeltaE()+Edep);// Add the energy deposited in the segment
	       h->set_cx(h->get_cx()+cx);
           h->set_cy(h->get_cy()+cy); 
           h->set_cz(h->get_cz()+cz);
	       h->set_xout(x);
           h->set_yout(y);   
           h->set_zout(z);
		   h->set_x((h->get_xin() + x)/2.);
		   h->set_y((h->get_yin() + y)/2.);
		   h->set_z((h->get_zin() + z)/2.);
	       nL++;//Increment the number of segments that compose the hit
	  }//if
         else //if this segment of track is the first part of the hit
	  {
	   //Fill the data at the entrance of the particle in the tracker layer
	   h->set_mregMC(mreg); 
	   //set layer number (for MC and data)
           h->set_L((int)mreg%11);
           h->set_mtrackMC(mtrack);
           h->set_typeMC((int)type);
           h->set_eMC(pMC); 				//total energy of particle crossing region
           h->set_flag(flag);
           h->set_k(nh);
           h->set_age(age);
	   h->set_cx(cx);
           h->set_cy(cy); 
           h->set_cz(cz);
	   h->set_xin(x);
           h->set_yin(y);   
           h->set_zin(z);
	   h->set_xout(x);
           h->set_yout(y);   
           h->set_zout(z); 
	   h->set_DeltaE(Edep);
           nL=1;
          } //else
         iL=i; 
         ii=7;
	}  //if
      }//ii
    prevreg=mreg;
   }//i

  //Fill the last event
  if(nFoam>0)
     {
     e->add_EneIsofoam(EneFoam); e->add_timeIsofoam(timeFoam);
     }  	 
   
 if(nShell>0)
     {
     e->add_EneShell(EneShell); e->add_timeShell(timeShell);
       }
     
 if(nT1>0){
	 e->add_EneT1(EneT1);e->add_timeT1(timeT1);
	 if (EneT1 > TrigThresh[0]) e->set_T1(true);
 	}
 if(nT2>0){
         e->add_EneT2(EneT2);
	     e->add_timeT2(timeT2);
         if (EneT2 > TrigThresh[1]) e->set_T2(true);
        }

 if(nT3>0){
	 e->add_EneT3(EneT3);e->add_timeT3(timeT3);
	 if (EneT3 > TrigThresh[2]) e->set_T3(true);
 	}
 if(nT4>0){
	 e->add_EneT4(EneT4);e->add_timeT4(timeT4);
	 if(EneT4 > TrigThresh[3]) e->set_T4(true);
    }
 if(ng>0){
	 e->add_Eneg(Eneg);e->add_timeg(timeg);
	 if (Eneg > TrigThresh[4]) e->set_guard(true);
 	}
 if(nL>0){h->set_cx(h->get_cx()/nL);h->set_cy(h->get_cy()/nL);h->set_cz(h->get_cz()/nL);}//add the last hit
 //Check the layer for internal trigger
 for (int ij=0;ij<7;ij++)
   {
    if (h->get_mregMC()== TckReg[ij]) {Titmp[ij]=1;ij=7;}
   }
 int tt=0;
 for (int ij=0;ij<7;ij++) tt+=Titmp[ij]*(int)TMath::Power(2,ij);
 e->set_Ti(tt);	  
 tree->Fill(); 
 delete e;
 
 //Write tree in output file
 tree->Write();
 //Close files
 fileout->Close();
 file->Close();
 return 1;
}


int MakeRawEventMCDisc(int typeT,int Ene,int seed,int cycle,string Inppath,string Inppath2,string Outpath,string startfile,string endfile)
{
 cout << "Calling MakeRawEventMCDisc" << endl;   
 float OffsetLL=0;
 float OffsetRL=0;
    
//Load region numbers used in the MC geometry 
 int*TckReg=new int[7];
 int*TrigReg=new int[4];
 int*GReg=new int[1];
 int*ShellReg=new int[2];
 float*TckZPos=new float[7];
 float*TrigThresh=new float[4];
 float*GuardThresh=new float[1];
 for(int i=0;i<7;i++)TckReg[i]=0;
 for(int i=0;i<4;i++)TrigReg[i]=0;
 for(int i=0;i<2;i++)ShellReg[i]=0;
 for(int i=0;i<1;i++)GReg[i]=0;
 for(int i=0;i<7;i++)TckZPos[i]=0;
 for(int i=0;i<4;i++)TrigThresh[i]=0;
 for(int i=0;i<1;i++)GuardThresh[i]=0;
	
 string MCparamfile="../src/ALSim/MCparameters.dat"; 
 
 LoadMCparameters(MCparamfile,TckReg,TrigReg,GReg,TckZPos,TrigThresh,GuardThresh,ShellReg);
 cout << "Isofoam region: " << ShellReg[0] <<endl;
 cout << "Shell region: " << ShellReg[1] <<endl;
      
 TFile *file;
 TFile *fileout;
 //Input file 
 if(typeT==3 || typeT==4)
  {
   file=new TFile(Form("%s/%d/%s/%s_%d_%dMeV%d%03d%s.root",Inppath.c_str(),typeT,Inppath2.c_str(),startfile.c_str(),typeT,Ene,seed,cycle,endfile.c_str()),"READ");
   cout << "Input file is open" <<endl;
   //Output file 
   fileout=new TFile(Form("%s/%d/RawEvent_%s_%d_%dMeV%d%03d%s.root",Outpath.c_str(),typeT,startfile.c_str(),typeT,Ene,seed,cycle,endfile.c_str()),"RECREATE");
   cout << "Output file " << Form("%s/%d/RawEvent_%s_%d_%dMeV%d%03d%s.root",Outpath.c_str(),typeT,startfile.c_str(),typeT,Ene,seed,cycle,endfile.c_str()) << " is created" <<endl;
  }
 else
  {
   file=new TFile(Form("%s/%d/%s/%s_%d_%dGeV%03d%s.root",Inppath.c_str(),typeT,Inppath2.c_str(),startfile.c_str(),typeT,Ene,cycle,endfile.c_str()),"READ");
   cout << "Input file is open" << endl;
   fileout=new TFile(Form("%s/%d/RawEvent_%s_%d_%dGeV%03d%s.root",Outpath.c_str(),typeT,startfile.c_str(),typeT,Ene,cycle,endfile.c_str()),"RECREATE"); 
   cout << "Output file is created " << endl;
   } 
   
 //Get ntuple from the input file
 TNtuple*ntuple=(TNtuple*)file->Get("Track");
 cout << "Got the ntuple from the input file" <<endl;

 //Define variables to read ntuple
 float ncase=0; 
 float mreg=0;
 float mtrack=0;
 float type=0;
 float EkMC=0;
 float pMC=0; 
 float x=0;
 float y=0;
 float z=0;
 float age=0; 
 float cx=0;
 float cy=0;
 float cz=0; 
 float Edep=0; 
 float flag=0; 

 //Set addresses to access ntuple data
 
 ntuple->SetBranchAddress("ncase",&ncase); 
 ntuple->SetBranchAddress("mreg",&mreg);
 ntuple->SetBranchAddress("mtrack",&mtrack);
 ntuple->SetBranchAddress("type",&type);
 ntuple->SetBranchAddress("p",&pMC); 
 ntuple->SetBranchAddress("x",&x);
 ntuple->SetBranchAddress("y",&y);
 ntuple->SetBranchAddress("z",&z);
 ntuple->SetBranchAddress("age",&age); 
 ntuple->SetBranchAddress("cx",&cx);
 ntuple->SetBranchAddress("cy",&cy);
 ntuple->SetBranchAddress("cz",&cz); 
 ntuple->SetBranchAddress("flag",&flag); 
 ntuple->SetBranchAddress("Edep",&Edep); 

 
 cout << "Set Branch Address done" <<endl;

 
 // Get number of entries in ntuple (number of lines in the output fluka file) 
 int nentries=ntuple->GetEntries();
 cout << "Number  of entries: " << nentries << endl;
 int ievt=-1; // Identify the event
 // Create a TTree
 TTree *tree = new TTree("MC","Raw event MC");
 ALEvent *e = new ALEvent();
 // Create a branch with event
 tree->Branch("event",&e);  

 int j=0;
 int nh=0;
 double EneFoam=0;
 double EneShell=0;
 double EneT1=0;
 double EneT2=0; 
 double EneT3=0; 
 double EneT4=0; 
 double Eneg=0; 
 double timeFoam=0;
 double timeShell=0;
 double timeT1=0;
 double timeT2=0; 
 double timeT3=0; 
 double timeT4=0;
 int prevreg=0; 
 double timeg=0; 
 double nFoam=0;				//number of segments of track in insulation foam
 double nShell=0;				//number of segments of track in shell
 double nT1=0; 						//number of segments of track in T1
 double nT2=0;
 double nT3=0; 						
 double nT4=0; 
 double ng=0; 
 int iFoam=0;
 int iShell=0;
 int iT1=0; 
 int iT2=0;
 int iT3=0; 
 int iT4=0; 
 int ig=0;
 
 float oldx=0;
 float oldy=0;
 float oldz=0;
 
//Temporary variable for internal triggers 
 int*Titmp=new int[7];
 for (int i=0;i<7;i++)Titmp[i]=0;
 
 int nL=0;
 int iL=0;
 ALTckhit* h=new ALTckhit();

 //Vector of segments to discretize
 //Clear every time there is a new layer crossed
 vector<float> X;
 vector<float> Y;
 vector<float> Z;
 vector<float> CZ;
 vector<float> T;
 vector<float> EDEP;
 float CoordDisc=-999;
 int fstrip=-1;
 int fstripID=-1;
 int nstrip=0;
 int chip=-1;
 //Loop over the entries and create the events
 // This is the main loop
 for (int i=0;i<nentries;i++)
   {
    ntuple->GetEntry(i); //Load the entry i in the variables  
    if(i%100000==0) cout << "entry: " << i <<endl;
//    cout << "i = " << i << ", ncase " << ncase << ", mreg " << mreg << " type = " << type << " p = " << pMC*1000 <<" MeV,  x=" << x << ", y=" << y << ", z = " << z << endl; 
    e->add_posX(x);
    e->add_posY(y);
    e->add_posZ(z);
	e->add_posType(type);
	e->add_posAge(age);
	e->add_posP(pMC);
   
    ////////////////////////////////// 
    // Triggers
    ////////////////////////////////// 
    if(prevreg==ShellReg[0] && prevreg!=mreg)//Foam insulation
     {
      if(nFoam>0)
       {
        e->add_EneIsofoam(EneFoam); e->add_timeIsofoam(timeFoam);
		EneFoam=0;timeFoam=0;nFoam=0;iFoam=0;
       }
     }	 
    if(prevreg==ShellReg[1] && prevreg!=mreg)//Aluminium Shell
     {
      if(nShell>0)
       {
        e->add_EneShell(EneShell); e->add_timeShell(timeShell);
		EneShell=0;timeShell=0;nShell=0;iShell=0;
       }
     }	
    if(prevreg==TrigReg[0] && (prevreg!=mreg || abs(oldx-x)>0.1 || abs(oldy-y)>0.1 || abs(oldz-z)>0.1))//T1
     {
      if(nT1>0)
       {
        e->add_EneT1(EneT1); e->add_timeT1(timeT1);
        if (EneT1 > TrigThresh[0]) e->set_T1(true);
        EneT1=0;timeT1=0;nT1=0;iT1=0;
       }
     }
   if(prevreg==TrigReg[1] && (prevreg!=mreg || abs(oldx-x)>0.1 || abs(oldy-y)>0.1 || abs(oldz-z)>0.1))//T2
     {
      if(nT2>0)
       {
        e->add_EneT2(EneT2);
        e->add_timeT2(timeT2);
        if (EneT2 > TrigThresh[1]) e->set_T2(true);
        EneT2=0;timeT2=0;nT2=0;iT2=0;
       }
     }
    if(prevreg==TrigReg[2] && (prevreg!=mreg || abs(oldx-x)>0.1 || abs(oldy-y)>0.1 || abs(oldz-z)>0.1))//T3
     {
      if(nT3>0)
       {
        e->add_EneT3(EneT3); e->add_timeT3(timeT3);
        if(EneT3 > TrigThresh[2]) e->set_T3(true);
        EneT3=0;timeT3=0;nT3=0;iT3=0;
       }
     }
    if(prevreg==TrigReg[3] && (prevreg!=mreg || abs(oldx-x)>0.1 || abs(oldy-y)>0.1 || abs(oldz-z)>0.1))//T4
     {
      if(nT4>0)
       {
        e->add_EneT4(EneT4); e->add_timeT4(timeT4);
        if (EneT4 > TrigThresh[3]) e->set_T4(true);
        EneT4=0;timeT4=0;nT4=0;iT4=0;
       }
     }
    if(prevreg==GReg[0] && (prevreg!=mreg || abs(oldx-x)>0.1 || abs(oldy-y)>0.1 || abs(oldz-z)>0.1))//guard

    {
      if(ng>0)
       {
        e->add_Eneg(Eneg); e->add_timeg(timeg);
        if (Eneg > TrigThresh[4]) e->set_guard(true);
        Eneg=0;timeg=0;ng=0;ig=0;
       }
     }

    ////////////////////////////////// 
    // Tracker
    ////////////////////////////////// 
     
    if(prevreg>=TckReg[0] && prevreg<=TckReg[6]&& (prevreg!=mreg || abs(oldx-x)>0.1 || abs(oldy-y)>0.1 || abs(oldz-z)>0.1))
     {       
      
      //APPLY DISCRETISATION
      CoordDisc=-999.;
      fstrip=-1;
      nstrip=0;
      chip=-1;
//      cout << "size of vector x = " << (int)X.size() << endl;      
      if((int)X.size()!=0) 
       {
//	cout << "applying discretization" << endl;
        CoordDisc=Discretize(h->get_L(),X,Y,Z,CZ,T,EDEP,&chip,&fstrip,&fstripID,&nstrip,OffsetLL,OffsetRL,true);
//        cout << "CoordDisc=" <<  CoordDisc <<endl; 
        if(CoordDisc!=-999.)
         {   
          //Update parameters of the hit after discretisation
          if(h->get_L()==0||h->get_L()==4||h->get_L()==6) h->set_x(CoordDisc); //Non bending plane
          if(h->get_L()==1||h->get_L()==2||h->get_L()==3||h->get_L()==5) h->set_y(CoordDisc);//Bending plane
          h->set_fstripID(fstripID);//from 0 to 767
          h->set_fstrip(fstrip);//from 0 to 63
          h->set_nstrips(nstrip);//number of strips in cluster
          h->set_chip(chip);//chip 0 to 11 of the first strip (strip lowest index)
          if(nL>0){h->set_cx(h->get_cx()/nL);h->set_cy(h->get_cy()/nL);h->set_cz(h->get_cz()/nL);}
          //End of the hit.
          //Add the hit to the event
	  //index k of hit starts at 0
          h->set_k(nh);
	  nh++;
          if(h->get_x()==-999) cout << "X coord is -999. Event " << i <<endl;
          if(h->get_y()==-999) cout << "Y coord is -999. Event " << i <<endl;
          if(h->get_z()==-999) cout << "Z coord is -999. Event " << i <<endl;

          e->add_hit(h);
          Titmp[(int)h->get_L()]=1;
         } 
        X.clear();
        Y.clear();
        Z.clear();
        CZ.clear();
        T.clear();
        EDEP.clear();
       //Check the layer for internal trigger       
       }
 

      //Reset the hit structure
      h=new ALTckhit();
      nL=0;
      iL=0;
     }  
    if(ncase!=ievt)//Start a new event
     {
      j++;
      nh=0;

      if(i==0)//if first event
       {
        e->add_posX(x);
        e->add_posY(y);
        e->add_posZ(z);
		e->add_posType(type);
		e->add_posAge(age);
		e->add_posP(pMC);
        ievt=ncase; 
        delete e;
        e=new ALEvent();
        e->set_eventnumber(j);
        e->set_ncase(ncase);
        //First line  with a given ncase contains the source particle info
        e->set_typeMC(type);
        e->set_pMC(pMC);
        //Calculate kinetic energy from momentum at point of injection
        float mass=0.000511;								//electron mass in GeV
        if(type==11 || type ==10)	mass=0.10566;			//muon mass in GeV
		if(type==1)					mass=0.93827;			//proton mass in GeV
		if(type==-6)				mass=3.72739;			//alpha mass in GeV 
        EkMC = TMath::Sqrt(pMC*pMC+mass*mass);
        e->set_EkMC(EkMC);
        e->set_X0MC(x);
        e->set_Y0MC(y);
        e->set_Z0MC(z);
        e->set_CX0MC(cx);
        e->set_CY0MC(cy);
        e->set_CZ0MC(cz);
        //First attempt at filling T2 
       // e->set_T2(false);
	for(int ij=0;ij<7;ij++)Titmp[ij]=0;
       }//if
      else
       {
        //Fill internal trigger information of the previous event
        int t=0;
        for (int ij=0;ij<7;ij++) t+=Titmp[ij]*(int)TMath::Power(2,ij);
        e->set_Ti(t);  
        //Fill the previous event in the output file  
        tree->Fill(); 
        delete e;
        e=new ALEvent();
        h=new ALTckhit();
        nh=0;
        ievt=ncase; 
        e->set_eventnumber(j);
        e->set_ncase(ncase);
        //First line with a given ncase contains the source particle info
        e->set_typeMC(type);
        e->set_pMC(pMC);
        //Calculate kinetic energy from momentum at point of injection
        float mass=0.000511;								//electron mass in GeV
        if(type==11 || type ==10)	mass=0.10566;			//muon mass in GeV
		if(type==1)					mass=0.93827;			//proton mass in GeV
		if(type==-6)				mass=3.72739;			//alpha mass in GeV 
        float EkMC = TMath::Sqrt(pMC*pMC+mass*mass);
        e->set_EkMC(EkMC);
        e->set_X0MC(x);
        e->set_Y0MC(y);
        e->set_Z0MC(z);
//cout << " x0 = " << x << " y0 = " << y << "  z0 = " << z << endl;
        e->set_CX0MC(cx);
        e->set_CY0MC(cy);
        e->set_CZ0MC(cz);
        //first attempt to fill T2
       // e->set_T2(false);
        //Reset internal trigger
        for (int ij=0;ij<7;ij++)Titmp[ij]=0;
       }     //else
     }  //if
    
    ////////////////////////////////// 
    //Check to see if the particle has crossed one of the scintallator 
    //(r19=air, r1=T1, r6=T3, r7=guard, r11=Tracker1 ,..., r17=Tracker7, r18=T4)
    ////////////////////////////////// 
	//Shell & Insulation
    if (mreg == ShellReg[0]) //Insulating foam
      {
       if(i==iFoam+1){EneFoam+=Edep;nFoam++;}
       else{nFoam=1;EneFoam=Edep;timeFoam=age;}
       iFoam=i;
      }
    if (mreg == ShellReg[1])//Aluminium Shell
      {
       if(i==iShell+1){EneShell+=Edep;nShell++;}
       else{nShell=1;EneShell=Edep;timeShell=age;}
       iShell=i;
      }	     
    //Triggers
    if (mreg == TrigReg[0])//T1
      {
       if(i==iT1+1){EneT1+=Edep;nT1++;}
       else{nT1=1;EneT1=Edep;timeT1=age;}
       iT1=i;
      }
   if (mreg == TrigReg[1])//T2
      {
       if(i==iT2+1){EneT2+=Edep;nT2++;}
       else{nT2=1;EneT2=Edep;timeT2=age;}
       iT2=i;
      }

    if (mreg == TrigReg[2])//T3
      {
       if(i==iT3+1){EneT3+=Edep;nT3++;}
       else{nT3=1;EneT3=Edep;timeT3=age;}
       iT3=i;
      }
    if (mreg == TrigReg[3])//T4
      {
       if(i==iT4+1){EneT4+=Edep;nT4++;}
       else{nT4=1;EneT4=Edep;timeT4=age;}
       iT4=i;
      }
    if (mreg == GReg[0])//guard
      {
       if(i==ig+1){Eneg+=Edep;ng++;}
       else{ng=1;Eneg=Edep;timeg=age;}
       ig=i;
      }
    
    for (int ii=0;ii<7;ii++)
      {
       if(mreg == TckReg[ii])
        {
         if(i==iL+1 && abs(oldx-x)<=0.1 && abs(oldy-y)<=0.1 && abs(oldz-z)<=0.1 ) //if the segment of track is not the first of the hit 
          {
           h->set_DeltaE(h->get_DeltaE()+Edep);// Add the energy deposited in the segment
           h->set_cx(h->get_cx()+cx);
           h->set_cy(h->get_cy()+cy); 
           h->set_cz(h->get_cz()+cz);
           h->set_xout(x);
           h->set_yout(y);   
           h->set_zout(z);
           h->set_x((h->get_xin() + x)/2.);
           h->set_y((h->get_yin() + y)/2.);
           h->set_z((h->get_zin() + z)/2.);	
           X.push_back(x);
           Y.push_back(y);
           Z.push_back(z);
          // cout << "  x = " << x << "  y = " << y << "  z = " << z << endl;
           T.push_back(type);
           CZ.push_back(cz);
           EDEP.push_back(Edep);
           nL++;//Increment the number of segments that compose the hit
          }//if
         else //if this segment of track is the first part of the hit
          {
           //Fill the data at the entrance of the particle in the tracker layer
           h->set_mregMC(mreg); 
           //set layer number (for MC and data)
           h->set_L((int)mreg%11);
           h->set_mtrackMC(mtrack);
           h->set_typeMC((int)type);
           h->set_eMC(pMC); 				//total energy of particle crossing region
           h->set_flag(flag);
           h->set_k(nh);
           h->set_age(age);
           h->set_cx(cx);
           h->set_cy(cy); 
           h->set_cz(cz);
           h->set_xin(x);
           h->set_yin(y);   
           h->set_zin(z);
           h->set_xout(x);
           h->set_yout(y);   
           h->set_zout(z);
	  // cout << "Layer L = " << (int)mreg%11 << "  x = " << x << "  y = " << y << "  z = " << z << endl;
           h->set_x(x);
           h->set_y(y);
           h->set_z(z);
           X.push_back(x);
           Y.push_back(y);
           Z.push_back(z);
           CZ.push_back(cz);
           T.push_back(type);
           EDEP.push_back(Edep);
           h->set_DeltaE(Edep);
           nL=1;
          } //else
         iL=i; 
         ii=7;
        }  //if
      }//ii
    prevreg=mreg;
    oldx=x;
    oldy=y;
    oldz=z;
   }//i

  //Fill the last event
  if(nFoam>0) {
     e->add_EneIsofoam(EneFoam); e->add_timeIsofoam(timeFoam); 
  }  	 
 if(nShell>0) {
     e->add_EneShell(EneShell); e->add_timeShell(timeShell);
       }
 if(nT1>0){
	 e->add_EneT1(EneT1);e->add_timeT1(timeT1);
	 if (EneT1 > TrigThresh[0]) e->set_T1(true);
 	}
 if(nT2>0){
         e->add_EneT2(EneT2);
         //e->add_timeT2(timeT2);
         if (EneT2 > TrigThresh[1]) e->set_T2(true);
        }

 if(nT3>0){
	 e->add_EneT3(EneT3);e->add_timeT3(timeT3);
	 if (EneT3 > TrigThresh[2]) e->set_T3(true);
 	}
 if(nT4>0){
	 e->add_EneT4(EneT4);e->add_timeT4(timeT4);
	 if(EneT4 > TrigThresh[3]) e->set_T4(true);
    }
 if(ng>0){
	 e->add_Eneg(Eneg);e->add_timeg(timeg);
	 if (Eneg > TrigThresh[4]) e->set_guard(true);
 	}
 
 
 //APPLY DISCRETISATION
 CoordDisc=-999.;
 fstrip=-1;
 fstripID=-1;
 nstrip=0;
 chip=-1;
 
 if((int)X.size()!=0) 
  {
   CoordDisc=Discretize(h->get_L(),X,Y,Z,CZ,T,EDEP,&chip,&fstrip,&fstripID,&nstrip,OffsetLL,OffsetRL,true);
  cout << "CoordDisc2=" <<  CoordDisc <<endl; 
   
   if(CoordDisc!=-999.)
    {   
     //Update parameters of the hit after discretisation
     if(h->get_L()==0||h->get_L()==4||h->get_L()==6) h->set_x(CoordDisc); //Non bending plane
     if(h->get_L()==1||h->get_L()==2||h->get_L()==3||h->get_L()==5) h->set_y(CoordDisc);//Bending plane
     h->set_fstripID(fstripID);//from 0 to 767
     h->set_fstrip(fstrip);//from 0 to 62
     h->set_nstrips(nstrip);//number of strips in cluster
     h->set_chip(chip);//chip 0 to 11
     if(nL>0){h->set_cx(h->get_cx()/nL);h->set_cy(h->get_cy()/nL);h->set_cz(h->get_cz()/nL);}//add the last hit 
     //End of the hit.
     //Add the hit to the event
     h->set_k(nh);
     nh++;
     e->add_hit(h);
    //Check the layer for internal trigger
     Titmp[(int)h->get_L()]=1;
    } 
   X.clear();
   Y.clear();
   Z.clear();
   CZ.clear();
   T.clear();
   EDEP.clear();
      
  }

 int tt=0;
 for (int ij=0;ij<7;ij++) tt+=Titmp[ij]*(int)TMath::Power(2,ij);
 e->set_Ti(tt);	  
 tree->Fill(); 
 delete e; 

 //Write tree in output file 
 tree->Write();
 //Close files
 fileout->Close();
 file->Close();
 return 1;
}

/*

float Discretize(int L,vector<float> x, vector<float> y,vector<float> z,vector<float> cz,vector<float>type,vector<float> Edep,int*chip,int* fstrip,int* fstripID,int*nstrip,float offsetLL, float offsetRL,bool MCflag)
{
 //cout << "In Discretize" <<endl;
 int Nn=(int )x.size(); 
 // cout << "Number of segments: "  <<Nn << " in Layer "<< L << endl;

 float* X=new float[Nn];
 float* Y=new float[Nn];
 float* Z=new float[Nn];
 float* CZ=new float[Nn];
 float* T=new float[Nn];
 float* E=new float[Nn];
 float* xx=new float[Nn];
 float* Ss=new float[Nn];
 //Output
 float Xout=-999.;
 //Strip pitch in cm
 float  strippitch=0.0228; 

 for(int i=0;i<Nn;i++)
  { 
   if(L==0||L==4||L==6)
    {
     X[i]=x.at(i);
     Y[i]=y.at(i); 
    }
   else if(L==1||L==2||L==3||L==5)
    {
     X[i]=y.at(i);
     Y[i]=x.at(i); //Inverse coordinated for bending plane	
    }
   E[i]=Edep.at(i);  
   Z[i]=z.at(i); 
   CZ[i]=cz.at(i); 
   T[i]=type.at(i); 
   Ss[i]=-1; //Set to defaults
   //cout << "X"<< i << "= " << X[i] << ", " ;
   //cout << "Y"<< i << "= " << Y[i] << ", "  ;
   //cout << "Z"<< i << "= " << Z[i] <<endl;
  }//i  

//ASSUMPTION: the four frames are perfectly aligned without space between each of them

//Get stripID for each segment of the track in the layer
for(int i=0;i<Nn;i++)
 {  
  Ss[i]=CoordtoStrip(X[i],Y[i],offsetLL,offsetRL,MCflag);
  //cout << "StripID="<< Ss[i] <<endl ;
 }//i

//Check if there is any strip touched
bool nostrip=true;
for(int i=0;i<Nn;i++)
 {  
  if(Ss[i]>-1) {nostrip=false;i=Nn;}
 } //i
 
if(nostrip)
 { 
 // cout << "No strip hitted" <<endl;
  return Xout;
 }

//Get first and last a segment in one strip
int kfirst=-1;
int klast=-1;

for(int i=0;i<Nn;i++)
  {
   cout << "i = " << i << ", y = " << Y[i] << " , Ss[i] = " << Ss[i] << endl;
   if(kfirst==-1 && Ss[i]>-1) {kfirst=i;klast=i;}    
   if(kfirst>-1 && Ss[i]>-1) klast=i;    
  }
  
if(kfirst<0 || klast<0)
 { 
  return Xout; 
 }

 Xout= (StriptoCoord(Ss[kfirst],offsetLL,offsetRL,MCflag)+StriptoCoord(Ss[klast],offsetLL,offsetRL,MCflag))/2.;

 if(Ss[klast]>=Ss[kfirst]) *fstripID=Ss[kfirst];
 else *fstripID=Ss[klast];

 *chip=(int)Ss[kfirst]%12;
 *fstrip=(int)Ss[kfirst]%64;
 
 *nstrip=(int) abs(Ss[klast]-Ss[kfirst])+1;
 
// cout << "Xout="<< Xout <<endl ;

 
 if(*nstrip>100)
  {
   cout << "Number of segments: "  <<Nn << " in Layer "<< L << endl;
   for(int i=0;i<Nn;i++)
     { 
      cout << "X"<< i << "= " << X[i] << ", " ;
      cout << "Y"<< i << "= " << Y[i] << ", "  ;
      cout << "Z"<< i << "= " << Z[i] << ", ";
      cout << "CZ"<< i << "= " << CZ[i] << ", ";
      cout << "Type"<< i << "= " << T[i] << endl;
      //cout << "Edep"<< i << "= " << Edep[i] <<endl;
     }//i
  }
 return Xout;
}


float StriptoCoord(int strip,float OffsetLL,float OffsetRL,bool MCflag)
{
 //Determine the coordinates from the strips
 //Equations from Sarah's email of September 4 2017.
 
 //Center of X,Y from alignment PIN in cm from Robert (January 24th 2018)
 float Xo=9.74745;   
 float Yo=4.48005; 
 
 float coord=-999;
 
 //Strip pitch in cm
 float  strippitch=0.0228;

 //for MC: 
 float N=384;                                //Number of strip in one module
 float Offset1=0.1088;                       //In cm. Offset from the center of the first/last strip to the edge of the module.
 float Offset2=0.0964;                       //In cm. Offset from the start/end of the strips to the edge of the module.

 float offsetLL=OffsetLL;//Left ladder offset
 float offsetRL=OffsetRL;//Right ladder offset   

 if(MCflag)//Assume perfect alignmentof the four pads
  {
   offsetLL=Xo-Offset1-(N-1)*strippitch; 
   offsetRL=Xo+Offset1; 
  }
 
 //First 6 chips: 0 to 5; strip number 0 to 383
 if(strip>=0 &&strip<N)
  {
   coord=offsetLL-Xo+strip*strippitch; 
  }
 //Last 6 chips: 6 to 11; strip number 384 to 767
 if(strip>=N)
  {
   coord=offsetRL-Xo+(strip-N)*strippitch;
  }
 
 return coord;
}


int CoordtoStrip(float Coord,float SecCoord,float OffsetLL,float OffsetRL,bool MCflag)
{
 //Geometry parameters from schematics
 //Used for discretisation
 float N=384;                                //Number of strip in one module
 float Offset1=0.1088;                       //In cm. Offset from the center of the first/last strip to the edge of the module.
 float Offset2=0.0964;                       //In cm. Offset from the start/end of the strips to the edge of the module.
 float sw=0.0056;                            //In cm. Width of a strip
 float sl=8.7572;                            //In cm. Length of a strip in one module = Size of the frame (8.95cm) - 2*Offset2  ??????????????????????????????
 float sh=0.04;                              //In cm. Height of a strip
 float percsh=0.2;                           //fraction of sh crossed by track to trigger a strip: To be fine tuned later
 float DeltaMax=0.0114;                      //In cm. Maximum allowed distance to the center of a strip
 float strippitch=0.0228;                    //In cm. Strip pitch 
 int strip=-1;
 
 //Center of X,Y from alignment PIN in cm from Robert (January 24th 2018)
 float Xo=9.74745;   
 float Yo=4.48005; 
 float offsetLL=OffsetLL;//Left ladder offset
 float offsetRL=OffsetRL;//Right ladder offset   
 if(MCflag)//Assume perfect alignmentof the four pads
  {
   offsetLL=Xo-Offset1-(N-1)*strippitch; 
   offsetRL=Xo+Offset1; 
  }
 float strip0=offsetLL;       //In cm. Position of the center of strip 0 

 //First: Remove hits outside the frame along Y
 if(SecCoord>sl+Offset2+DeltaMax) return strip;
 if(SecCoord<-sl-Offset2-DeltaMax) return strip;

 //Second: Remove hits outside the frame along X
 if(Coord>offsetRL-Xo+(N-1)*strippitch+DeltaMax) return strip;
 if(Coord<offsetLL-Xo+strippitch-DeltaMax) return strip;

 //Third: Remove Inner Cross
 if(Coord<offsetRL-Xo-DeltaMax && Coord>offsetLL-Xo+(N-1)*strippitch+DeltaMax) return strip;
 if(abs(SecCoord)<Offset2-DeltaMax) return strip;
 
 //Fourth: Determine strip ID from 0 to 767 
 float relX=0;
 int IrelX=0; 
 float RrelX=0;

 if(Coord<0)
  {
   relX=Coord+Xo-offsetLL;   
   IrelX=(int)(relX/strippitch);
   RrelX=relX-(float)IrelX*strippitch;  
   if(RrelX>DeltaMax)IrelX++;
   else if(RrelX<-DeltaMax)IrelX--;
   if(IrelX>=0 && IrelX<N) strip=IrelX; 
  }  
  
 if(Coord>0) 
  {
   relX=Coord+Xo-offsetRL;   
   IrelX=(int)(relX/strippitch);
   RrelX=relX-(float)IrelX*strippitch;  
   if(RrelX>DeltaMax)IrelX++;
   else if(RrelX<-DeltaMax)IrelX--;
   if(IrelX>=0 && IrelX<N) strip=IrelX+N; 
  }    
 return strip;
    
}
*/



