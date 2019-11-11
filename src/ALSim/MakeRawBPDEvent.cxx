//////////////////////////////////////////////////////////////////////////////////////////
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, September 8, 2017
//////////////////////////////////////////////////////////////////////////////////////////

#include "MakeRawBPDEvent.h"

int MakeRawBPDEvent(string filename)
{
 //Get input filename and ASIC Length
 vector<string> input;
 //Split the line
 input=split(&filename,' ');
 int ASILength =stoi(input.at(1));
 int EVTLength =stoi(input.at(2));
 //Open input file
 ifstream filestr;
 filestr.open(input.at(0));
 if(!filestr.is_open())
  {
   cout << "The file " <<input.at(0) << " is not open ... "<< endl ;
   return -1;
  }




 //Load configuration parameter
/* int*TckReg=new int[7];
 int*TrigReg=new int[4];
 int*GReg=new int[1];for(int i=0;i<1;i++)GReg[i]=0;
 for(int i=0;i<7;i++)TckReg[i]=0;
 for(int i=0;i<4;i++)TrigReg[i]=0;
 for(int i=0;i<1;i++)GReg[i]=0;
 string MCparamfile="../src/ALSim/MCparameters.dat";

 LoadMCparameters(MCparamfile,TckReg,TrigReg,GReg);
*/

 //Output root file
 TFile*fileout=new TFile(Form("%s.root",input.at(0).c_str()),"RECREATE");
 cout << "Output file is created" <<endl;

 // Create a TTree
 TTree *tree = new TTree("BPD","BPD Raw event");
 ALEvent *e = new ALEvent();

 // Create a branch with event
 tree->Branch("event",&e);

//Temporary variable for internal triggers
 int*Titmp=new int[7];
 for (int i=0;i<7;i++)Titmp[i]=0;
 int*Tictmp=new int[7];
 for (int i=0;i<7;i++)Tictmp[i]=0;

 int nL=0;
 int iL=0;
 vector<ALTckhit*> Hh;


///////////////////////////////////////////////
//Define variables needed during the reading of the file
///////////////////////////////////////////////

 //Line number
 int kLine=0;
 //Event number
 int kEvent=0;

 //vector of the last "PHA" lines
 vector <string> prevPHAline;
 //Number of "PHA" lines kept in memory
 int NPHAMem=5;
 //Last "EVT" lines
 string prevEVTline;

 //number of words per line
 //int NPHAwords=13; //PHA line for MIPFlit3.30
 //int NPHAwords=11; //PHA line for MIPFlit3.31
 int NPHAwords=12; //PHA line for MIPFlit3.72
 int NEVTwords=0;
 if(EVTLength==0)  NEVTwords=4; //EVT line from MIPFlit3.30 to 3.35
 if(EVTLength==1)  NEVTwords=16; //Was increased at the same time of ASIC Length
 if(EVTLength==2)  NEVTwords=18; //Was increased at the same time of ASIC Length from NL0155
 int NASIwords= 4; //ASI line


//Number of Layers read
 int kASI=0;

//Start to read the file
 for (string line; getline(filestr,line);)
   {
   // if(kEvent>3) return -1;
    //Here the full line is read and stored in the string object "line"
    kLine++;
    //if(kLine%1000==0)
    //      cout<< "Line " << kLine << " is read" <<endl;
///////////////////////////////////////////////
//Compare the first 3 char of the line to "PHA"
///////////////////////////////////////////////

    if (line.compare(0,3,"PHA",3)==0)
     {
      //cout << "In PHA line" <<endl;
      //Each PHA Line creates one event
      //If not the first PHA Line fill the previous event in the tree
      if(kEvent!=0)
       {
        int t=0;
        int tc=0;
        for (int ij=0;ij<7;ij++) t+=Titmp[ij]*(int)TMath::Power(2,ij);
        for (int ij=0;ij<7;ij++) tc+=Tictmp[ij]*(int)TMath::Power(2,ij);
        e->set_Ti(t);
        e->set_Tic(tc);
        tree->Fill();
        //reset event
        delete e;
       }

      //Define new event
      e=new ALEvent();
      //Increment the number of read events
      kEvent++;
      e->set_eventnumber(kEvent);

      //Define internal tracker trigger
      Titmp=new int[7];
      for(int ij=0;ij<7;ij++)Titmp[ij]=0;
      Tictmp=new int[7];
      for(int ij=0;ij<7;ij++)Tictmp[ij]=0;
      ///////////////////////////////////////////////
      //Extract the data from the last PHA line
      ///////////////////////////////////////////////
      vector<string> datalinePHA;

      //Split the line
      datalinePHA=split(&line,' ');

      if((int)datalinePHA.size()!=NPHAwords)
        {
         cout << "Number of words in the PHA line is wrong : " << (int)datalinePHA.size() << " (Should be " << NPHAwords  << "): "<< endl;
         cout << line <<endl;
         continue;
        }
      int year=0;
      int month=0;
      int day=0;
      int hour=0;
      int minute=0;
      int second=0;

      //Extract PHA date
      year=month=day=0;
      //cout << "extract PHA date" <<endl;
      //cout << dataline.at(1) <<endl;
      extractdate(&year,&month,&day,&datalinePHA.at(1));
      //cout << "extract PHA date is done" <<endl;
      //Extract PHA time
      hour=minute=second=0;
      //cout << "extract PHA time" <<endl;
      extracttime(&hour,&minute,&second,&datalinePHA.at(2));
      //cout << "extract PHA time is done" <<endl;
      //Fill the event structure with PHA time variables
      e->set_yPHA(year);
      e->set_mPHA(month);
      e->set_dPHA(day);
      e->set_hPHA(hour);
      e->set_miPHA(minute);
      e->set_sPHA(second);
      //cout << "PHA time was added to the event" <<endl;
      //Fill the trigger information from the PHA line
      e->add_EneT1(s2lf(&datalinePHA.at(3)));
      e->add_EneT2(s2lf(&datalinePHA.at(4)));
      e->add_EneT3(s2lf(&datalinePHA.at(5)));
      e->add_EneT4(s2lf(&datalinePHA.at(6)));
      e->add_Eneg(s2lf(&datalinePHA.at(7)));
      e->add_PHA6(s2lf(&datalinePHA.at(8)));
      e->set_GoPHA(s2lf(&datalinePHA.at(9)));
      e->set_tPHA(s2lf(&datalinePHA.at(10)));
     }//end read line starting with PHA

///////////////////////////////////////////////
//Compare the first 3 char of the line to "EVT"
///////////////////////////////////////////////
    if (line.compare(0,3,"EVT",3)==0)
     {
      //cout << "In EVT line" <<endl;
      //An EVT line is always associated to the previous PHA line
      //If several EVT lines follow one PHA line then the last EVT line is recorded
      ///////////////////////////////////////////////
      //Extract the data from the EVT line
      ///////////////////////////////////////////////
      vector<string> datalineEVT;

      //Split the line
      datalineEVT=split(&line,' ');
      if((int)datalineEVT.size()!=NEVTwords)
       {
        cout << "Number of words in the EVT line is wrong : " << (int)datalineEVT.size() << " (Should be " << NEVTwords  << "): "<< endl;
        cout << line <<endl;
       }
      int year=0;
      int month=0;
      int day=0;
      int hour=0;
      int minute=0;
      int second=0;
      //Extract EVT date
      //cout << "extract EVT date" <<endl;
      //cout << dataline.at(1) <<endl;
      extractdate(&year,&month,&day,&datalineEVT.at(1));
      //cout << "extract EVT date is done" <<endl;
      //Extract PHA time
      hour=minute=second=0;
      //cout << "extract EVT time" <<endl;
      extracttime(&hour,&minute,&second,&datalineEVT.at(2));
      //cout << "extract EVT time is done" <<endl;
      //Fill the event structure with EVT time variables
      e->set_yEVT(year);
      e->set_mEVT(month);
      e->set_dEVT(day);
      e->set_hEVT(hour);
      e->set_miEVT(minute);
      e->set_sEVT(second);
      //cout << "EVT time was added to the event" <<endl;
      //Fill data
      e->set_EVT(datalineEVT.at(3));
      if(EVTLength==1)
       {
        e->set_GoEVT(s2lf(&datalineEVT.at(11)));
        e->set_tEVT(s2lf(&datalineEVT.at(13)));
       }
      if(EVTLength==2)
       {
        e->set_nHitLEVT(s2i(&datalineEVT.at(5)));
        e->set_CCEVT(s2i(&datalineEVT.at(7)));
        e->set_PatternEVT(s2i(&datalineEVT.at(9)));
        e->set_Q1EVT(s2i(&datalineEVT.at(11)));
        e->set_GoEVT(s2lf(&datalineEVT.at(13)));
        e->set_tEVT(s2lf(&datalineEVT.at(15)));
        e->set_TrigEVT(s2lf(&datalineEVT.at(17)));
       }

      // cout << "EVT data was added to the event" <<endl;
     }//end read line starting with EVT

///////////////////////////////////////////////
//Compare the first 3 char of the line to "ASI"
///////////////////////////////////////////////
    if (line.compare(0,3,"ASI",3)==0)
     {
      //cout << "In ASI line" <<endl;
      //Define vector of string "dataline"
      //The vector will contains the words of the line
      vector<string> dataline;
      dataline=split(&line,' '); //split the line in a vector of words

      if((int)dataline.size()!=NASIwords)
       {
        cout << "Number of words in the ASI line is wrong : " << (int)dataline.size() << " (Should be " << NASIwords  << "): "<< endl;
        cout << line <<endl;
        return -1;
       }

      //Extract ASI date
      int yearASI=0;
      int monthASI=0;
      int dayASI=0;
      //cout << "extract ASI date" <<endl;
      //cout << dataline.at(1) <<endl;
      extractdate(&yearASI,&monthASI,&dayASI,&dataline.at(1));
      //cout << "extract ASI date is done" <<endl;
      //Extract ASI time
      int hourASI=0;
      int minuteASI=0;
      int secondASI=0;
      //cout << "extract ASI time" <<endl;
      extracttime(&hourASI,&minuteASI,&secondASI,&dataline.at(2));
      //cout << "extract ASI time is done" <<endl;

      //Decode ASIC code
      if(ASILength==0)DecodeASIShort(dataline.at(3),&Hh,Titmp);
      int n=0;
      if(ASILength==1)DecodeASILong(dataline.at(3),&Hh,Titmp,Tictmp,&n);

      for(int ij=0;ij<(int)Hh.size();ij++)
        {
         //e->set_L(Hh.at(ij)->getL(),dataline.at(3));
         e->set_flagL(Hh.at(ij)->get_L(),1);
         Hh.at(ij)->set_year(yearASI);
         Hh.at(ij)->set_m(monthASI);
         Hh.at(ij)->set_d(dayASI);
         Hh.at(ij)->set_hour(hourASI);
         Hh.at(ij)->set_mi(minuteASI);
         Hh.at(ij)->set_s(minuteASI);
         e->add_hit(Hh.at(ij));
        }
      //Clear hits
      Hh.clear();
     }//end read line starting with ASI
   }//end read file


 //Fill last event
 int t=0;
 int tc=0;
 for (int ij=0;ij<7;ij++) t+=Titmp[ij]*(int)TMath::Power(2,ij);
 for (int ij=0;ij<7;ij++) tc+=Tictmp[ij]*(int)TMath::Power(2,ij);
 e->set_Ti(t);
 e->set_Tic(tc);
 tree->Fill();
 //reset event
 delete e;


 //Write tree in output file
 tree->Write(0,TObject::kOverwrite);

 //Close the file
 filestr.close();
 cout << "The file "<<filename  << " is closed"<< endl;
 fileout->Close();
 cout << "The output file "<<Form("%s.root",input.at(0).c_str())  << " is closed"<< endl;

 return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////
//Test readout with internal trigger on
//////////////////////////////////////////////////////////////////////////////////////////
int MakeRawBPDEventIT(string filename)
{

  cout << "In MakeRawBPDEvent"  <<endl;


 //Get input filename and ASIC Length and EVT Length
 vector<string> input;
 //Split the line
 input=split(&filename,' ');
 int ASILength =stoi(input.at(1));
 int EVTLength =stoi(input.at(2));
 //Open input file
 ifstream filestr;
 filestr.open(input.at(0));
 if(!filestr.is_open())
  {
   cout << "The file " <<input.at(0) << " is not open ... "<< endl ;
   return -1;
  }

 //Output root file
 TFile*fileout=new TFile(Form("%s.root",input.at(0).c_str()),"RECREATE");
 cout << "Output file is created" <<endl;

 // Create a TTree
 TTree *tree = new TTree("BPD","BPD Raw event");
 ALEvent *e = new ALEvent();

 // Create a branch with event
 tree->Branch("event",&e);

//Temporary variable for internal triggers
 int*Titmp=new int[7];
 for (int i=0;i<7;i++)Titmp[i]=0;
 int*Tictmp=new int[7];
 for (int i=0;i<7;i++)Tictmp[i]=0;

 int nL=0;
 int iL=0;
 vector<ALTckhit*> Hh;


///////////////////////////////////////////////
//Define variables needed during the reading of the file
///////////////////////////////////////////////

 //Line number
 int kLine=0;
 //Event number
 int kEvent=0;

 //vector of the last "PHA" lines
 vector <string> prevPHAline;
 //Number of "PHA" lines kept in memory
 int NPHAMem=5;
 //Last "EVT" lines
 string prevEVTline=" ";
 double prevGo=-999;


 //LAST CT1 and CT3 line info
  //HOUSEKEEPING FROM COUNTERS 1 AND 3
 //Data FROM "CT1" LINE
 int yCT1=-1;//Year from CT1 line linked to the event (last read CT1 line)
 int mCT1=-1;//Month from CT1 line linked to the event (last read CT1 line)
 int dCT1=-1;//Day from CT1 line linked to the event (last read CT1 line)
 int hCT1=-1;//Hour from CT1 line linked to the event (last read CT1 line)
 int miCT1=-1;//Minute from CT1 line linked to the event (last read CT1 line)
 int sCT1=-1;//Second from CT1 line linked to the event (last read CT1 line)
 float TempCT1=-999; //Temperature measured on the board CT1

 int OnTimeCT1=-1;//1/second counter which now gives time since power on (the on-chip batteries have failed; this used to keep incrementing with power off)
 int LastCT1=-1;//The last command received by the payload, expressed as a decimal number (it is in HeX on the GUI display)
 int CountCT1=-1;//Count of commands received by the payload since power on

 //Barometer information: NOT INTERPRETED from line CT1
 float Baro1T=-999;//Barometer 1 Temperature
 float Baro1P=-999;//Barometer 1 Pressure
 float Baro2T=-999;//Barometer 2 Temperature
 float Baro2P=-999;//Barometer 2 Pressure
 //Barometer information: INTERPRETED from line CT1
 float TempB1=-999;//Barometer 1 Temperature
 float TempB2=-999;//Barometer 2 Temperature
 float PressB1=-999;//Barometer 1 Pressure
 float PressB2=-999;//Barometer 2 Pressure

 float GOCT1=-999;
 float coinCT1=-999;

 //Voltages
 float Volt5VCT1=-999;  // Positive 5V from line CT1
 float Volt15VCT1=-999; // Positive 15V from line CT1

 //Data FROM "CT3" LINE
 int yCT3=-1;//Year from CT3 line linked to the event (last read CT3 line)
 int mCT3=-1;//Month from CT3 line linked to the event (last read CT3 line)
 int dCT3=-1;//Day from CT3 line linked to the event (last read CT3 line)
 int hCT3=-1;//Hour from CT3 line linked to the event (last read CT3 line)
 int miCT3=-1;//Minute from CT3 line linked to the event (last read CT3 line)
 int sCT3=-1;//Second from CT3 line linked to the event (last read CT3 line)
 float TempCT3=-999; //Temperature measured on the board CT3

 int OnTimeCT3=-1;//1/second counter which now gives time since power on (the on-chip batteries have failed; this used to keep incrementing with power off)
 int LastCT3=-1;//The last command received by the payload, expressed as a decimal number (it is in HeX on the GUI display)
 int CountCT3=-1;//Count of commands received by the payload since power on
 //Voltages
 float Volt5VCT3=-999;  // Positive 5V from line CT3
 float Volt15VCT3=-999; // Positive 15V from line CT3
 //TRIGGER RATES (PHA AND LOGIC) from CT3
 float T1L=-1;
 float T1A=-1;
 float T2L=-1;
 float T2A=-1;
 float T3L=-1;
 float T3A=-1;
 float T4L=-1;
 float T4A=-1;
 float GRDL=-1;
 float GRDA=-1;

  //HOUSEKEEPING FROM POW
 int yPOW=-1;
 int mPOW=-1;
 int dPOW=-1;
 int hPOW=-1;
 int miPOW=-1;
 int sPOW=-1;
 int OnTimePOW=-1;
 float MainC=-999;
 float MainV=-999;
 float HeatC=-999;
 float HeatV=-999;
 float TrackC=-999;
 float TrackV=-999;

 //From VCI line
 int yVCI=-1;
 int mVCI=-1;
 int dVCI=-1;
 int hVCI=-1;
 int miVCI=-1;
 int sVCI=-1;
 int* LrVCI=new int[7];
 for (int ijk=0;ijk<7;ijk++)LrVCI[ijk]=0;

 int NCT1words=25; //CT1 line
 int NCT3words=21; //CT3 line
 int NPOWwords=16; //POW line
 int NVCIwords=33; //VCI line

 //number of words per line
 //int NPHAwords=13; //PHA line for MIPFlit3.30
 //int NPHAwords=11; //PHA line for MIPFlit3.31
 int NPHAwords=12; //PHA line for MIPFlit3.72

 int NEVTwords=0;
 if(EVTLength==0)  NEVTwords=4; //EVT line from MIPFlit3.30 to 3.35
 if(EVTLength==1)  NEVTwords=16; //Was increased at the same time of ASIC Length
 if(EVTLength==2)  NEVTwords=18; //Was increased at the same time of ASIC Length from NL0155

 int NASIwords= 4; //ASI line

//Number of Layers read
 int kASI=0;
 int kEVT=0;
 int kPHA=0;
 int flagGoodGo=0;
 int NoisyClus=0;


//Start to read the file
 for (string line; getline(filestr,line);)
   {
    // if(kEvent>3) return -1;
    //Here the full line is read and stored in the string object "line"
    //if(kEvent>10) continue;

    kLine++;
    // if(kLine%1000==0)
    //  cout<< "Line " << kLine << " is read" <<endl;

///////////////////////////////////////////////
//Compare the first 3 char of the line to "CT1"
///////////////////////////////////////////////

    if (line.compare(0,3,"CT1",3)==0)
     {
      ///////////////////////////////////////////////
      //Extract the data from the last CT1 line
      ///////////////////////////////////////////////
      vector<string> datalineCT1;
      //Split the line
      datalineCT1=split(&line,' ');

      if((int)datalineCT1.size()!=NCT1words)
        {
         cout << "Number of words in the CT1 line is wrong : " << (int)datalineCT1.size() << " (Should be " << NCT1words  << "): "<< endl;
         cout << line <<endl;
         continue;
        }
      //Extract CT1 date
      //cout << "extract CT1 date" <<endl;
      //cout << dataline.at(1) <<endl;
      extractdate(&yCT1,&mCT1,&dCT1,&datalineCT1.at(1));
      //cout << "extract CT1 date is done" <<endl;
      //cout << "extract CT1 time" <<endl;
      extracttime(&hCT1,&miCT1,&sCT1,&datalineCT1.at(2));

      //Format from Paul's email of Feb. 8, 2018

      OnTimeCT1=s2i(&datalineCT1.at(3));
      LastCT1=s2i(&datalineCT1.at(4));
      CountCT1= s2i(&datalineCT1.at(5));
      Baro1T=s2f(&datalineCT1.at(8));
      Baro1P=s2f(&datalineCT1.at(9));
      Baro2T=s2f(&datalineCT1.at(10));
      Baro2P=s2f(&datalineCT1.at(11));
      GOCT1=s2f(&datalineCT1.at(14));
      coinCT1=s2f(&datalineCT1.at(15));
      PressB1=s2f(&datalineCT1.at(16));
      TempB1=s2f(&datalineCT1.at(17));
      PressB2=s2f(&datalineCT1.at(18));
      TempB2=s2f(&datalineCT1.at(19));
      TempCT1=s2f(&datalineCT1.at(20));
      Volt5VCT1=s2f(&datalineCT1.at(21));
      Volt15VCT1=s2f(&datalineCT1.at(22));
     }
///////////////////////////////////////////////
//Compare the first 3 char of the line to "CT3"
///////////////////////////////////////////////

    if (line.compare(0,3,"CT3",3)==0)
     {
      ///////////////////////////////////////////////
      //Extract the data from the last CT3 line
      ///////////////////////////////////////////////
      vector<string> datalineCT3;
      //Split the line
      datalineCT3=split(&line,' ');

      if((int)datalineCT3.size()!=NCT3words)
        {
         cout << "Number of words in the CT3 line is wrong : " << (int)datalineCT3.size() << " (Should be " << NCT3words  << "): "<< endl;
         cout << line <<endl;
         continue;
        }

      //Extract CT3 date
      //cout << "extract CT3 date" <<endl;
      //cout << dataline.at(1) <<endl;
      extractdate(&yCT1,&mCT1,&dCT3,&datalineCT3.at(1));
      //cout << "extract CT3 date is done" <<endl;
      //cout << "extract CT3 time" <<endl;
      extracttime(&hCT3,&miCT3,&sCT3,&datalineCT3.at(2));

      //Format from Paul's email of Feb. 8, 2018

      OnTimeCT3=s2i(&datalineCT3.at(3));
      LastCT3=s2i(&datalineCT3.at(4));
      CountCT3= s2i(&datalineCT3.at(5));
      T1L=s2f(&datalineCT3.at(6));
      T1A=s2f(&datalineCT3.at(7));
      T2L=s2f(&datalineCT3.at(8));
      T2A=s2f(&datalineCT3.at(9));
      T3L=s2f(&datalineCT3.at(10));
      T3A=s2f(&datalineCT3.at(11));
      T4L=s2f(&datalineCT3.at(12));
      T4A=s2f(&datalineCT3.at(13));
      GRDL=s2f(&datalineCT3.at(14));
      GRDA=s2f(&datalineCT3.at(15));
      TempCT3=s2f(&datalineCT3.at(16));
      Volt5VCT3=s2f(&datalineCT3.at(17));
      Volt15VCT3=s2f(&datalineCT3.at(18));

     }
///////////////////////////////////////////////
//Compare the first 3 char of the line to "POW"
///////////////////////////////////////////////

    if (line.compare(0,3,"POW",3)==0)
     {
      ///////////////////////////////////////////////
      //Extract the data from the last POW line
      ///////////////////////////////////////////////
      vector<string> datalinePOW;
      //Split the line
      datalinePOW=split(&line,' ');

      if((int)datalinePOW.size()!=NPOWwords)
        {
         cout << "Number of words in the POW line is wrong : " << (int)datalinePOW.size() << " (Should be " << NPOWwords  << "): "<< endl;
         cout << line <<endl;
         continue;
        }

      //Extract datalinePOW date
      //cout << "extract datalinePOW date" <<endl;
      //cout << dataline.at(1) <<endl;
      extractdate(&yPOW,&miPOW,&dPOW,&datalinePOW.at(1));
      //cout << "extract POW date is done" <<endl;
      //cout << "extract POW time" <<endl;
      extracttime(&hPOW,&miPOW,&sPOW,&datalinePOW.at(2));

      //Format from Paul's email of Feb. 8, 2018

      OnTimePOW=s2i(&datalinePOW.at(3));
      MainV=s2f(&datalinePOW.at(4));
      MainC=s2f(&datalinePOW.at(5));
      HeatV=s2f(&datalinePOW.at(6));
      HeatC=s2f(&datalinePOW.at(7));
      TrackV=s2f(&datalinePOW.at(8));
      TrackC=s2f(&datalinePOW.at(9));
     }

///////////////////////////////////////////////
//Compare the first 3 char of the line to "VCI"
///////////////////////////////////////////////

    if (line.compare(0,3,"VCI",3)==0)
     {
      ///////////////////////////////////////////////
      //Extract the data from the last VCI line
      ///////////////////////////////////////////////
      vector<string> datalineVCI;
      //Split the line
      datalineVCI=split(&line,' ');

      if((int)datalineVCI.size()!=NVCIwords)
       {
        cout << "Number of words in the VCI line is wrong : " << (int)datalineVCI.size() << " (Should be " << NVCIwords  << "): "<< endl;
        cout << line <<endl;
        continue;
       }

      //Extract datalineVCI date
      //cout << "extract datalineVCI date" <<endl;
      //cout << dataline.at(1) <<endl;
      extractdate(&yVCI,&miVCI,&dVCI,&datalineVCI.at(1));
      //cout << "extract POW date is done" <<endl;
      //cout << "extract POW time" <<endl;
      extracttime(&hVCI,&miVCI,&sVCI,&datalineVCI.at(2));

      //Format from Paul's email of April 26, 2019
      for(int ijk=0;ijk<7;ijk++) LrVCI[ijk]=s2i(&datalineVCI.at(19+ijk));

     }


///////////////////////////////////////////////
//Compare the first 3 char of the line to "PHA"
///////////////////////////////////////////////

    if (line.compare(0,3,"PHA",3)==0)
     {
      //cout << "In PHA line" <<endl;
      //Each PHA Line creates one event
      //If not the first PHA Line fill the previous event in the tree
      if(kEvent!=0)
       {
        int t=0;
        int tc=0;
        for (int ij=0;ij<7;ij++) t+=Titmp[ij]*(int)TMath::Power(2,ij);
        for (int ij=0;ij<7;ij++) tc+=Tictmp[ij]*(int)TMath::Power(2,ij);
        e->set_Nhnoisy(NoisyClus);
        e->set_Ti(t);
        e->set_Tic(tc);
        tree->Fill();
        //reset event
        kEVT=0;
        kASI=0;
        kPHA=0;
        flagGoodGo=0;
        NoisyClus=0;
        delete e;
       }

      //Define new event
      e=new ALEvent();
      //Increment the number of read events
      kEvent++;
      kPHA=1;
      e->set_eventnumber(kEvent);

      //Define internal tracker trigger
      Titmp=new int[7];
      for(int ij=0;ij<7;ij++)Titmp[ij]=0;
      Tictmp=new int[7];
      for(int ij=0;ij<7;ij++)Tictmp[ij]=0;
      ///////////////////////////////////////////////
      //Extract the data from the last PHA line
      ///////////////////////////////////////////////
      vector<string> datalinePHA;

      //Split the line
      datalinePHA=split(&line,' ');

      if((int)datalinePHA.size()!=NPHAwords)
        {
         cout << "Number of words in the PHA line is wrong : " << (int)datalinePHA.size() << " (Should be " << NPHAwords  << "): "<< endl;
         cout << line <<endl;
         continue;
        }
      int year=0;
      int month=0;
      int day=0;
      int hour=0;
      int minute=0;
      int second=0;

      //Extract PHA date
      year=month=day=0;
      //cout << "extract PHA date" <<endl;
      //cout << dataline.at(1) <<endl;
      extractdate(&year,&month,&day,&datalinePHA.at(1));
      //cout << "extract PHA date is done" <<endl;
      //Extract PHA time
      hour=minute=second=0;
      //cout << "extract PHA time" <<endl;
      extracttime(&hour,&minute,&second,&datalinePHA.at(2));
      //cout << "extract PHA time is done" <<endl;
      //Fill the event structure with PHA time variables
      e->set_yPHA(year);
      e->set_mPHA(month);
      e->set_dPHA(day);
      e->set_hPHA(hour);
      e->set_miPHA(minute);
      e->set_sPHA(second);
      //cout << "PHA time was added to the event" <<endl;
      //Fill the trigger information from the PHA line
      e->add_EneT1(s2lf(&datalinePHA.at(3)));
      e->add_EneT2(s2lf(&datalinePHA.at(4)));
      e->add_EneT3(s2lf(&datalinePHA.at(5)));
      e->add_EneT4(s2lf(&datalinePHA.at(6)));
      e->add_Eneg(s2lf(&datalinePHA.at(7)));
      e->add_PHA6(s2lf(&datalinePHA.at(8)));
      e->set_GoPHA(s2lf(&datalinePHA.at(9)));
      e->set_tPHA(s2lf(&datalinePHA.at(10)));


      //Set the housekeeping variables (It can be defautl value at the beginning of a run)
      //CT1
      e->set_yCT1(yCT1);
      e->set_mCT1(mCT1);
      e->set_dCT1(dCT1);
      e->set_hCT1(hCT1);
      e->set_miCT1(miCT1);
      e->set_sCT1(sCT1);
      e->set_OnTimeCT1(OnTimeCT1);
      e->set_LastCT1(LastCT1);
      e->set_CountCT1(CountCT1);
      e->set_TempCT1(TempCT1);
      e->set_Volt5VCT1(Volt5VCT1);
      e->set_Volt15VCT1(Volt15VCT1);
      e->set_Baro1T(Baro1T);
      e->set_Baro1P(Baro1P);
      e->set_Baro2T(Baro2T);
      e->set_Baro2P(Baro2P);
      e->set_GOCT1(GOCT1);
      e->set_coinCT1(coinCT1);
      e->set_PressB1(PressB1);
      e->set_TempB1(TempB1);
      e->set_PressB2(PressB2);
      e->set_TempB2(TempB2);


      //CT3
      e->set_yCT3(yCT3);
      e->set_mCT3(mCT3);
      e->set_dCT3(dCT3);
      e->set_hCT3(hCT3);
      e->set_miCT3(miCT3);
      e->set_sCT3(sCT3);
      e->set_OnTimeCT3(OnTimeCT3);
      e->set_LastCT3(LastCT3);
      e->set_CountCT3(CountCT3);
      e->set_TempCT3(TempCT3);
      e->set_Volt5VCT3(Volt5VCT3);
      e->set_Volt15VCT3(Volt15VCT3);
      e->set_T1L(T1L);
      e->set_T1A(T1A);
      e->set_T2L(T2L);
      e->set_T2A(T2A);
      e->set_T3L(T3L);
      e->set_T3A(T3A);
      e->set_T4L(T4L);
      e->set_T4A(T4A);
      e->set_GRDL(GRDL);
      e->set_GRDA(GRDA);

      //POW
      e->set_yPOW(yPOW);
      e->set_mPOW(mPOW);
      e->set_dPOW(dPOW);
      e->set_hPOW(hPOW);
      e->set_miPOW(miPOW);
      e->set_sPOW(sPOW);
      e->set_OnTimePOW(OnTimePOW);
      e->set_MainC(MainC);
      e->set_MainV(MainV);
      e->set_HeatC(HeatC);
      e->set_HeatV(HeatV);
      e->set_TrackC(TrackC);
      e->set_TrackV(TrackV);

      //VCI
      e->set_yVCI(yVCI);
      e->set_mVCI(mVCI);
      e->set_dVCI(dVCI);
      e->set_hVCI(hVCI);
      e->set_miVCI(miVCI);
      e->set_sVCI(sVCI);

      for(int ijk=0;ijk<7;ijk++)  e->set_Lrate(ijk,LrVCI[ijk]);

     }//end read line starting with PHA

///////////////////////////////////////////////
//Compare the first 3 char of the line to "EVT"
///////////////////////////////////////////////
    if (line.compare(0,3,"EVT",3)==0)
     {
     // cout << "In EVT line" <<endl;

      //An EVT line is always associated to the previous PHA line
      //If several EVT lines follow one PHA line then the last EVT line is recorded
      ///////////////////////////////////////////////
      //Extract the data from the EVT line
      ///////////////////////////////////////////////
      vector<string> datalineEVT;
      kEVT++;
      //Split the line
      datalineEVT=split(&line,' ');
      double oldGo=prevGo;
      string oldline=prevEVTline;
      prevEVTline=line;
      if((int)datalineEVT.size()!=NEVTwords)
       {
        cout << "Number of words in the EVT line is wrong : " << (int)datalineEVT.size() << " (Should be " << NEVTwords  << "): "<< endl;
        cout << line <<endl;
       }

       if(EVTLength==1)
        {
         prevGo=s2lf(&datalineEVT.at(11));
        }
       if(EVTLength==2)
        {
         prevGo=s2lf(&datalineEVT.at(13));
        }

     // if(prevGo!=(oldGo+1))continue;
      flagGoodGo=1;
      //cout << oldline << endl;
     // cout << line <<endl;
      //cout << "Line: "<< kLine<<", oldGO=" << oldGo << ", prevGo=" <<prevGo<<endl;

      if(kPHA!=1)continue;
      if(kEVT!=1) {continue;} //Does not fill if the EVT line doesnt follow directly a PHA line

      int year=0;
      int month=0;
      int day=0;
      int hour=0;
      int minute=0;
      int second=0;
      //Extract EVT date
      //cout << "extract EVT date" <<endl;
      //cout << dataline.at(1) <<endl;
      extractdate(&year,&month,&day,&datalineEVT.at(1));
      //cout << "extract EVT date is done" <<endl;
      //Extract PHA time
      hour=minute=second=0;
      //cout << "extract EVT time" <<endl;
      extracttime(&hour,&minute,&second,&datalineEVT.at(2));

      //cout << "extract EVT time is done" <<endl;
      //Fill the event structure with EVT time variables
      e->set_yEVT(year);
      e->set_mEVT(month);
      e->set_dEVT(day);
      e->set_hEVT(hour);
      e->set_miEVT(minute);
      e->set_sEVT(second);
      //cout << "EVT time was added to the event" <<endl;
      //Fill data
      e->set_EVT(datalineEVT.at(3));
      if(EVTLength==1)
       {
        e->set_GoEVT(s2lf(&datalineEVT.at(11)));
        e->set_tEVT(s2lf(&datalineEVT.at(13)));
       }
      if(EVTLength==2)
       {
        e->set_nHitLEVT(s2i(&datalineEVT.at(5)));
        e->set_CCEVT(s2i(&datalineEVT.at(7)));
        e->set_PatternEVT(s2i(&datalineEVT.at(9)));
        e->set_Q1EVT(s2i(&datalineEVT.at(11)));
        e->set_GoEVT(s2lf(&datalineEVT.at(13)));
        e->set_tEVT(s2lf(&datalineEVT.at(15)));
        e->set_TrigEVT(s2lf(&datalineEVT.at(17)));
       }
        //cout << "GoEVT "<<s2f(&datalineEVT.at(11)) << endl;
      // cout << "EVT data was added to the event" <<endl;
     }//end read line starting with EVT

///////////////////////////////////////////////
//Compare the first 3 char of the line to "ASI"
///////////////////////////////////////////////
    if (line.compare(0,3,"ASI",3)==0)
     {
      //cout << "In ASI line" <<endl;
      if (kASI==7) continue;//Skip the ASI lines after 7 lines
      //Define vector of string "dataline"
      //The vector will contains the words of the line
      vector<string> dataline;
      dataline=split(&line,' '); //split the line in a vector of words
      kASI++;
      if((int)dataline.size()!=NASIwords)
       {
        cout << "Number of words in the ASI line is wrong : " << (int)dataline.size() << " (Should be " << NASIwords  << "): "<< endl;
        cout << line <<endl;
        return -1;
       }
      if(kPHA!=1)continue;
      if(flagGoodGo==0) continue;
      //Extract ASI date
      int yearASI=0;
      int monthASI=0;
      int dayASI=0;
      //cout << "extract ASI date" <<endl;
      //cout << dataline.at(1) <<endl;
      extractdate(&yearASI,&monthASI,&dayASI,&dataline.at(1));
      //cout << "extract ASI date is done" <<endl;
      //Extract ASI time
      int hourASI=0;
      int minuteASI=0;
      int secondASI=0;
      //cout << "extract ASI time" <<endl;
      extracttime(&hourASI,&minuteASI,&secondASI,&dataline.at(2));
      //cout << "extract ASI time is done" <<endl;
      int tmpNoisyClus=0;

      //Decode ASIC code
      if(ASILength==0)DecodeASIShort(dataline.at(3),&Hh,Titmp);
      if(ASILength==1)DecodeASILong(dataline.at(3),&Hh,Titmp,Tictmp,&tmpNoisyClus);

      for(int ij=0;ij<(int)Hh.size();ij++)
        {
         //e->set_L(Hh.at(ij)->getL(),dataline.at(3));
         e->set_flagL(Hh.at(ij)->get_L(),1);
         Hh.at(ij)->set_year(yearASI);
         Hh.at(ij)->set_m(monthASI);
         Hh.at(ij)->set_d(dayASI);
         Hh.at(ij)->set_hour(hourASI);
         Hh.at(ij)->set_mi(minuteASI);
         Hh.at(ij)->set_s(minuteASI);
         e->add_hit(Hh.at(ij));
        }
      NoisyClus+=tmpNoisyClus;


      //Clear hits
      Hh.clear();
     }//end read line starting with ASI
   }//end read file


 //Fill last event
 int t=0;
 int tc=0;
 for (int ij=0;ij<7;ij++) t+=Titmp[ij]*(int)TMath::Power(2,ij);
 for (int ij=0;ij<7;ij++) tc+=Tictmp[ij]*(int)TMath::Power(2,ij);
 e->set_Nhnoisy(NoisyClus);
 e->set_Ti(t);
 e->set_Tic(tc);
 tree->Fill();
 //reset event
 delete e;


 //Write tree in output file
 tree->Write(0,TObject::kOverwrite);

 //Close the file
 filestr.close();
 cout << "The file "<<filename  << " is closed"<< endl;
 fileout->Close();
 cout << "The output file "<<Form("%s.root",input.at(0).c_str())  << " is closed"<< endl;

 return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////
//Decode the hexadecimal ASI data
//////////////////////////////////////////////////////////////////////////////////////////

int DecodeASIShort(string data,vector<ALTckhit*>* Hh,int*Ti)
{

//For data of August 10, 2017

//Get string length

 int datalength=data.length();
 //Check the number of character
 //It must be 5 + a multiple of 3
 //The last character is always a " " so we remove 1 to the datalength
//cout <<endl;
 //cout << "ASI Data length " << data.length();
 //cout << " ASI Data: " << data <<endl;
 if((datalength-5)%3!=0)
  {
   cout << "The length of the data string does not pass the requirement: (datalength-5)%3==0 : Length=" <<datalength <<endl;
   cout << data <<endl;
   //for(int i=0;i<data.size();i++)
   //    cout << data[i]<<endl;
   cout << "No hit is added"<< endl;
   return 0;
  }

///////////////////////////////////////
//The first 5 digits have a static format
///////////////////////////////////////
//Read the first 2 digits: 1 byte (8 bits)
//They corresponds to the FPGA address of the layer

 string FPGAadd=data.substr(0,2);
 int LL=0;
 int L=0;
 LL=stoi(FPGAadd,0,16);
//Convert FPGA address to Layer number from 1 to 7
 switch(LL)
  {
   case 8: L=0;break;
   case 1: L=1;break;
   case 2: L=2;break;
   case 3: L=3;break;
   case 4: L=4;break;
   case 5: L=5;break;
   case 6: L=6;break;
  }//switch
 // cout << "Layer: " << L <<endl;

//Read the 3rd and 4th digit: 1 byte (8 bits)
//From left to right
//Event tag first part: 5 bits: clock counter running in the trackerboard FPGA
//The counters are synched at 0 when a reset/synch signal is received by the FPGA.
//Event tag second part: 2 bits: ASIC trigger tag
//Error flag: 1 bit: is set to 1 in case of trigger tag mismatch

 string Eventtag=data.substr(2,2); //Not used for now

//Read the 5th digit: 4 bits
//Number of chips reporting cluster data

 string sNchip=data.substr(4,1);
 //convert to integer
 int Nchip=stoi(sNchip,0,16);

 //cout << "Nchip: " << Nchip <<endl;
 // If there is not cluster in the layer
 if (Nchip==0) return 1;
///////////////////////////////////////
//From now the format is dynamical
//However it always comes 3 digits per 3 digits
///////////////////////////////////////
 //Define the digits index to read
 int index=0;
 //Start with 5 for the first chip header
 index=5;
 if(index>=datalength||index+3>=datalength) return 0;//Number of bytes incorrect

 //Variables to combine with the chip boundary clusters

 int* laststr=new int [12];
 int* firststr=new int [12];
 for(int i=0;i<12;i++) laststr[i]=-1;  //will be equal to the j index  of the cluster that has a strip 63
 for(int i=0;i<12;i++) firststr[i]=-1;  //will be equal to the j index  of the cluster that has a strip 0

 vector<ALTckhit*>tmpHV;
 ALTckhit*tmpH;

 for(int i=0;i<Nchip;i++)//Loop over the number of chips reporting cluster data
   {

    //Read the 3-digit ASIC header: 12 bits
    //1 bit: cluster overflow
    //1 bit: unused
    //4 bits: number of clusters
    //1 bit: chip error
    //1 bit: parity error
    //4 bits: chip address

    //Read the number of clusters
    //The number of clusters is 4 bits contained in 2 digits
    //The 2 lowest-order bits of the first digit (4 bits) and
    //the 2 highest-order bits of the second digit (4 bits)
    //First read the 2 digits
    if(index>=datalength||index+2>=datalength) return 0;//Number of bytes incorrect
    string sNclus=data.substr(index,2);
    //convert to 16-bit integer
    uint16_t intNclus= (uint16_t)stoi(sNclus,0,16);
    //Shift the 4 bits to the lowest position: shift of 2 is needed
    //then keep only the last 4 bits: Make a AND with 0x000F
    int Nclus=(int)((intNclus >> 2) &0x000F);

    //Read cluster overflow: Shift of 7 bits
    int Clusoverflow=(int)((intNclus >> 7) &0x0001);
    //Read chip error: Shift of 1 bit
    int chiperr=(int)((intNclus >> 1) &0x0001);
    //Read parity: Shift of 0 bit
    int parityerr=(int)(intNclus &0x0001);

    //cout << "Chip error: " << chiperr <<endl;
    //cout << "Parity error: " << parityerr <<endl;

    //if(chiperr==1) continue;
    //if(parityerr==1) continue;

    //Read the chip address  between 0 and 15
    if(index+2>=datalength||index+3>=datalength) return 0; //Number of bytes incorrect
    string sChipadd=data.substr(index+2,1);
    //cout << sChipadd << endl;
    //convert to integer
    int Chipadd=stoi(sChipadd,0,16);

   // cout << "Chip: "<< Chipadd << ", NClus: " << Nclus <<endl;

    index+=3;
    if(index>=datalength||index+3>=datalength) return 0;//Number of bytes incorrect
    for(int j=0;j<Nclus;j++)//Loop over the number of clusters
      {
       tmpH=new ALTckhit();
       //Read the 3-digit cluster data: 12 bits
       //6 bits: number of strips minus 1 (0 to 63)
       //6 bits: address of the first strip (0 to 63)
       if(index>=datalength||index+3>=datalength) return 0;//Number of bytes incorrect
       string sClus=data.substr(index,3);
       //convert to 16-bit integer
      //cout << sClus <<endl;
       if(sClus.empty() || sClus==" ")
        {
         tmpH->set_L(L);
         tmpH->set_chip(Chipadd);
         tmpH->set_nstrips(-1);
         tmpH->set_fstrip(-1);
         tmpH->set_fstripID(-1);
         tmpH->set_overflow(Clusoverflow,-1.);
         tmpH->set_chiperr(chiperr,-1.);
         tmpH->set_parityerr(parityerr,-1.);
         //Fill up the vector
         Hh->push_back(tmpH);
         //Increment the index to the next 3-digit word
         index+=3;
         continue;
        }
       uint16_t Clus;

       Clus=(uint16_t)stoi(sClus,0,16);

       //Mask to read the address of the first strip
       //First shift of 6 bits to the right
       //then read the last 6 bits
       int Nstrip=(int)((Clus >>6) & 0x003F);

       //Mask to read the address of the first strip(last 6 bits): 0x003F
       int Stripadd=(int)(Clus & 0x003F);

       //cout << "Nstrip: "<< Nstrip+1 << ", Stripadd: " << Stripadd <<endl;

       //Determine the strip ID
       //Equation from Sarah's email of September 4 2017.
       //From  0 to 767 : Lowest index strip of the cluster
       //int  firstStripNumber=64*(Chipadd+1)-Stripadd-1-Nstrip;
       int  firstStripNumber=64*(Chipadd+1)-Stripadd-1;

       tmpH->set_L(L);
       tmpH->set_chip(Chipadd);
       tmpH->set_nstrips(Nstrip);
       tmpH->set_fstrip(Stripadd);
       tmpH->set_fstripID(firstStripNumber);

       tmpH->set_overflow(Clusoverflow,-1.);
       tmpH->set_chiperr(chiperr,-1.);
       tmpH->set_parityerr(parityerr,-1.);

       //First  strip from left to right is hit
       //if(Stripadd+Nstrip==63) firststr[Chipadd]=tmpHV.size();
       if(Stripadd==63) firststr[Chipadd]=tmpHV.size();
       //Last strip from left to right is hit
       //if(Stripadd==0) laststr[Chipadd]=tmpHV.size();
       if(Stripadd-Nstrip==0) laststr[Chipadd]=tmpHV.size();

       //Internal trigger for Layer L set to 1
       Ti[L]=1;
       tmpH->set_noisy(0);
       //Mask of bad strips
       for(int k=0;k<Nstrip;k++)
         {
          int tmpstrip=firstStripNumber+k;
          //if(L==3&&tmpstrip==9){tmpH->set_noisy(1);}
          if(L==4&&tmpstrip==357){tmpH->set_noisy(1);}
          if(L==4&&tmpstrip==358){tmpH->set_noisy(1);}
          //if(L==5&&tmpstrip==691){tmpH->set_noisy(1);}
          //if(L==6&&tmpstrip==576){tmpH->set_noisy(1);}
         }

       //Fill up the vector
       tmpHV.push_back(tmpH);

       //Increment the index to the next 3-digit word
       index+=3;
      }//j
   }//i
 //Loop over the clusters and merge at the boundaries if needed
 for(int i=0;i<(int) tmpHV.size();i++)
   {
    int chip=tmpHV.at(i)->get_chip();
    int fstrip=tmpHV.at(i)->get_fstrip();
    int nstrips=tmpHV.at(i)->get_nstrips();
    int noisy=tmpHV.at(i)->get_nstrips();

    if(fstrip==0 && chip<11 && chip!=5) //The last strip is touched not chip 5 nor 11
     {
      int ij= firststr[chip+1];
      //There is a cluster on the first strip of the next chip
      if(ij!=-1)
       {
        if (noisy==0)
          tmpHV.at(i)->set_noisy(tmpHV.at(ij)->get_noisy());

        tmpHV.at(i)->set_nstrips(nstrips+tmpHV.at(ij)->get_nstrips()+1);
        tmpHV.at(i)->set_nstripsNC(tmpHV.at(ij)->get_nstrips());
        tmpHV.at(i)->set_overflow((unsigned int)1,tmpHV.at(ij)->get_overflow(0));
        tmpHV.at(i)->set_chiperr((unsigned int)1,tmpHV.at(ij)->get_chiperr(0));
        tmpHV.at(i)->set_parityerr((unsigned int)1,tmpHV.at(ij)->get_parityerr(0));

       }
      Hh->push_back(tmpHV.at(i));
     }
     else if(fstrip+nstrips==63 && chip>0 && chip!=6) //The first strip is touched not chip 0 nor 6
     {
      int ij= laststr[chip-1];
      //There is not a cluster on the last strip of the previous chip: This is a real hit
      if(ij==-1)
       {
        if(!(noisy==1 && nstrips==0))  //Do not add cluster with 1 single inner strip that is noisy
          Hh->push_back(tmpHV.at(i));
       }
      //If there is a cluster on the last strip of the previous we don't record the information.
      //The information is filled  with the other cluster.
     }
    else //Fill the non boudary clusters
     {
       if(!(noisy==1 && nstrips==0))  //Do not add cluster with 1 single inner strip that is noisy
        Hh->push_back(tmpHV.at(i));
     }
   }//i


 return 1;
}



int DecodeASILong(string data,vector<ALTckhit*>* Hh,int*Ti,int*Tic,int* Nhitnoisy)
{

//For data from November, 2017

//Get string length
 int noisyhit=0;
 int datalength=data.length();
 //Check the number of characters
 //It must be 9 + a multiple of 3
 //The last character is always a " " so we remove 1 to the datalength
//cout <<endl;
 //cout << "ASI Data length " << data.length();
 //cout << " ASI Data: " << data <<endl;
 if((datalength-9)%3!=0)
  {
   cout << "The length of the data string does not pass the requirement: (datalength-9)%3==0 : Length=" <<datalength <<endl;
   cout << data <<endl;
   //for(int i=0;i<data.size();i++)
   //    cout << data[i]<<endl;
   cout << "No hit is added"<< endl;
   return 0;
  }

///////////////////////////////////////
//The first 9 digits have a static format
///////////////////////////////////////

//Read the first 2 digits: 1 byte (8 bits)
//They correspond to length of the following bytes
// It is the lentgh of the ASIC string + 1
 string ASICstrlgth=data.substr(0,2);

//Read the 3rd and 4th digits: 1 byte (8 bits)
//It should always be "E7"
 string Indentifier=data.substr(2,2);


//Read the 5th and 6th digits: 1 byte (8 bits)
//They correspond to the FPGA address of the layer

 string FPGAadd=data.substr(4,2);
 int LL=0;
 int L=0;
 LL=stoi(FPGAadd,0,16);
//Convert FPGA address to Layer number from 1 to 7
 switch(LL)
  {
   case 8: L=0;break;
   case 1: L=1;break;
   case 2: L=2;break;
   case 3: L=3;break;
   case 4: L=4;break;
   case 5: L=5;break;
   case 6: L=6;break;
  }//switch
  //cout << "Layer: " << L <<endl;

//Read the 7th and 8th digits: 1 byte (8 bits)
//From left to right
//Event tag first part: 5 bits: clock counter running in the trackerboard FPGA
//The counters are synched at 0 when a reset/synch signal is received by the FPGA.
//Event tag second part: 2 bits: ASIC trigger tag
//Error flag: 1 bit: is set to 1 in case of trigger tag mismatch

 string Eventtag=data.substr(6,2); //Not used for now

//Read the 9th digit: 4 bits
//Number of chips reporting cluster data

 string sNchip=data.substr(8,1);
 //convert to integer
 int Nchip=stoi(sNchip,0,16);

// cout << "Nchip: " << Nchip <<endl;
 // If there is not cluster in the layer
 if (Nchip==0) return 1;
///////////////////////////////////////
//From now the format is dynamical
//However it always comes 3 digits per 3 digits
///////////////////////////////////////
 //Define the digits index to read
 int index=0;
 //Start with 9 for the first chip header
 index=9;


 if(index>=datalength||index+3>=datalength) return 0;//Number of bytes incorrect


 //Variables to combine with the chip boundary clusters

 int* laststr=new int [12];
 int* firststr=new int [12];
 for(int i=0;i<12;i++) laststr[i]=-1;  //will be equal to the j index  of the cluster that has a strip 63
 for(int i=0;i<12;i++) firststr[i]=-1;  //will be equal to the j index  of the cluster that has a strip 0

 vector<ALTckhit*>tmpHV;
 ALTckhit*tmpH;
 for(int i=0;i<Nchip;i++)//Loop over the number of chips reporting cluster data
   {

    //Read the 3-digit ASIC header: 12 bits
    //1 bit: cluster overflow
    //1 bit: unused
    //4 bits: number of clusters
    //1 bit: chip error
    //1 bit: parity error
    //4 bits: chip address

    //Read the number of clusters
    //The number of clusters is 4 bits contained in 2 digits
    //The 2 lowest-order bits of the first digit (4 bits) and
    //the 2 highest-order bits of the second digit (4 bits)
    //First read the 2 digits
    if(index>=datalength||index+2>=datalength) return 0;//Number of bytes incorrect
    string sNclus=data.substr(index,2);
    //convert to 16-bit integer
    uint16_t intNclus= (uint16_t)stoi(sNclus,0,16);
    //Shift the 4 bits to the lowest position: shift of 2 is needed
    //then keep only the last 4 bits: Make a AND with 0x000F
    int Nclus=(int)((intNclus >> 2) &0x000F);

    //Read cluster overflow: Shift of 7 bits
    int Clusoverflow=(int)((intNclus >> 7) &0x0001);
    //Read chip error: Shift of 1 bit
    int chiperr=(int)((intNclus >> 1) &0x0001);
    //Read parity overflow: Shift of 0 bit
    int parityerr=(int)(intNclus &0x0001);

    //cout << "Chip error: " << chiperr <<endl;
    //cout << "Parity error: " << parityerr <<endl;

    //if(chiperr==1) continue;
    //if(parityerr==1) continue;

    //Read the chip address  between 0 and 15
    if(index+2>=datalength||index+3>=datalength) return 0; //Number of bytes incorrect
    string sChipadd=data.substr(index+2,1);
    //cout << sChipadd << endl;
    //convert to integer
    int Chipadd=stoi(sChipadd,0,16);

    //cout << "Chip: "<< Chipadd << ", NClus: " << Nclus <<endl;

    index+=3;
    if(index>=datalength||index+3>=datalength) return 0;//Number of bytes incorrect


    for(int j=0;j<Nclus;j++)//Loop over the number of clusters
      {
       tmpH=new ALTckhit();
       //Read the 3-digit cluster data: 12 bits
       //6 bits: number of strips minus 1 (0 to 63)
       //6 bits: address of the first strip (0 to 63)
       if(index>=datalength||index+3>=datalength) return 0;//Number of bytes incorrect
       string sClus=data.substr(index,3);
       //convert to 16-bit integer
       //cout << sClus <<endl;
       if(sClus.empty() || sClus==" ")
        {
         tmpH->set_L(L);
         tmpH->set_chip(Chipadd);
         tmpH->set_nstrips(-1);
         tmpH->set_fstrip(-1);
         tmpH->set_fstripID(-1);
         tmpH->set_overflow(Clusoverflow,-1.);
         tmpH->set_chiperr(chiperr,-1.);
         tmpH->set_parityerr(parityerr,-1.);
        //Fill up the vector
         //Hh->push_back(tmpH);
         //Increment the index to the next 3-digit word
         index+=3;
         continue;
        }
       uint16_t Clus;

       Clus=(uint16_t)stoi(sClus,0,16);

       //Mask to read the address of the first strip
       //First shift of 6 bits to the right
       //then read the last 6 bits
       int Nstrip=(int)((Clus >>6) & 0x003F);

       //Mask to read the address of the first strip(last 6 bits): 0x003F
       int Stripadd=(int)(Clus & 0x003F);

      // cout << "Nstrip: "<< Nstrip+1 << ", Stripadd: " << Stripadd <<endl;

       //Determine the strip ID
       //Equation from Sarah's email of September 4 2017.
       //From 0 to 767
       int  firstStripNumber=64*(Chipadd+1)-Stripadd-1-Nstrip;
       //int  firstStripNumber=64*(Chipadd+1)-Stripadd-1;

       tmpH->set_L(L);
       tmpH->set_chip(Chipadd);
       tmpH->set_nstrips(Nstrip);
       tmpH->set_fstrip(Stripadd);
       tmpH->set_fstripID(firstStripNumber);

       tmpH->set_overflow(Clusoverflow,-1.);
       tmpH->set_chiperr(chiperr,-1.);
       tmpH->set_parityerr(parityerr,-1.);

       //First  strip from left to right is hit
       if(Stripadd+Nstrip==63) firststr[Chipadd]=tmpHV.size();
       //if(Stripadd==63) firststr[Chipadd]=tmpHV.size();
       //Last strip from left to right is hit
       if(Stripadd==0) laststr[Chipadd]=tmpHV.size();
       //if(Stripadd-Nstrip==0) laststr[Chipadd]=tmpHV.size();

       //Internal trigger for Layer L set to 1
       Ti[L]=1;
       tmpH->set_noisy(0);

       //Flag is 1: Don't use the hit
       int flagN=0;
       if(L==4&&firstStripNumber==357&&Nstrip<=1)flagN=1;
       if(L==4&&firstStripNumber==358&&Nstrip==0)flagN=1;

       //Mask of bad strips
       for(int k=0;k<Nstrip+1;k++)
         {
          //For flight configuration 2018
	        int tmpstrip=firstStripNumber+k;
          if(L==4&&tmpstrip==357){tmpH->set_noisy(1);}
          if(L==4&&tmpstrip==358){tmpH->set_noisy(1);}
         }

       //Fill up the vector
       if(flagN==0) tmpHV.push_back(tmpH);

       //Increment the index to the next 3-digit word
       index+=3;
      }//j
   }//i

 //Loop over the clusters and merge at the boundaries if needed
 for(int i=0;i<(int) tmpHV.size();i++)
   {
    int chip=tmpHV.at(i)->get_chip();
    int fstrip=tmpHV.at(i)->get_fstrip();
    int noisy=tmpHV.at(i)->get_noisy();
    int nstrips=tmpHV.at(i)->get_nstrips();

    if(fstrip==0 && chip<11 && chip!=5) //The last strip is touched not chip 5 nor 11
    //if(fstrip-nstrips==0 && chip<11 && chip!=5) //The last strip is touched not chip 5 nor 11
     {
      int ij= firststr[chip+1];
      //There is a cluster on the first strip of the next chip
      if(ij!=-1)
       {
        if(noisy==0){ tmpHV.at(i)->set_noisy(noisy+tmpHV.at(ij)->get_noisy());noisy=noisy+tmpHV.at(ij)->get_noisy();}
        nstrips=nstrips+tmpHV.at(ij)->get_nstrips()+1;
        tmpHV.at(i)->set_nstrips(nstrips);
        tmpHV.at(i)->set_nstripsNC(tmpHV.at(ij)->get_nstrips());
        tmpHV.at(i)->set_overflow((unsigned int)1,tmpHV.at(ij)->get_overflow(0));
        tmpHV.at(i)->set_chiperr((unsigned int)1,tmpHV.at(ij)->get_chiperr(0));
        tmpHV.at(i)->set_parityerr((unsigned int)1,tmpHV.at(ij)->get_parityerr(0));
       }
      Hh->push_back(tmpHV.at(i));
      if(nstrips<3)Tic[L]=1;
      noisyhit+=noisy;
     }
    else if(fstrip+nstrips==63 && chip>0 && chip!=6) //The first strip is touched not chip 0 nor 6
     //else if(fstrip==63 && chip>0 && chip!=6) //The first strip is touched not chip 0 nor 6
     {
      int ij= laststr[chip-1];
      //There is not a cluster on the last strip of the previous chip: This is a real hit
      if(ij==-1)
       {
        Hh->push_back(tmpHV.at(i));
        if(nstrips<3)Tic[L]=1;
        noisyhit+=noisy;
       }

      //If there is a cluster on the last strip of the previous we don't record the information.
      //The information is filled  with the other cluster.
     }
    else //Fill the non boudary clusters
     {
      Hh->push_back(tmpHV.at(i));
      if(nstrips<3)Tic[L]=1;
      noisyhit+=noisy;

     }
   }//i
 *Nhitnoisy=noisyhit;
 return 1;
}



//*******************************************************************
// Split the string str into a vector of string
// The separator is the character between each word of the string str
// Each string of the vector contains a word of the string
// The vector is ordered from the first to the last word
// The vector index starts from 0
//*******************************************************************
vector<string> split (string* str,char separator)
{

  vector<string> output;
  string::size_type prev_pos=0, pos=0;

  while ((pos =str->find(separator,pos))!=string::npos)
    {
      string substring(str->substr(prev_pos,pos-prev_pos));
      if(!substring.empty())output.push_back(substring);//Check that substring is not empty before filling the vectors
      prev_pos=++pos;
    }
  output.push_back(str->substr(prev_pos,pos-prev_pos)); //Last word
  return output;
}

//*******************************************************************
// Transform the string str into a double
//*******************************************************************
double s2lf(string* str)
{
  istringstream buffer(*str);
  double temp;
  buffer >> temp;
  return temp;
}

//*******************************************************************
// Transform the string str into a float
//*******************************************************************
float s2f(string* str)
{
  istringstream buffer(*str);
  float temp;
  buffer >> temp;
  return temp;
}

//*******************************************************************
// Transform the string str to int
//*******************************************************************
int s2i(string* str)
{
  istringstream buffer(*str);
  int temp;
  buffer >> temp;
  return temp;

}


//*******************************************************************
// Two functions for sorting in ascending order a array
//*******************************************************************

/*Function for partitioning the array*/

int Partition(int low,int high,float arr[])
{
 float high_vac,low_vac,pivot;
 pivot=arr[low];
 while(high>low)
   {
    high_vac=arr[high];
  //  cout << "high_vac : " << high_vac << " pivot " << pivot<< endl;
    while(pivot<high_vac)
     {
      if(high<=low) break;
      high--;
      high_vac=arr[high];
     }
    arr[low]=high_vac;
    low_vac=arr[low];
    while(pivot>=low_vac)// modification
     {
      if(high<=low) break;
      low++;
      low_vac=arr[low];
     }
    arr[high]=low_vac;
   }
 arr[low]=pivot;
 return low;
}

void Quick_sort(int low,int high,float arr[])
{
 //cout << "in Quick_sort " << low <<" " << high<<endl;
 int Piv_index;
 if(low<high)
  {
   Piv_index=Partition(low,high,arr);
   // cout << "Piv_index " << Piv_index <<endl;

   Quick_sort(low,Piv_index-1,arr);
   Quick_sort(Piv_index+1,high,arr);
  }
}



//*******************************************************************
// Extract year,month and day from the string str
//*******************************************************************
void extractdate(int *y,int*m,int*d,string*str)
{
	  vector<string> date;
	  date=split(str,'/'); //separator is '/'
	  *y=s2i(&date.at(2));
	  *m=s2i(&date.at(0));
	  *d=s2i(&date.at(1));
}

//*******************************************************************
// Extract hour,minute and seconde from the string str
//*******************************************************************
void extracttime(int *h,int*m,int*s,string*str)
{
	  vector<string> time;
	  time=split(str,':');// separator is ':'
	  *h=s2i(&time.at(0));
	  *m=s2i(&time.at(1));
	  *s=s2i(&time.at(2));
}
