////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, September 8, 2017
////////////////////////////////////////////////////////////////////////////////////////// 

#include "MakeRawBPDEvent.h"

int MakeRawBPDEvent(string filename)
{

 //Open input file
 ifstream filestr;
 filestr.open(filename);
 if(!filestr.is_open())
  {
   cout << "The file " <<filename << " is not open ... "<< endl ;
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
 TFile*fileout=new TFile(Form("%s.root",filename.c_str()),"RECREATE");
 cout << "Output file is created" <<endl;

 // Create a TTree
 TTree *tree = new TTree("BPD","BPD Raw event");
 ALEvent *e = new ALEvent();
 
 // Create a branch with event
 tree->Branch("event",&e);  

//Temporary variable for internal triggers 
 int*Titmp=new int[7];
 for (int i=0;i<7;i++)Titmp[i]=0;
 
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
 int NPHAwords=11; //PHA line for MIPFlit3.31
 int NEVTwords= 4; //EVT line
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
//      cout << "In PHA line" <<endl;
      //Check if the number of PHA line in memory is lower than NPHAmem
      if((int)prevPHAline.size()==NPHAMem)//Remove the first item of the vector as it is the oldest line
        prevPHAline.erase(prevPHAline.begin());
      //Add the current line to the vectors
      prevPHAline.push_back(line);                    
     }//end read line starting with PHA  
         
///////////////////////////////////////////////
//Compare the first 3 char of the line to "EVT"
///////////////////////////////////////////////
    if (line.compare(0,3,"EVT",3)==0) 
     {
//      cout << "In EVT line" <<endl;
      prevEVTline=line;
     }//end read line starting with EVT  

///////////////////////////////////////////////
//Compare the first 3 char of the line to "ASI"
///////////////////////////////////////////////
    if (line.compare(0,3,"ASI",3)==0&&(int)prevPHAline.size()>0) 
     {
//      cout << "In ASI line" <<endl;
      //Define vector of string "dataline"
      //The vector will contains the words of the line 
      vector<string> dataline;
      dataline=split(&line,' '); //split the line in a vector of words
      
      if((int)dataline.size()!=NASIwords)
       {
        cout << "Number of words in the PHA line is wrong : " << (int)dataline.size() << " (Should be " << NASIwords  << "): "<< endl;
        cout << line <<endl;
        return -1;
       }    
      if(kASI==0)//First line of the event
       {
        //Define the object for a possible new event
        e=new ALEvent();
        //Define internal tracker trigger 
        Titmp=new int[7];
        for(int ij=0;ij<7;ij++)Titmp[ij]=0;
        //Extract ASI date
        int yearASI=0;
        int monthASI=0;
        int dayASI=0;
//	cout << "extract ASI date" <<endl;
//	cout << dataline.at(1) <<endl;
        extractdate(&yearASI,&monthASI,&dayASI,&dataline.at(1));
//	cout << "extract ASI date is done" <<endl;
        //Extract ASI time
        int hourASI=0;
        int minuteASI=0;
        int secondASI=0;
//	cout << "extract ASI time" <<endl;
        extracttime(&hourASI,&minuteASI,&secondASI,&dataline.at(2));
//	cout << "extract ASI time is done" <<endl;

        ///////////////////////////////////////////////
        //Extract the data from the last PHA line
        ///////////////////////////////////////////////  
        vector<string> datalinePHA;
        //Get the last index of the vector of PHA lines
        int PHAindex=(int)prevPHAline.size();
        if (PHAindex<1)
         {
          cout << "There is no PHA line in memory" <<endl;    
          cout << "Check the file before the ASI line"<<endl;
          cout << "The ASI line is not taken into account"<<endl;
          cout << line << endl;         
          return -1;
        }
        string tmpPHAline=prevPHAline.at(PHAindex-1);
        //Split the last line of the list of PHA lines
        datalinePHA=split(&tmpPHAline,' '); 
      
        if((int)datalinePHA.size()!=NPHAwords)
         {
          cout << "Number of words in the PHA line is wrong : " << (int)datalinePHA.size() << " (Should be " << NPHAwords  << "): "<< endl;
          cout << prevPHAline.at(PHAindex-1) <<endl;
          return -1;
         }    
        int year=0;
        int month=0;
        int day=0;
        int hour=0;
        int minute=0;
        int second=0;

        //Extract PHA date
        year=month=day=0;
//	cout << "extract PHA date" <<endl;
//	cout << dataline.at(1) <<endl;
        extractdate(&year,&month,&day,&datalinePHA.at(1));
//	cout << "extract PHA date is done" <<endl;
        //Extract PHA time
        hour=minute=second=0;
//	cout << "extract PHA time" <<endl;
        extracttime(&hour,&minute,&second,&datalinePHA.at(2));
//	cout << "extract PHA time is done" <<endl;
        //Fill the event structure with PHA time variables
        e->set_yPHA(year);
        e->set_mPHA(month);
        e->set_dPHA(day);
        e->set_hPHA(hour);
        e->set_miPHA(minute);
        e->set_sPHA(second);
        //cout << "PHA time was added to the event" <<endl;
        //Fill the trigger information from the PHA line
        e->add_EneT1(s2i(&datalinePHA.at(3)));
        e->add_EneT2((double)s2i(&datalinePHA.at(4)));
        e->add_EneT3((double)s2i(&datalinePHA.at(5)));
        e->add_EneT4((double)s2i(&datalinePHA.at(6)));
        e->add_Eneg((double)s2i(&datalinePHA.at(7)));
        e->set_GoPHA((int)s2i(&datalinePHA.at(9)));
        e->set_tPHA((int)s2i(&datalinePHA.at(10)));

        //cout << "Trigger PHA height were added to the event" <<endl;
        ///////////////////////////////////////////////
        //Extract the data from the EVT line
        ///////////////////////////////////////////////  
        vector<string> datalineEVT;
        datalineEVT=split(&prevEVTline,' '); 
      
        if((int)datalineEVT.size()!=NEVTwords)
         {
          cout << "Number of words in the PHA line is wrong : " << (int)datalineEVT.size() << " (Should be " << NEVTwords  << "): "<< endl;
          cout << prevEVTline <<endl;
          return -1;
         }    
        //Extract PHA date
        year=month=day=0;
//	cout << "extract EVT date" <<endl;
//	cout << dataline.at(1) <<endl;
        extractdate(&year,&month,&day,&datalineEVT.at(1));
//	cout << "extract EVT date is done" <<endl;
        //Extract PHA time
        hour=minute=second=0;
//	cout << "extract EVT time" <<endl;
        extracttime(&hour,&minute,&second,&datalineEVT.at(2));
//	cout << "extract EVT time is done" <<endl;
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
       // cout << "EVT data was added to the event" <<endl;

       }//if kASI==0
       
      //Fill data from the layer       
      e->set_L(kASI,dataline.at(3));
      DecodeASI(dataline.at(3),&Hh,Titmp);   
            
      kASI++;
      if(kASI==7)// The seven ASI lines have been read. Let's make the output
       {
        //Increment the number of read events 
        kEvent++;
        e->set_eventnumber(kEvent);
        //Fill internal trigger
        int t=0;
        for (int ij=0;ij<7;ij++) t+=Titmp[ij]*(int)TMath::Power(2,ij);
        e->set_Ti(t);
        //Fill hits
        cout << "Nhits: " << (int)Hh.size() << endl;
        for(int ij=0;ij<(int)Hh.size();ij++)e->add_hit(Hh.at(ij));
        tree->Fill(); 
        delete e;
        Hh.clear();
	//PrintEvent(event,&outfile);
	//Reset the variable kASI 
        kASI=0;

       }       
     }//end read line starting with ASI
   }//end read file 

 
 
 //Write tree in output file
 tree->Write();

 //Close the file
 filestr.close();
 cout << "The file "<<filename  << " is close"<< endl;  
 fileout->Close();

 return 0;
}




////////////////////////////////////////////////////////////////////////////////////////// 
//Decode the hexadecimal ASI data
////////////////////////////////////////////////////////////////////////////////////////// 

int DecodeASI(string data,vector<ALTckhit*>* Hh,int*Ti)
{
 
//For data of August 10, 2017
  
 ALTckhit*tmpH=new ALTckhit();      
    
//Get string length

 int datalength=data.length(); 
 //Check the number of character
 //It must be 5 + a multiple of 3
 //The last character is always a " " so we remove 1 to the datalength
cout <<endl;
 cout << "ASI Data length " << data.length();
 cout << " ASI Data: " << data <<endl;
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
    cout << sChipadd << endl;
    //convert to integer
    int Chipadd=stoi(sChipadd,0,16); 

    cout << "Chip: "<< Chipadd << ", NClus: " << Nclus <<endl;

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
       cout << sClus <<endl;
       if(sClus.empty() || sClus==" ") 
        {
         tmpH->set_L(L);         
         tmpH->set_chip(Chipadd);      
         tmpH->set_nstrips(-1);   
         tmpH->set_fstrip(-1);    
         tmpH->set_fstripID(-1);  
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
       
       cout << "Nstrip: "<< Nstrip+1 << ", Stripadd: " << Stripadd <<endl;
       
       //Determine the strip ID
       //Equation from Sarah's email of Spetember 4 2017.
       int  firstStripNumber=64*(Chipadd+1)-Stripadd;
       
       tmpH->set_L(L);         
       tmpH->set_chip(Chipadd);      
       tmpH->set_nstrips(Nstrip);   
       tmpH->set_fstrip(Stripadd);    
       tmpH->set_fstripID(firstStripNumber);  

       //Internal trigger for Layer L set to 1
       Ti[L]=1;
       tmpH->set_noisy(0);
       //Mask of bad strips
       for(int k=0;k<Nstrip;k++)
         {
          int tmpstrip=firstStripNumber+k;   
          if(L==4&&tmpstrip==359){tmpH->set_noisy(1);}        
          if(L==6&&tmpstrip==576){tmpH->set_noisy(1);}        
         }  
               
       //Fill up the vector 
       Hh->push_back(tmpH);
       
       //Increment the index to the next 3-digit word
       index+=3;
      }//j
   }//i
 

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
	  *m=s2i(&date.at(1));
	  *d=s2i(&date.at(0));	
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




