#include "MainRawEventMC.h"
using namespace std;


int main() 
{
 //This function transforms the rootfiles from the output of write99toroot to RawEvent 
 //This function reconstructs the raw events file per file 

 //Input files 
 string Inppath="/data/psmangeard/AESOPLite/Fluka/UniB/V1";
 string Inppath2="rootfiles";
 string startfile="aesopliteUniB_V1";
 string endfile="_fort.99";
 
 //Output files
 string Outpath="/home/psmangeard/MCproduction/AESOPLITE/rootfiles/UniB/V1";
 
 
 
 //Fluka type of particle
 int type=3; //3 for electrons
 //Number of energies
 int Nene=12;
 //Energies
 int Ene[12]={10,20,30,40,50,60,70,80,90,100,200,300};
 //Number of cycles per energy
 int* Ncycles=new int[Nene]; 
 for(int i=0;i<Nene;i++)Ncycles[i]=100;  
 
 
 
 for(int i=0;i<Nene;i++)//Energies
   {
    for(int j=0;j<Ncycles[i];j++)//Number of cycles
      {
       cout << Form("%s/%d/%s/%s_%d_%dMeV%03d%s.root",Inppath.c_str(),type,Inppath2.c_str(),startfile.c_str(),type,Ene[i],j+1,endfile.c_str()) <<endl;
       //Create the RawEvent
       MakeRawEventMC(type,Ene[i],j+1,Inppath,Inppath2,Outpath,startfile,endfile);
       //       ntuple =new TNtuple("Track","Track","ncase:mreg:mtrack:type:age:e:x:y:z:cx:cy:cz:Edep:flag");
      }//j
    }//i

 
 
 return 0;
}