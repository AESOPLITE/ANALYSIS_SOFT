#include "MainRawEventAtmoMC.h"
#include "headers.h"
using namespace std;


//This function transforms the rootfiles from the output of write99toroot to RawEvent 
//This function reconstructs the raw events file per file 

int main(int argc, char*argv[]) 
{

 if(argc!=5)
  {
   cout << "Wrong number of parameters!!!!" << endl;
   cout << "The program needs 4 input parameters:"<< endl;
   cout << "First is Fluka type of simulated particles" <<endl;
   cout << "Second is the atmospheric layer simulated" <<endl;
   cout << "Third is first cycle to reconstruct (Starts at 1)" <<endl;
   cout << "Fourth is the number of cycle to reconstruct" <<endl;
   return -1;
  }
 //Fluka type of particle
 int type=(int) atoi(argv[1]); //3 for electrons
 int layer=(int) atoi(argv[2]);     //atmospheric layer
 int Ncycles=(int) atoi(argv[3]);   //first cycle to reconstruct
 int Ncycles2=(int) atoi(argv[4]);   //last cycle

 //Input files 
 string Inppath="/data/smechbal/Fluka/NonUniB/V4SP";
 string Inppath2="rootfiles";
 string startfile="aesopliteNonUniB_V4SP";
 string endfile="_fort.99";
 
 //Output files
 string Outpath="/home/smechbal/MCproduction/AESOPLITE/rootfiles/NonUniB/V4SP";
 

//for(int j=0;j<1;j++)//Number of cycles
 for(int j=Ncycles;j<Ncycles+Ncycles2;j++)//Number of cycles
      {
       cout << Form("%s/%d/%s/%s_%d_%d_%03d%s.root",Inppath.c_str(),type,Inppath2.c_str(),startfile.c_str(),type,layer,j,endfile.c_str()) <<endl;
	
       	
       //Create the RawEvent
     MakeRawEventAtmoMC(type,layer,j,Inppath,Inppath2,Outpath,startfile,endfile);

      }//j

 
 return 0;
}

