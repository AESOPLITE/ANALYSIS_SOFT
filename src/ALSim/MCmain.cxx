#include "TriggerRE.h"
#include "MakeRawEventMC.h"

////////////////////////////////////////////////////////////////////////////////////////// 
int main(int argc, char* argv[]);
////////////////////////////////////////////////////////////////////////////////////////// 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// Apply Dead time on Counts
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
#ifndef __CINT__
int main(int argc, char* argv[])
{
  
 if(argc!=1)
  {
   cout << "Wrong number of parameters!!!!" << endl;
   cout << "The program needs 0 input parameter:"<< endl;
   return -1;
  }
// MakeRawEventMC();
 //TriggerRE(true);

 return 0;
}
#endif