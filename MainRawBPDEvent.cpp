////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, September 8, 2017
////////////////////////////////////////////////////////////////////////////////////////// 





#include "MainRawBPDEvent.h"
using namespace std;


int main(int argc, char* argv[]);

    
//This function read and parse the BPD files
//This function reconstructs the raw events file per file 
//The first paramter contains the file paths 
//This will create one root file per BPD file
//The x,y,z coordinates won't be calculated here
int main(int argc, char* argv[])
{
//Check number of parameters 
 if(argc!=2)
  {
   cout << "Wrong number of parameters!!!!" << endl;
   cout << "The program needs 1 input parameter:"<< endl;
   cout << "The first parameter is the filename"<< endl;
   cout << "that contains the BPD files paths"<< endl;
   return -1;
  }
 
//Open input file
 ifstream filestr;
 filestr.open(argv[1]);
 if(!filestr.is_open())
  {
   cout << "The file " << argv[1] << " is not open ... "<< endl ;
   return -1;
  }   
 else
   cout << "The file "<<argv[1] << " is open"<< endl;
 

  //Start to read the file
 for (string line; getline(filestr,line);)
   {
    //Here the full line is read and stored in the string object "line" 
    cout << line <<endl; //Print the file path
    
    //First step
    //Read the text file .BPD
    //Parse the information to create events
    //The events are stored in a root file line.root 
    MakeRawBPDEvent(line);
    
    //Temporary Second Step
    //Read the rootfile previously created
    //Determine the coordinates of the clusters
    //in the X,Y,Z coordinates
    //this is where misalignement must be taken into account
    //MakeEventData(line);


    
    

   }//line 
 
 filestr.close();
 
 return 0;
}
