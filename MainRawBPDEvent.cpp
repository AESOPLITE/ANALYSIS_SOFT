////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, September 8, 2017
////////////////////////////////////////////////////////////////////////////////////////// 

#include "LoadDataparameters.h"
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
 if(argc!=6)
  {
   cout << "Wrong number of parameters!!!!" << endl;
   cout << "The program needs 5 input parameters:"<< endl;
   cout << "The first parameter is the filename that contains the BPD files paths"<< endl;
   cout << "The second parameter is the geometry configuration file version"<< endl;
   cout << "The third parameter is an integer to set the KF reconstruction in uni of non-uni B-field (uni = 1, non-uni= 0)" << endl;
   cout << "The fourth parameter is an integer  to set the number of iterations in the KF (true = two iterations, false = one iteration) " << endl;
   cout << "The fifth parameter is the reconstruction ID flag" << endl;
   return -1;
  }
 

//Open input file
 ifstream filestr;
 filestr.open(argv[1]);
 int geoconf=stoi(argv[2]); 
 int FieldConf = stoi(argv[3]);
 bool TwoIter = stoi(argv[4]);
 string RecoID = argv[5];


 if(!filestr.is_open())
  {
   cout << "The file " << argv[1] << " is not open ... "<< endl ;
   return -1;
  }   
 else
   cout << "The file "<<argv[1] << " is open"<< endl;
 
 //Load configuration parameter
 float* zL=new float[7];
 float*OffsetLL=new float[7];
 float*OffsetRL=new float[7];
 float*TrigThresh=new float[5];
 for(int i=0;i<7;i++)zL[i]=OffsetLL[i]=OffsetRL[i]=0;
 for(int i=0;i<5;i++)TrigThresh[i]=0;
 string paramfile=Form("../src/ALSim/Dataparameters%d.dat",geoconf); 
 LoadDataparameters(paramfile,zL,OffsetLL,OffsetRL,TrigThresh);
 for(int i=0;i<7;i++)
   {
    cout << "L"<<i <<", zL:" << zL[i] ;
    cout << ", OffsetLL:" << OffsetLL[i] ;
    cout << ", OffsetRL:" << OffsetRL[i] << endl;
   }  
 cout << "T1 threshold: " << TrigThresh[0] <<endl;
 cout << "T2 threshold: " << TrigThresh[1] <<endl;
 cout << "T3 threshold: " << TrigThresh[2] <<endl;
 cout << "T4 threshold: " << TrigThresh[3] <<endl;
 cout << "Guard threshold: " << TrigThresh[4] <<endl;
  //Start to read the file
 for (string line; getline(filestr,line);)
   {
    //Here the full line is read and stored in the string object "line" 
    cout << line <<endl; //Print the file path
    
    //First step
    //Read the text file .BPD
    //Parse the information to create events
    //The events are stored in a root file line.root 
    MakeRawBPDEventIT(line);
    
    //Temporary Second Step
    //Read the rootfile previously created
    //Determine the coordinates of the clusters
    //in the X,Y,Z coordinates
    //this is where misalignement must be taken into account
    MakeEventData(line,geoconf,zL,OffsetLL,OffsetRL,TrigThresh,FieldConf,TwoIter,RecoID);


   }//line 
 
	cout << "Done processing data" << endl;
 filestr.close();
 
 return 0;
}
