#include<iostream>
#include<fstream>
#include<string>

#include"enumtypes.h"

using namespace std;

int readtest2(testtypedescr &testtype, timeintmethod &timeint, string &meshname, double &CFL, double &it, double &et, int &deg, string testname)
{
  int timeno,testno;
  ifstream InFile(testname.c_str(), ios::in);
  InFile >> testno;
  switch (testno)
    {    
    case 0 : testtype = Advection;
      break;
    case 1 : testtype = Burgers;
      break;
    case 2 : testtype = SWEBurgers;
      break;
    case 3 : testtype = SWELinearWaveSolution;
      break;
    case 4 : testtype = SWERiemannProblemLeftRarefactionRightShock;
      break;
    case 5 : testtype = SWERiemannProblemLeftShockRightRarefaction; 
      break;
    case 6 : testtype = SWERiemannProblemLeftShockRightShock; 
      break;
    case 7 : testtype = SWERiemannProblemLeftRarefactionRightRarefaction;
      break;
    case 8 : testtype = SWEContinuousTopography;
      break;
    case 9 : testtype = SWEDiscontinuousTopography;
      break;
    case 10 : testtype = SWEFlowOverIsolatedContinuousRidgeI;
      break;
    case 11 : testtype = SWEFlowOverIsolatedDiscontinuousRidgeI;
      break;
    case 12 : testtype = SWEFlowOverIsolatedContinuousRidgeIV;
      break;
    case 13 : testtype = SWEFlowOverIsolatedDiscontinuousRidgeIV;
      break;
    case 14 : testtype = BurgersDiffusive;
      break;
    case 15 : testtype = SedimentTransport;
      break;
    case 16 : testtype = HLLGrass;
      break;
    case 17 : testtype = LFGrass;
      break;
    case 18 : testtype = LFGrassMomentum;
      break;
    case 19 : testtype = DecoupledGrassSedimentOnly;
      break;
    case 20 : testtype = SuspendedSediment;
      break;
    case 21 : testtype = LFGrassMomentumDiffusion;
      break;
    case 22 : testtype = SWEFlowOverSinyBedI;
      break;
    }
  InFile >> timeno;
  switch (timeno)
    {    
    case 0 : timeint = EulerForward;
      break;
    case 1 : timeint = RungeKutta3;
      break;
    case 2 : timeint = CrankNicolsonWithPredictor;
      break;
    }
  InFile >> meshname;
  InFile >> CFL;
  InFile >> it;
  InFile >> et;
  InFile >> deg;
  InFile.close();
}

int writetest(testtypedescr &testtype, timeintmethod &timeint, string &meshname, double &CFL, double &it, double &et, int &deg)
{
  ofstream OutFile("testfile.temp", ios::out);
  OutFile << testtype;
  OutFile << timeint;  
  OutFile << meshname;
  OutFile << CFL;
  OutFile << it;
  OutFile << et;
  OutFile << deg;
  OutFile.close();
}

int readtest(testtypedescr &testtype, timeintmethod &timeint, string &meshname, double &CFL, double &it, double &et, int &deg)
{
  string testfilename;
  cout << "Enter the testfile name: \n";
  cin  >> testfilename;
  readtest2(testtype,timeint,meshname,CFL,it,et,deg,testfilename.c_str());
  writetest(testtype,timeint,meshname,CFL,it,et,deg);
}




