#include <vector>
#include <iostream>
#include <math.h>
#include <stdlib.h>

#include "state.h"
#include "element.h"
#include "node.h"
#include "flux.h"
#include "testfcn.h"

//#include "output.h"
#include "updatenodes.h"
#include "enumtypes.h"

using namespace std;

int updatenodesu(vector<node> &Nodes, vector<element> &Elements, int bcl, int bcr, testtypedescr testtype, double t)
{
  int i,j;

  for (i = 0; i < Elements.size(); i++){
    //!Update interpolation points in the Element.
    Elements[i].updateu();
  }
  //Loop over internal nodes.
  for (i = 1; i < Nodes.size()-1; i++){  
    //Loop over variables in first system.
    for (j = 0; j < element::systemsize1; j++){ 
      //Get left and right values at the node and write to uL and uR.
      Nodes[i].uL[j]=Elements[i-1].U[j].get(1.0);   
      Nodes[i].uR[j]=Elements[i].U[j].get(-1.0);
    }
  }    

  //Prescribed boundary condtions 
  double etaE, hE, uE, mE;
  
  //Get left and right boundary face values from Elements.
  for (j = 0; j < element::systemsize1; j++){
    Nodes[0].uR[j]=Elements[0].U[j].get(-1.0);
    Nodes[Nodes.size()-1].uL[j]=Elements[Nodes.size()-2].U[j].get(1.0);
  }
  
  //Periodic Boundary conditions 
  if (bcl == 1){
    if (bcr != 1){
      cout << "Right boundary condition should be periodic as well\n";
      exit (-1);
    }
    for (j = 0; j < element::systemsize1; j++){
      Nodes[0].uL[j]=Nodes[Nodes.size()-1].uL[j];
      Nodes[Nodes.size()-1].uR[j]=Nodes[0].uR[j];
    }
  }
  
  //Transmissive Boundary conditions;
  if (bcl == 2){
    for (j = 0; j < element::systemsize1; j++){
      Nodes[0].uL[j]=Nodes[0].uR[j];
    }
  }
  if (bcr == 2){
    for (j = 0; j < element::systemsize1; j++){
      Nodes[Nodes.size()-1].uR[j]=Nodes[Nodes.size()-1].uL[j];
    }
  }
  //Boundary Conditions for Flow over Isolated Ridge.
  if (bcl == 4){
    for (j = 0; j < element::systemsize1; j++){
      Nodes[0].uL[j]=testfcn(testtype,j,Nodes[0].x,t);
    }
  }
  
//     switch (testtype)
//       {
//       case SWERiemannProblemLeftRarefactionRightShock:
//       case DecoupledGrassSedimentOnly:
// 	for (j = 0; j < element::systemsize1; j++){
// 	  Nodes[0].uL[j]=testfcn(testtype,j,Nodes[0].x,t);
// 	}
// 	break;
//       case SWEFlowOverIsolatedContinuousRidgeI:
//       case SWEFlowOverIsolatedDiscontinuousRidgeI:
// 	hE = 0.4;
// 	uE = 0.3962;
// 	etaE = -0.2;
// 	Nodes[0].uL[0]=hE;
// 	Nodes[0].uL[1]=hE*uE;
// 	Nodes[0].uL[2]=etaE;
// 	break;
//       case SWEFlowOverIsolatedContinuousRidgeIV:
//       case SWEFlowOverIsolatedDiscontinuousRidgeIV:
// 	hE = 0.4;
// 	uE = 3.7637;
// 	etaE = -0.2;
// 	Nodes[0].uL[0]=hE;
// 	Nodes[0].uL[1]=hE*uE;
// 	Nodes[0].uL[2]=etaE;
// 	break;
//       case BurgersDiffusive:
// 	uE = testfcn(testtype,0,Nodes[0].x,t);
// 	Nodes[0].uL[0]=uE;
// 	break;
//       case SedimentTransport:
//       case HLLGrass:
//       case LFGrass:
//       case LFGrassMomentum:
//       case LFGrassMomentumDiffusion:
// 	hE = testfcn(testtype,0,Nodes[0].x);
// 	mE = testfcn(testtype,1,Nodes[0].x);
// 	etaE = testfcn(testtype,2,Nodes[0].x);
// 	Nodes[0].uL[0]=hE;
// 	Nodes[0].uL[1]=mE;
// 	Nodes[0].uL[2]=etaE;
// 	break;	
//       default:
// 	cout << "U: Left prescribed contitions not defined for this testcase.\n";
// 	exit(-1);
// 	break;
//       }
  
  if (bcr == 4){
    for (j = 0; j < element::systemsize1; j++){
      Nodes[Nodes.size()-1].uR[j]=testfcn(testtype,j,Nodes[Nodes.size()-1].x,t);
    }
  }

//     switch (testtype)
//       {
//       case SWERiemannProblemLeftRarefactionRightShock:
//       case DecoupledGrassSedimentOnly:
// 	for (j = 0; j < element::systemsize1; j++){
// 	  Nodes[Nodes.size()-1].uR[j]=testfcn(testtype,j,Nodes[Nodes.size()-1].x,t);
// 	}
// 	break;
//       case SWEFlowOverIsolatedContinuousRidgeI:
//       case SWEFlowOverIsolatedDiscontinuousRidgeI:
// 	hE = 0.4;
// 	uE = 0.3962;
// 	etaE = -0.2;
// 	Nodes[Nodes.size()-1].uR[0]=hE;
// 	Nodes[Nodes.size()-1].uR[1]=hE*uE;
// 	Nodes[Nodes.size()-1].uR[2]=etaE;
// 	break;
//       case SWEFlowOverIsolatedContinuousRidgeIV:
//       case SWEFlowOverIsolatedDiscontinuousRidgeIV:
// 	hE = 0.4;
// 	uE = 3.7637;
// 	etaE = -0.2;
// 	Nodes[Nodes.size()-1].uR[0]=hE;
// 	Nodes[Nodes.size()-1].uR[1]=hE*uE;
// 	Nodes[Nodes.size()-1].uR[2]=etaE;
// 	break;
//       case BurgersDiffusive:
// 	uE = testfcn(testtype,0,Nodes[Nodes.size()-1].x,t);
// 	Nodes[Nodes.size()-1].uR[0]=uE;
// 	break;
//       case SedimentTransport:
//       case HLLGrass:
//       case LFGrass:
//       case LFGrassMomentum:
//       case LFGrassMomentumDiffusion:
// 	hE = testfcn(testtype,0,Nodes[Nodes.size()-1].x);
// 	mE = testfcn(testtype,1,Nodes[Nodes.size()-1].x);
// 	etaE = testfcn(testtype,2,Nodes[Nodes.size()-1].x);
// 	Nodes[Nodes.size()-1].uR[0]=hE;
// 	Nodes[Nodes.size()-1].uR[1]=mE;
// 	Nodes[Nodes.size()-1].uR[2]=etaE;
// 	break;	
//       default:
//  	cout << "U: Right prescribed contitions not defined for this testcase.\n";
//  	exit(-1);
//  	break;
//       }

}

int updatenodesq(vector<node> &Nodes, vector<element> &Elements, int bcl, int bcr, testtypedescr testtype, double t)
{
  if (element::systemsize2 > 0){
    int i,j;
    double qE;
    for (i = 0; i < Elements.size(); i++){
      //Update interpolation points in the Element.
      Elements[i].updateq();
    }
    //Loop over internal nodes.
    for (i = 1; i < Nodes.size()-1; i++){  
      //Loop over variables in second system.
      for (j = 0; j < element::systemsize2; j++){ 
	//Get left and right values at the node and write to uL and uR.
	Nodes[i].qL[j]=Elements[i-1].Q[j].get(1.0);   
	Nodes[i].qR[j]=Elements[i].Q[j].get(-1.0);
      }
    }    

    //Get left and right boundary face values from Elements.
    for (j = 0; j < element::systemsize2; j++){
      Nodes[0].qR[j]=Elements[0].Q[j].get(-1.0);
      Nodes[Nodes.size()-1].qL[j]=Elements[Nodes.size()-2].Q[j].get(1.0);
    }
    
  //Periodic Boundary conditions 
    if (bcl == 1){
      if (bcr != 1){
	cout << "Right boundary condition should be periodic as well\n";
	exit (-1);
      }
      for (j = 0; j < element::systemsize2; j++){
      Nodes[0].qL[j]=Nodes[Nodes.size()-1].qL[j];
      Nodes[Nodes.size()-1].qR[j]=Nodes[0].qR[j];
      }
    }
    
    //Transmissive Boundary conditions;
    if (bcl >= 2){
      for (j = 0; j < element::systemsize2; j++){
	Nodes[0].qL[j]=Nodes[0].qR[j];
      }
    }
    if (bcr >= 2){
      for (j = 0; j < element::systemsize2; j++){
	Nodes[Nodes.size()-1].qR[j]=Nodes[Nodes.size()-1].qL[j];
      }
    }
    //Prescribed Boundary conditions;                                                                                                
    if (bcl == 4){
      switch (testtype)
      {
      case SWEFlowOverSinyBedI:
        //Update discharge 
        Nodes[0].qL[1]=testfcn(testtype,1,Nodes[0].x,0);
      default:
        break;
      }  
    }
    if (bcr == 4){
      switch (testtype)
	{
        //Update water depth 
	case SWEFlowOverSinyBedI:
	  Nodes[0].qL[0]=testfcn(testtype,1,Nodes[Nodes.size()-1].x,0);
	default:
	  break;
	}
    }

    //Boundary Conditions.
//     if (bcl == 4){
//       switch (testtype)
// 	{
// 	case BurgersDiffusive:
// 	  qE = 0.0;
// 	  Nodes[0].qL[0]=qE;
// 	  break;
// 	default:
// 	  cout << "Q: Left prescribed contitions not defined for this testcase.\n";
// 	  exit(-1);
// 	  break;
// 	}
//     }
//     if (bcr == 4){
//       switch (testtype)
// 	{
// 	case BurgersDiffusive:
// 	  qE = 0.0;
// 	  Nodes[Nodes.size()-1].qR[0]=qE;
// 	  break;
//       default:
// 	cout << "Q: Right prescribed contitions not defined for this testcase.\n";
// 	exit(-1);
// 	break;
// 	}
//     } 
  }
}
