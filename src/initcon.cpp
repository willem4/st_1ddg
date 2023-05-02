#include <vector>
#include <iostream>
#include <cmath>

#include "initcon.h"
#include "testfcn.h"

#include "enumtypes.h"
#include "element.h"
#include "node.h" 
#include "state.h"

using namespace std;

int initcon(vector<node> Nodes, vector<element>& Elements, testtypedescr testtype, double it){
  double lambda;
  double a, b, xi1, xi2;
  int i,j;
  state temp;
  switch (testtype){
  case SWEContinuousTopography: 
  case SWEFlowOverIsolatedContinuousRidgeI: 
  case SWEFlowOverIsolatedContinuousRidgeIV:
  case SWEFlowOverSinyBedI:
    lambda = 1.0;
    break;
  default:    
    lambda =(sqrt(3.0)-1.0)/(2.0*sqrt(3.0));//.5773502690
    break;
  }
  for(i=0; i<Elements.size(); i++){
    a = Nodes[Elements[i].Node1].x;      
    b = Nodes[Elements[i].Node2].x;
    xi1 = a+lambda*(b-a);
    xi2 = b-lambda*(b-a);
    //    cout << "check: " << xi1 << "\t" << xi2 << "\n";
    //    cout << "check: " << testfcn(testtype,0,xi1) << "\t" << testfcn(testtype,0,xi2) << "\n";
    //    cout << "check: " << testfcn(testtype,1,xi1) << "\t" << testfcn(testtype,1,xi2) << "\n";
    //    cout << "--------------------------------------\n";
    for (j=0; j<element::systemsize1; j++){
      temp.coefficient[0]=(testfcn(testtype,j,xi1,it)+testfcn(testtype,j,xi2,it))/2;
      //      cout << "checkit: " << j <<  "  " << temp.coefficient[0] << " \n";
      if (state::statesize == 2)
	{
	  temp.coefficient[1]=(testfcn(testtype,j,xi2,it)-testfcn(testtype,j,xi1,it))/(xi2-xi1)*(b-a)/2;
	}
      Elements[i].U[j] = temp;
    }
    Elements[i].Kriv = 1.0;
  }
}
