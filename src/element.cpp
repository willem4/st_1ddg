#include <iostream>
#include <math.h>

#include "state.h"
#include "element.h"

using namespace std;

int element::systemsize1;
int element::systemsize2;
double element::ep = 1.0/sqrt(3.0);
double element::CFL;

element::element(){
}

element::element(int value1, int value2, double value3){
  Node1 = value1;
  Node2 = value2;
  length = value3;
  U = new state [element::systemsize1];
  RHSofU = new state [element::systemsize1];
  U_old = new state [element::systemsize1];
  RHSofU_old = new state [element::systemsize1];  
  Q = new state [element::systemsize2];
  RHSofQ = new state [element::systemsize2];
  u1 = new double [element::systemsize1];
  u2 = new double [element::systemsize1];
  u1slope = new double [element::systemsize1];
  u2slope = new double [element::systemsize1];
  q1 = new double [element::systemsize1];
  q2 = new double [element::systemsize1];
  state StabOp;
  state StabOp_old;
  double ws;
  double dds;
  double Kriv;
  int KrivCount;
}

void element::updateu(){
  int j;
  //! Loop over the variables.
  for (j = 0; j < element::systemsize1; j++){
    //! Determine the new values at +ep and -ep.
    u1[j]=U[j].get(-ep);
    u2[j]=U[j].get(ep);
    u1slope[j]=2*U[j].getslope(-ep)/length;
    u2slope[j]=2*U[j].getslope(ep)/length;
  }
}

void element::updateq(){
  int j;
  //! Loop over the variables.
  for (j = 0; j < element::systemsize2; j++){
    //! Determine the new values at +ep and -ep.
    q1[j]=Q[j].get(-ep);
    q2[j]=Q[j].get(ep);
  }
}



element::~element(){
}
