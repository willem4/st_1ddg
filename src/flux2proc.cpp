#include <vector>
#include <iostream>
#include <math.h>

#include "state.h"
#include "element.h"
#include "node.h"
#include "flux2.h"

#include "enumtypes.h"
#include "flux2proc.h"
#include "updatenodes.h"

using namespace std;

int flux2proc(vector<node> &Nodes, vector<element>& Elements, int bcl, int bcr, testtypedescr testtype, flux2base * flux2type, double t)
{  
  int i,j;
  //   
  //    |<-Element-->|<-Element-->|
  //    |            |            | 
  //    |            |            |
  //   L|R 1      2 L|R  1    2  L|R
  // ======*======*======*====*=====  
  //   r l          r l          r l
  //    |            |            |
  //    |            |            |
  //    ^node        ^node        ^node
  
  double ds;
  double * F1 = new double [element::systemsize1];
  double * F2 = new double [element::systemsize1];
  double * F  = new double [element::systemsize1];

  //Fill boundary conditions and fill u1,u2,uL,uR for each node.
  updatenodesu(Nodes, Elements, bcl, bcr, testtype, t);
  
  //Initialise flux vectors to prevent garbage output when no fluxprocedure is called.
  for (j = 0; j < element::systemsize2; j++){
    F[j]=0;
    F1[j]=0;
    F2[j]=0;
  }

  //  Calculate integral of F in each element and add to RHS.
  if (state::statesize == 2){
    //Loop over Elements
    for (i = 0; i < Elements.size(); i++){
      //Get the function values and write to F1 and F2;
      flux2type->eval(F1,Elements[i].u1,Elements[i].u1);
      flux2type->eval(F2,Elements[i].u2,Elements[i].u2);
      for (j = 0; j < element::systemsize2; j++){
	//Add 3*integral of F in element to RHS vector.  
	Elements[i].RHSofQ[j].coefficient[1]+=(F1[j]+F2[j])*3.0;
      }
    }
  }

  // Calculate fluxes for internal faces.
  // Loop over the internal nodes. 
  for (i = 1; i < Nodes.size()-1; i++){
    //Calculate the conservative flux at the node i
    ds = flux2type->eval(F,Nodes[i].uL,Nodes[i].uR);
    if (ds > Elements[i-1].dds) {
      Elements[i-1].dds=ds;
    }
    if (ds > Elements[i].dds) {
      Elements[i].dds=ds;
    }

    //Add F to RHSofQ vector.
    for (j = 0; j < element::systemsize2; j++){
      Elements[i-1].RHSofQ[j].coefficient[0]-=F[j];
      Elements[i].RHSofQ[j].coefficient[0]+=F[j];
    }
    if (state::statesize == 2){
      for (j = 0; j < element::systemsize2; j++){
	Elements[i-1].RHSofQ[j].coefficient[1]-=F[j]*3.0;
	Elements[i].RHSofQ[j].coefficient[1]-=F[j]*3.0;
      }
    }
  }
  
  // Now calculate flux at Nodes[0] (left boundary).
  i = 0; 
  //Calculate the conservative flux at the node i
  ds = flux2type->eval(F,Nodes[i].uL,Nodes[i].uR);
  if (ds > Elements[i].dds) {
    Elements[i].dds=ds;
  }

  //  cout << Nodes[i].uL[0] << Nodes[i].uR[0] << Nodes[i].qL[0] << Nodes[i].qR[0] "\n";
  //Add F, Fl and Fr to RHSofQ vector.
  for (j = 0; j < element::systemsize2; j++){
    Elements[i].RHSofQ[j].coefficient[0]+=F[j];
  }
  if (state::statesize == 2){
    for (j = 0; j < element::systemsize2; j++){
      Elements[i].RHSofQ[j].coefficient[1]-=F[j]*3.0;
    }
  }

  // Now calculate flux at Node[Nodes.size()] (right boundary).
  i = Nodes.size()-1; 
  //Calculate the conservative flux at the node i
  ds = flux2type->eval(F,Nodes[i].uL,Nodes[i].uR);
  if (ds > Elements[i-1].dds) {
    Elements[i-1].dds=ds;
  }
  //Add F to RHSofQ vector.
  for (j = 0; j < element::systemsize2; j++){
    Elements[i-1].RHSofQ[j].coefficient[0]-=F[j];
  }
  if (state::statesize == 2){
    for (j = 0; j < element::systemsize2; j++){
      Elements[i-1].RHSofQ[j].coefficient[1]-=F[j]*3.0;
    }
  }
  delete [] F;
  delete [] F2;
  delete [] F1;
  
}

