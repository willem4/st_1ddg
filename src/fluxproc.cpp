#include <vector>
#include <iostream>
#include <cmath>

#include "state.h"
#include "element.h"
#include "node.h"
#include "flux.h"

//#include "output.h"
#include "fluxproc.h"
#include "enumtypes.h"
#include "updatenodes.h"
#include "physpara.h"

using namespace std;

int fluxproc(vector<node> &Nodes, vector<element>& Elements, int bcl, int bcr, testtypedescr testtype, fluxbase * fluxtype, double t){
  int i,j;
  double max;
  double ws;
  
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

  double * F1 = new double [element::systemsize1];
  double * F2 = new double [element::systemsize1];
  double * F  = new double [element::systemsize1];
  double * Fl = new double [element::systemsize1];
  double * Fr = new double [element::systemsize1];
  double * FL = new double [element::systemsize1];
  double * FR = new double [element::systemsize1];
   
  //Fill boundary conditions and fill q1,q2,qL,qR for each node.
  //u1,u2,uL,uR remain unchanged after flux2proc.
  updatenodesq(Nodes, Elements, bcl, bcr, testtype, t);

  //Initialise flux vectors to prevent garbage output when no fluxprocedure is called.
  for (j = 0; j < element::systemsize1; j++){
    Fl[j]=0.0;  
    Fr[j]=0.0;
    F[j]=0.0;
    F1[j]=0.0;
    F2[j]=0.0;
    FL[j]=0.0;
    FR[j]=0.0;

  }

  //  Calculate integral of F in each element and add to RHS.
  if (state::statesize == 2){
    //Loop over Elements
    for (i = 0; i < Elements.size(); i++){
      //Get the function values and write to F1 and F2;
      fluxtype->eval(F1,Elements[i].u1,Elements[i].u1,Elements[i].q1,Elements[i].q1);
      fluxtype->eval(F2,Elements[i].u2,Elements[i].u2,Elements[i].q2,Elements[i].q2);
      for (j = 0; j < element::systemsize1; j++){
	//Add 3*integral of F in element to RHS vector.  
	Elements[i].RHSofU[j].coefficient[1]+=(F1[j]+F2[j])*3.0;
      }
      //Initialise the Krivodonva Discontinuity Detector;
      Elements[i].Kriv = 0.0;
      Elements[i].KrivCount = 0;
    }
  }
  // Calculate fluxes for internal faces.
  // Loop over the internal nodes. 
  for (i = 1; i < Nodes.size()-1; i++){
    //Calculate the conservative flux at the node i
    ws = fluxtype->eval(F,Nodes[i].uL,Nodes[i].uR,Nodes[i].qL,Nodes[i].qR);
    //    cout << ws << "\n";
//     if (ws > Elements[i-1].ws) {
//       Elements[i-1].ws=ws;
//     }
//     if (ws > Elements[i].ws) {
//       Elements[i].ws=ws;
//     }

    //For the Stabilisation operator.
    fluxtype->eval(FL,Nodes[i].uL,Nodes[i].uL,Nodes[i].qL,Nodes[i].qL);
    fluxtype->eval(FR,Nodes[i].uR,Nodes[i].uR,Nodes[i].qR,Nodes[i].qR);
    //Calculate the flux at the node at the left side of the right element.
    fluxtype->eval_l(Fl,Nodes[i].uL,Nodes[i].uR,Nodes[i].qL,Nodes[i].qR);
    //Calculate the flux at the node at the right side of the left element.
    fluxtype->eval_r(Fr,Nodes[i].uL,Nodes[i].uR,Nodes[i].qL,Nodes[i].qR);
    //Add F, Fl and Fr to RHSofU vector.
    for (j = 0; j < element::systemsize1; j++){
      Elements[i-1].RHSofU[j].coefficient[0]-=(F[j]+Fr[j]);
      Elements[i].RHSofU[j].coefficient[0]+=(F[j]+Fl[j]);
    }
    if (state::statesize == 2){
      for (j = 0; j < element::systemsize1; j++){
	Elements[i-1].RHSofU[j].coefficient[1]-=(F[j]+Fr[j])*3.0;
	Elements[i].RHSofU[j].coefficient[1]-=(F[j]+Fl[j])*3.0;

	//Stabilisation operator.
	Elements[i-1].StabOp.coefficient[1]+=abs((FL[j]+Fr[j])-(FR[j]+Fl[j]));
	Elements[i].StabOp.coefficient[1]+=abs((FL[j]+Fr[j])-(FR[j]+Fl[j]));
      }
      //Krivodonova discontinuity detector. 
      //      cout << (FL[0]+Fr[0])-(FR[0]+Fl[0]) <<".";
      if (F[0]+Fr[0] < 0) {
	Elements[i-1].Kriv+=Nodes[i].uL[0]-Nodes[i].uR[0];
	Elements[i-1].KrivCount+=1;
      }
      if (F[0]+Fl[0] > 0) {
	Elements[i].Kriv+=Nodes[i].uR[0]-Nodes[i].uL[0];
	Elements[i].KrivCount+=1;
      }
    }
  }
  
  // Now calculate flux at Nodes[0] (left boundary).
  i = 0; 
  //Calculate the conservative flux at the node i
  fluxtype->eval(F,Nodes[i].uL,Nodes[i].uR,Nodes[i].qL,Nodes[i].qR);
//   if (ws > Elements[i].ws) {
//     Elements[i].ws=ws;
//   }
  //For the Stabilisation operator.
  fluxtype->eval(FL,Nodes[i].uL,Nodes[i].uL,Nodes[i].qL,Nodes[i].qL);
  fluxtype->eval(FR,Nodes[i].uR,Nodes[i].uR,Nodes[i].qR,Nodes[i].qR);
  
  //Calculate the flux at the node at the left side of the right element.
  fluxtype->eval_l(Fl,Nodes[i].uL,Nodes[i].uR,Nodes[i].qL,Nodes[i].qR);
  //Calculate the flux at the node at the right side of the left element.
  fluxtype->eval_r(Fr,Nodes[i].uL,Nodes[i].uR,Nodes[i].qL,Nodes[i].qR);
  //Add F, Fl and Fr to RHS vector.
  for (j = 0; j < element::systemsize1; j++){
    Elements[i].RHSofU[j].coefficient[0]+=(F[j]+Fl[j]);
  }
  if (state::statesize == 2){
    for (j = 0; j < element::systemsize1; j++){
      Elements[i].RHSofU[j].coefficient[1]-=(F[j]+Fl[j])*3.0;
      
      //Stabilisation operator.
      Elements[i].StabOp.coefficient[1]+=abs((FL[j]+Fr[j])-(FR[j]+Fl[j]));
    }
    //Krivodonova discontinuity detector. 
    if (F[0]+Fl[0] > 0) {
      Elements[i].Kriv+=Nodes[i].uR[0]-Nodes[i].uL[0];
      Elements[i].KrivCount+=1;
      //      Elements[i].KrivDiv+=F[0]+Fl[0];
    }    
  }

  // Now calculate flux at Node[Nodes.size()] (right boundary).
  i = Nodes.size()-1; 
  //Calculate the conservative flux at the node i
  fluxtype->eval(F,Nodes[i].uL,Nodes[i].uR,Nodes[i].qL,Nodes[i].qR);
  //   if (ws > Elements[i-1].ws) {
  //     Elements[i-1].ws=ws;
  //   }
  //For the Stabilisation operator.
  fluxtype->eval(FL,Nodes[i].uL,Nodes[i].uL,Nodes[i].qL,Nodes[i].qL);
  fluxtype->eval(FR,Nodes[i].uR,Nodes[i].uR,Nodes[i].qR,Nodes[i].qR);

  //Calculate the flux at the node at the left side of the right element.
  fluxtype->eval_l(Fl,Nodes[i].uL,Nodes[i].uR,Nodes[i].qL,Nodes[i].qR);
  //Calculate the flux at the node at the right side of the left element.
  fluxtype->eval_r(Fr,Nodes[i].uL,Nodes[i].uR,Nodes[i].qL,Nodes[i].qR);
  //Add F, Fl and Fr to RHS vector.
  for (j = 0; j < element::systemsize1; j++){
    Elements[i-1].RHSofU[j].coefficient[0]-=(F[j]+Fr[j]);
  }
  if (state::statesize == 2){
    for (j = 0; j < element::systemsize1; j++){
      Elements[i-1].RHSofU[j].coefficient[1]-=(F[j]+Fr[j])*3.0;
      
      //Stabilisation operator.
      Elements[i-1].StabOp.coefficient[1]+=abs((FL[j]+Fr[j])-(FR[j]+Fl[j]));
    } 
    //Krivodonova discontinuity detector. 
    if (F[0]+Fr[0] < 0) {
      Elements[i-1].Kriv+=Nodes[i].uL[0]-Nodes[i].uR[0];
      Elements[i-1].KrivCount+=1;
      //      Elements[i-1].KrivDiv-=(F[0]+Fr[0]);
    }
  }
  if (state::statesize == 2){
    //Loop over Elements
    for (i = 0; i < Elements.size(); i++){
      //Scale the stabilisation operator.
      Elements[i].StabOp.coefficient[1]=dissipation_const*dissipation_beta*12.0*Elements[i].StabOp.coefficient[1]/pow(Elements[i].length,3.0-dissipation_alpha);
      //Scale the Krivodonova discontinuity detector.
      Elements[i].Kriv=abs(Elements[i].Kriv)/(pow(Elements[i].length/2.0,state::statesize/2.0)*Elements[i].U[0].coefficient[0]);
      if (Elements[i].KrivCount == 0) {
	Elements[i].Kriv = 0.0;
      }
      if (Elements[i].KrivCount == 2) {
	Elements[i].Kriv = Elements[i].Kriv/2.0;
      }
      //When solution is smooth, do not apply any stabilisation.
      if (Elements[i].Kriv <= 1.0) {
	Elements[i].StabOp.coefficient[1] = 0.0;
      }
      //      cout << Elements[i].KrivCount <<"";
      //      if (Elements[i].Kriv > 1.0) {
      //      cout << i << ", ";// << Elements[i].StabOp.coefficient[1] << "\n"; //<< Elements[i].Kriv << "  Krivodonova check \n";
      //      }
    }
  }
  
  delete [] FL;
  delete [] FR;
  delete [] Fr;
  delete [] Fl;
  delete [] F;
  delete [] F2;
  delete [] F1;
  return 0;
}

  

