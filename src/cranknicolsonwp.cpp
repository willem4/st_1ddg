#include <vector>
#include <iostream>

#include "node.h"
#include "state.h"
#include "element.h"
#include "source.h"
#include "flux.h"
#include "flux2.h"

#include "fluxproc.h"
#include "flux2proc.h"
#include "sourceproc.h"

#include "cranknicolsonwp.h"
#include "updatenodes.h"
using namespace std;

int cranknicolsonwp(double &dt, vector<node> &Nodes, vector<element> &Elements, int bcl, int bcr, testtypedescr testtype, fluxbase * fluxtype, sourcebase * sourcetype, flux2base * flux2type, double t)
{ 
  int i,j,k;
  double dttemp;
  //Q1 = R_q(U1) 
  //Update Q1 variables
  flux2proc(Nodes,Elements,bcl,bcr,testtype,flux2type,t);
  for (i = 0; i < Elements.size(); i++){
    for (j = 0; j < element::systemsize2; j++){
      Elements[i].Q[j]=Elements[i].RHSofQ[j]*(1.0/Elements[i].length);
      Elements[i].RHSofQ[j].clear();
    }
  }
  //U* = U1+dt*R(U1,Q1)/(1 + dt*D)
  //Update U* variables
  fluxproc(Nodes,Elements,bcl,bcr,testtype,fluxtype,t);

  //Determine dt. 
//   for (i = 0; i < Elements.size(); i++){
//     dttemp = element::CFL*Elements[i].length/Elements[i].ws;
//     if (dttemp < dt) {
//       dt = dttemp;
//     }
//     dttemp = element::CFL*Elements[i].length*Elements[i].length/Elements[i].dds;
//     if (dttemp < dt) {
//       dt = dttemp;
//     }  
//     //    if (dt>0.0001){
//       //      dt = 0.0001;   
//     //    } 
//   }  

  sourceproc(Nodes,Elements,sourcetype);
  for (i = 0; i < Elements.size(); i++){
    Elements[i].StabOp_old=Elements[i].StabOp;
    for (j = 0; j < element::systemsize1; j++){
      Elements[i].U_old[j]=Elements[i].U[j];
      Elements[i].U[j]+=Elements[i].RHSofU[j]*(dt/Elements[i].length);
      //Stabilisation Operator:
      for (k = 1; k < state::statesize; k++){
	Elements[i].U[j].coefficient[k]=Elements[i].U[j].coefficient[k]/(1.0+dt*Elements[i].StabOp.coefficient[k]);
      }
      //Elements[i].RHSofU[j].clear();
    }
    Elements[i].StabOp.clear();
  }
  //Q* = R_q(U*) 
  //Update Q* variables
  flux2proc(Nodes,Elements,bcl,bcr,testtype,flux2type,t); //+dt?
  for (i = 0; i < Elements.size(); i++){
    for (j = 0; j < element::systemsize2; j++){
      Elements[i].Q[j]=Elements[i].RHSofQ[j]*(1.0/Elements[i].length);
      Elements[i].RHSofQ[j].clear();
    }
  }
  //U2 = U1(1 - dt/2*D*) +dt/2*(R(U*,Q*) + R(U1,Q1))/(1 - dt/2*D)
  //Update U* variables
  fluxproc(Nodes,Elements,bcl,bcr,testtype,fluxtype,t); //+dt?
  sourceproc(Nodes,Elements,sourcetype);
  for (i = 0; i < Elements.size(); i++){
    for (j = 0; j < element::systemsize1; j++){
      Elements[i].U[j]=Elements[i].U_old[j];
      //Stabilisation Operator:
      for (k = 1; k < state::statesize; k++){
	Elements[i].U[j].coefficient[k]=Elements[i].U[j].coefficient[k]*(1.0-0.5*dt*Elements[i].StabOp.coefficient[k]);
      }
      Elements[i].U[j]+=Elements[i].RHSofU[j]*(0.5*dt/Elements[i].length);
      // Elements[i].U[j]+=Elements[i].RHSofU_old[j]*(0.5*dt/Elements[i].length);
      for (k = 1; k < state::statesize; k++){
	//Stabilisation Operator:
	Elements[i].U[j].coefficient[k]=Elements[i].U[j].coefficient[k]/(1.0+0.5*dt*Elements[i].StabOp.coefficient[k]);
      }
      Elements[i].RHSofU[j].clear();
      Elements[i].U_old[j].clear();
    }
    Elements[i].StabOp.clear();
    Elements[i].StabOp_old.clear();
  }
}
    
  
