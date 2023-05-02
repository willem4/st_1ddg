#include <vector>
#include <iostream>

#include "node.h"
#include "state.h"
#include "element.h"
#include "flux.h"
#include "source.h"
#include "flux2.h"

#include "fluxproc.h"
#include "flux2proc.h"
#include "sourceproc.h"

#include "rk3.h"

using namespace std;

int rk3(double &dt, vector<node> &Nodes, vector<element>& Elements, int bcl, int bcr, testtypedescr testtype, fluxbase * fluxtype, sourcebase * sourcetype, flux2base * flux2type, double t){ 
  int i,j,k;
  double dttemp;
  for (i = 0; i < Elements.size(); i++){
    Elements[i].dds=0;
    Elements[i].ws=0;
  }
  //  Q1 = R_q(U1)
  flux2proc(Nodes,Elements,bcl,bcr,testtype,flux2type,t);
  for (i = 0; i < Elements.size(); i++){
    for (j = 0; j < element::systemsize2; j++){
      Elements[i].Q[j]=Elements[i].RHSofQ[j]*(1.0/Elements[i].length);
      Elements[i].RHSofQ[j].clear();
    }
  }  
  //  U2 = U1+dt*R(U1,Q1) 
  fluxproc(Nodes,Elements,bcl,bcr,testtype,fluxtype,t);
  
  //Determine dt. 
  for (i = 0; i < Elements.size(); i++){
    //    dttemp = element::CFL/Elements[i].ws;
    //    //    dttemp = element::CFL/(Elements[i].ws*Elements[i].ws/Elements[i].length + Elements[i].ws/Elements[i].length/Elements[i].length);
    //    if (dttemp < dt) {
    //      dt = dttemp;
    //    }
    //    dttemp = element::CFL*Elements[i].length*Elements[i].length/Elements[i].dds;
    //    if (dttemp < dt) {
    //      dt = dttemp;
    //    }  
  }
  //  if (dt > 0.0005){
  //    dt = 0.0005;
  //  }

  sourceproc(Nodes,Elements,sourcetype);
  for (i = 0; i < Elements.size(); i++){
    for (j = 0; j < element::systemsize1; j++){
      Elements[i].U_old[j]=Elements[i].U[j];
      Elements[i].U[j]+=Elements[i].RHSofU[j]*(dt/Elements[i].length);
      //Stabilisation Operator:
      for (k = 1; k < state::statesize; k++){
	Elements[i].U[j].coefficient[k]=Elements[i].U[j].coefficient[k]/(1.0+dt*Elements[i].StabOp.coefficient[k]);
      }
      Elements[i].RHSofU[j].clear();
      Elements[i].StabOp.clear();
    }
  }
  //  Q2 = R_q(U2)
  flux2proc(Nodes,Elements,bcl,bcr,testtype,flux2type,t);
  for (i = 0; i < Elements.size(); i++){
    for (j = 0; j < element::systemsize2; j++){
      Elements[i].Q[j]=Elements[i].RHSofQ[j]*(1.0/Elements[i].length);
      Elements[i].RHSofQ[j].clear();
    }
  }
  //  U3 = 3/4*U1+1/4*U2+1/4*dt*R(U2,Q2) 
  fluxproc(Nodes,Elements,bcl,bcr,testtype,fluxtype,t);
  sourceproc(Nodes,Elements,sourcetype);
  for (i = 0; i < Elements.size(); i++){
    for (j = 0; j < element::systemsize1; j++){
      Elements[i].U[j]+=Elements[i].RHSofU[j]*(dt/Elements[i].length);
      Elements[i].U[j]=Elements[i].U[j]*(1.0/4.0);
      //Stabilisation Operator:
      for (k = 1; k < state::statesize; k++){
	Elements[i].U[j].coefficient[k]=Elements[i].U[j].coefficient[k]/(1.0+dt*Elements[i].StabOp.coefficient[k]);
      }
      Elements[i].U[j]+=Elements[i].U_old[j]*(3.0/4.0);
      Elements[i].RHSofU[j].clear();
      Elements[i].StabOp.clear();
    }
  }
  //  Q3 = R_q(U3)
  flux2proc(Nodes,Elements,bcl,bcr,testtype,flux2type,t);
  for (i = 0; i < Elements.size(); i++){
    for (j = 0; j < element::systemsize2; j++){
      Elements[i].Q[j]=Elements[i].RHSofQ[j]*(1.0/Elements[i].length);
      Elements[i].RHSofQ[j].clear();
    }
  }  
  //U4 = 1/3*U1+2/3*U3+2/3*dt*R(U3,Q3)
  fluxproc(Nodes,Elements,bcl,bcr,testtype,fluxtype,t);
  sourceproc(Nodes,Elements,sourcetype);
  for (i = 0; i < Elements.size(); i++){
    for (j = 0; j < element::systemsize1; j++){
      Elements[i].U[j]+=Elements[i].RHSofU[j]*(dt/Elements[i].length);
      Elements[i].U[j]=Elements[i].U[j]*(2.0/3.0);
      //Stabilisation Operator:
      for (k = 1; k < state::statesize; k++){
	Elements[i].U[j].coefficient[k]=Elements[i].U[j].coefficient[k]/(1.0+dt*Elements[i].StabOp.coefficient[k]);
      }
      Elements[i].U[j]+=Elements[i].U_old[j]*(1.0/3.0);
      Elements[i].RHSofU[j].clear();
      Elements[i].U_old[j].clear();
      Elements[i].StabOp.clear();
    }
  }
  //   cout << Elements[2].U1[0].coefficient[0] << ": check after flux\n";
  return 0;
}

