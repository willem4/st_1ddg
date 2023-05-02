#include <vector>
#include <iostream>
#include <math.h>

#include "state.h"
#include "element.h"
#include "node.h"

#include "source.h"
#include "sourceproc.h"

using namespace std;

int sourceproc(vector<node> Nodes, vector<element>& Elements, sourcebase * sourcetype)
{
  int i,j;

  double * F1 = new double [element::systemsize1];
  double * F2 = new double [element::systemsize1];
  for (j = 0; j < element::systemsize1; j++){
    F1[j] = 0;
    F2[j] = 0;
  }

  for (i = 0; i < Elements.size(); i++){
    sourcetype->eval(F1,Elements[i].u1,Elements[i].u1slope,Elements[i].q1);
    sourcetype->eval(F2,Elements[i].u2,Elements[i].u2slope,Elements[i].q2);
    for (j = 0; j < element::systemsize1; j++){
      Elements[i].RHSofU[j].coefficient[0]+=0.5*(F1[j]+F2[j])*Elements[i].length;
      if (state::statesize>1){
	Elements[i].RHSofU[j].coefficient[1]+=1.5*(F1[j]*(-element::ep)+F2[j]*element::ep)*Elements[i].length;
      }
    }
  }
  
  delete [] F1;
  delete [] F2;
}

