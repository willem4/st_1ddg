#include <iostream>
#include <cmath>

#include "element.h"
#include "physpara.h"
#include "flux2.h"

using namespace std;

double diffusion::eval(double* F, double* uL, double* uR){
  //LDG diffusion
  int j;
  double ds;
  for (j=0;j<element::systemsize2;j++){
    F[j] = uR[j];
    ds = abs(uR[j]);
    if (abs(uR[j])>abs(uL[j])) {
      ds = abs(uR[j]);
    }
    else {
      ds = abs(uL[j]);
    }
  }
  return ds;
}

double sedimentdiffusion::eval(double* F, double* uL, double* uR){
  //LDG diffusion
  int j;
  double ds;
  for (j=0;j<element::systemsize2;j++){
    F[j] = uR[2];
    ds = abs(uR[2]);
  }
  return ds;
}
