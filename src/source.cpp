#include <iostream>
#include <cmath>

#include "physpara.h"
#include "element.h"
#include "source.h"

using namespace std;

void topography::eval(double * F2, double * uL, double * uLslope, double * qL){ 
  int j;  
  double hL,mL,uuL,hbL,qqL,hLslope,mLslope,hbLslope;
  double g = gravity_const;  /// * H_0 / (U_0)^2;
  hL=uL[0];
  mL=uL[1];
  uuL=mL/hL;
  hbL= uL[2];    
  hLslope = uLslope[0];
  mLslope = uLslope[1];
  hbLslope = uLslope[2];
  for (j = 0; j < element::systemsize1; j++){
    F2[j] = 0.0;
  }
  F2[1] = -g*hL*hbLslope;
}  

void exchangesusp::eval(double * F2, double * uL, double * uLslope, double * qL){ 
  int j;  
  double hL,mL,uuL,hbL,qqL,hLslope,mLslope,hbLslope;
  double g = gravity_const;  /// * H_0 / (U_0)^2;
  double A = Grass_const;
  hL=uL[0];
  mL=uL[1];
  uuL=mL/hL;
  hbL= uL[2];    
  double chihL= uL[3];
  hLslope = uLslope[0];
  mLslope = uLslope[1];
  hbLslope = uLslope[2];
  //  cout << chihL << " " << uuL*uuL << "\n";
  double D = 0.15*chihL;
  double E = 0.05*3*A*abs(uuL*uuL*uuL);
  for (j = 0; j < element::systemsize1; j++){
    F2[j] = 0.0;
  }
  F2[0] = 0;
  F2[1] = -g*hL*hbLslope;
  F2[2] = D-E;
  F2[3] = E-D;
  //  cout << 0.1*(D-E)/hL << "," << chiL << "\t";
}  
