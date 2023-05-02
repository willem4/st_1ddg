#include "testfcn.h"
#include "readtest.h"
#include "physpara.h"
#include "enumtypes.h"
#include "exactriemann.h"

#include <math.h>
#include <iostream>




double testfcn(testtypedescr testtype, int j, double x)
{
  double t = 0.0;
  return testfcn(testtype,j,x,t);
}

double testfcn(testtypedescr testtype, int j, double x, double t)
{
  double f;
  double pi = 4.0*atan(1.0);
  double K = 2.0;
  double g = gravity_const;
  double delta, xp, hc;
  double u;
  double L, nu, alpha;
  switch (testtype) {
  case Advection:
    //Advection.
    switch (j){
    case 0:
      f=x*(x-3);
      f = f/2.25;
      f=-f*f*f;
      break;
    }
    break;
  case Burgers:
    //Burgers.
    switch (j){
    case 0:
      f = sin(2*pi*x);
      break;
    }
    break;
  case SWEBurgers: 
    //Shallow water equations - Burgers test. 
    switch (j){
    case 0:
      f = 0.4*sin(2*pi*x);
      f = (K - f)/(3*sqrt(g));
      f = f*f;
      break;
    case 1:
      f = testfcn(testtype, 0, x);
      f = (K - 2*sqrt(g)*sqrt(f))*f;
      break;
    }
    break;
  case SWELinearWaveSolution: 
    //Shallow water equations - Linear wave solution.
    switch (j){
    case 0:
      f = (50.0-x);
      f = -f*f;
      f = 0.01*f;
      f = 0.001*exp(f)+10;
      break;
    case 1:
      f = 0;
      break;
    }
    break;
  case SWERiemannProblemLeftRarefactionRightShock:
    f = exactriemann(2.0, 0.0, 1.0, 0.0, 0.5, x, t, j);
    break;
  case SWERiemannProblemLeftShockRightRarefaction: 
    f = exactriemann(1, 0, 2, -1, 0.5, x, t, j);
    break;
  case SWERiemannProblemLeftShockRightShock: 
    f = exactriemann(1, 2, 1, -2, 0.5, x, t, j);
    break;
  case SWERiemannProblemLeftRarefactionRightRarefaction:
    f = exactriemann(2, -1, 2, 1, 0.5, x, t, j);
    break;
  case SWEContinuousTopography:   case SWEDiscontinuousTopography: 
    //Rest = Rest case.
    switch (j){
    case 0:
      f = 10 - testfcn(testtype,2,x);
      break;
    case 1:
      f = 0;
      break;
    case 2: 
      if (x>=40 & x <=60){
	f = sin(5*pi*x*0.01);
      }
      else{
	f = 0;
      }
      break;
    }
    break;
  case SWEFlowOverIsolatedContinuousRidgeI: 
  case SWEFlowOverIsolatedDiscontinuousRidgeI: 
      xp = 10;
      delta = 0.05;
      hc=0.2;
    switch (j){
    case 0:
      f = 0.2 - testfcn(testtype,2,x);
      break;
    case 1:
      //      if ((x - xp)>=-sqrt(hc/delta) & (x - xp)<= 0){
      //	u = 1.0286 - (xp - x)/sqrt(hc/delta)*(1.0286-0.3962);
      //      }
      //      if ((x - xp)<=sqrt(hc/delta) & (x - xp) > 0) {
      // 	u = 1.0286 - (x - xp)/sqrt(hc/delta)*(1.0286-0.3962);
      //      }
      //      if ((x-xp) > sqrt(hc/delta) or (x-xp) < -sqrt(hc/delta)){
	u = 0.3962;
	//      }
      f = u * testfcn(testtype,0,x);
      //      cout << f << "\t";
      break;
    case 2: 
      if ((x - xp)>=-sqrt(hc/delta) & (x - xp)<= sqrt(hc/delta)){
	f = x - xp;
	f= f*f;
	f = -delta*f;
      }
      else{
	f = -hc;
      }
      break;
    }
    break;
  case SWEFlowOverSinyBedI: 
      xp = 10;
      delta = 0.05;
      hc=0.2;
    switch (j){
    case 0:
      f = 0.2 - testfcn(testtype,2,x);
      break;
    case 1:
      //      if ((x - xp)>=-sqrt(hc/delta) & (x - xp)<= 0){
      //	u = 1.0286 - (xp - x)/sqrt(hc/delta)*(1.0286-0.3962);
      //      }
      //      if ((x - xp)<=sqrt(hc/delta) & (x - xp) > 0) {
      // 	u = 1.0286 - (x - xp)/sqrt(hc/delta)*(1.0286-0.3962);
      //      }
      //      if ((x-xp) > sqrt(hc/delta) or (x-xp) < -sqrt(hc/delta)){
	u = 0.3962;
	//      }
      f = u * testfcn(testtype,0,x);
      //      cout << f << "\t";
      break;
    case 2: 
      if (x>=3 & x <=17){
	f = 0.05*sin(pi*x);
      }
      else{
	f = 0;
      }
      break;
    }
    break;
  case SWEFlowOverIsolatedContinuousRidgeIV:
  case SWEFlowOverIsolatedDiscontinuousRidgeIV:
    u = 3.7637;
    switch (j){
    case 0:
      f = 0.2 - testfcn(testtype,2,x);
      break;
    case 1:
      f = u * testfcn(testtype,0,x);
      break;
    case 2: 
       xp = 10;
       delta = 0.05;
       hc=0.2;
       if ((x - xp)>=-sqrt(hc/delta) & (x - xp)<= sqrt(hc/delta)){
	 f = x - xp;
	 f= f*f;
	 f = -delta*f;
       }
       else{
	 f = -hc;
       }
      break;
    }
    break;
  case BurgersDiffusive:
    switch (j){
    case 0:
      L = 4; 
      nu = 1;
      alpha = 4;
      f = -2*(alpha*nu/L)*(sinh(alpha*x/L)/(cosh(alpha*x/L)+exp(-(alpha*alpha*nu*t)/(L*L))));
      break;    
    }
    break;
    //case default:
    //    cout << "Unknown initial conditions. \n"; 
    //    f = -11;
    //    break;
  case SedimentTransport:    
  case LFGrassMomentum:
  case LFGrassMomentumDiffusion:
    //Hudson Test
//     switch (j){
//     case 0:
//       f = 10 - testfcn(testtype,2,x);
//       break;
//     case 1:
//       f = 10;
//       break;
//     case 2: 
//       if (x > 300 and x < 500){
// 	f = sin(sin(pi*(x-300)/200));
// 	  }
//       else{
// 	f = 0;
//       }
//       break;
//     }
    xp = 10;
    delta = 0.05;
    hc=0.2;
    switch (j){
    case 0:
      f = 1 - testfcn(testtype,2,x);//+sin(pi*x/5)*0.05;
      break;
    case 1:
      //      if ((x - xp)>=-sqrt(hc/delta) & (x - xp)<= 0){
      //	u = 1.0286 - (xp - x)/sqrt(hc/delta)*(1.0286-0.3962);
      //      }
      //      if ((x - xp)<=sqrt(hc/delta) & (x - xp) > 0) {
      // 	u = 1.0286 - (x - xp)/sqrt(hc/delta)*(1.0286-0.3962);
      //      }
      //      if ((x-xp) > sqrt(hc/delta) or (x-xp) < -sqrt(hc/delta)){
      //      u = 0.3962;
      //      }
      //      f = u * testfcn(testtype,0,x);
      f = 0.4754;
      //      cout << f << "\t";
      break;
    case 2: 
      //       if (x<xp-1){
      // 	f = -0.25;
      //       }
      //       if (x>=(xp-1) and x <= xp){
      // 	f = (x-xp)*(0.25-0.12)-0.12;
      //       }
      //       if (x>xp){
      // 	f = -0.12;
      //       }
      if ((x - xp)>=-sqrt(hc/delta) & (x - xp)<= sqrt(hc/delta)){
  	f = x - xp;
  	f= f*f;
  	f = -0.1*delta*f - 0.9*hc;
      }
      else{
	 f = -hc;
      }
      break;
    }
    break;
  case HLLGrass:
  case LFGrass:
    xp = 10;
    delta = 0.05;
    hc=0.2;
    switch (j){
    case 0:
      f = 1 - testfcn(testtype,2,x);//+sin(pi*x/5)*0.05;
      break;
    case 1:
      f = 0.4754;
      break;
    case 2: 
//       if (x<xp-1){
// 	f = -0.25;
//       }
//       if (x>=(xp-1) and x <= xp){
// 	f = (x-xp)*(0.25-0.12)-0.12;
//       }
//       if (x>xp){
// 	f = -0.12;
//       }
      if ((x - xp)>=-sqrt(hc/delta) & (x - xp)<= sqrt(hc/delta)){
 	f = x - xp;
 	f= f*f;
 	f = -0.1*delta*f - 0.9*hc;
      }
      else{
 	f = -hc;
      }
      break;
    }
    break;
  case DecoupledGrassSedimentOnly:
    xp = 10;
    delta = 0.05;
    hc=0.2;
    switch (j){
    case 0: 
      f = 0.4754/(1-testfcn(testtype,1,x));
      break;
    case 1:
      if ((x - xp)>=-sqrt(hc/delta) & (x - xp)<= sqrt(hc/delta)){
 	f = x - xp;
 	f= f*f;
 	f = -0.1*delta*f - 0.9*hc;
      }
      else{
 	f = -hc;
      }
      break;
    }
    break;
  case SuspendedSediment:
    xp = 10;
    delta = 0.05;
    hc=0.2;
    switch (j){
    case 0: 
      f = 1-testfcn(testtype,2,x);
      break;
    case 1:
      f = 0.4754;
      break;
    case 2:
      if ((x - xp)>=-sqrt(hc/delta) & (x - xp)<= sqrt(hc/delta)){
  	f = x - xp;
  	f= f*f;
  	f = -1.0*delta*f - 0.0*hc;
      }
      else{
 	f = -hc;
      }
      //      f=2*delta*sin(pi*x/5);
      break;
    case 3:
      //      if (x > xp){
      f = 0.05;//*testfcn(testtype,0,x);
      //      }
      //      else {
      //      f = 0;
      //      }
      break;
    }
    break;
  }
  return f;

}
