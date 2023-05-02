#include <iostream>
#include <cmath>

#include "physpara.h"
#include "element.h"
#include "flux.h"

using namespace std;

double advection::eval(double* F, double* uL, double* uR, double* qL, double* qR){
  //Upwind-flux for advection 
  int j;
  double a = advection_speed;
  for (j=0;j<element::systemsize1;j++){
    if (a>=0){
      F[j]=a*uL[j];
    }
    else {
      F[j]=a*uR[j];
    }
  }
  return abs(a);
}

double burgers::eval(double* F, double* uL, double* uR, double* qL, double* qR){
  //Enquist Osher Flux for Burger's equation.
  int j;
  double ws;
  for (j=0;j<element::systemsize1;j++){
    if ((uL[j]<0) & (uR[j]<0)){
      F[j]=0.5*uR[j]*uR[j];
      ws = abs(uR[j]);
    }
    if ((uL[j]>0) & (uR[j]>0)){
      F[j]=0.5*uL[j]*uL[j];
      ws = abs(uL[j]);
    }
    if ((uL[j]<0) & (uR[j]>0)){
      F[j]=0;
      ws = 0;
    }
    if ((uL[j]>0) & (uR[j]<0)){
      F[j]=0.25*(uL[j]*uL[j]+uR[j]*uR[j]);
      ws = abs(0.5*(uL[j]+uR[j]));
    }
  }
  return ws;
}

double burgersdiffusion::eval(double* F, double* uL, double* uR, double* qL, double* qR){
  //Enquist Osher Flux for Burger's equation.
  int j;
  double ws;
  for (j=0;j<element::systemsize1;j++){
    if ((uL[j]<0) & (uR[j]<0)){
      F[j]=0.5*uR[j]*uR[j];
      ws = abs(uR[j]);
    }
    if ((uL[j]>0) & (uR[j]>0)){
      F[j]=0.5*uL[j]*uL[j];
      ws = abs(uL[j]);
    }
    if ((uL[j]<0) & (uR[j]>0)){
      F[j]=0;
    }
    if ((uL[j]>0) & (uR[j]<0)){
      F[j]=0.25*(uL[j]*uL[j]+uR[j]*uR[j]);
      ws = abs(0.5*(uL[j]+uR[j]));
    }
    if (abs(uR[j])>abs(uL[j])) {
      ws = abs(uR[j]);
    }
    else {
      ws = abs(uL[j]);
    }
    F[j]+=qL[j];
  }
  return ws;
}

double swehllc::eval(double* F, double* uL, double* uR, double* qL, double* qR){
  //HLLC Flux for the shallow water equations 
  // including upwind flux for Grass bed updating equation. 
  int j; 
  double hL,mL,hR,mR,uuL,uuR,SL,SR,SM,hM,mM,hbL,hbR,uu;
  double g = gravity_const; // scaled = H_0 / (U_0)^2 g;
  for (j = 0; j < element::systemsize1; j++){  
    F[j]=0.0;
  }
  hL=uL[0];   
  mL=uL[1];   
  uuL=mL/hL;
  hR=uR[0];
  mR=uR[1];
  uuR=mR/hR;
  if (uuL-( sqrt(g*hL))<uuR-(sqrt(g*hR)))
    {
      SL=uuL-sqrt(g*hL);
    }
  else 
    {
      SL=uuR-sqrt(g*hR);
    }
  if (uuL+sqrt(g*hL)>uuR+sqrt(g*hR))
    {
      SR=uuL+sqrt(g*hL);
    }
  else 
    {
      SR=uuR+sqrt(g*hR);
    }   
  double ws = abs(SL); 
  if (abs(SR) > ws){
    ws = abs(SR);
  }
  SM=(0.5*g*hL*hL-0.5*g*hR*hR-uuL*hL*(SL-uuL)+uuR*hR*(SR-uuR))/(hR*(SR-uuR)-hL*(SL-uuL));
  if (SL > 0){
    F[0] = mL;
    F[1] = mL*mL/hL+0.5*g*hL*hL;
  }
  if (SR < 0){
      F[0] = mR;
      F[1] = mR*mR/hR+0.5*g*hR*hR;
    }
  if (SL <= 0 and SM > 0){
      hM = hL*(SL-uuL)/(SL-SM);         //sqrt(2/g*(0.5*g*hL*hL+hL*(SL-uuL)*(SM-uuL)));
      mM = SM*hM;
      F[0] = mL+SL*(hM-hL);
      F[1] = mL*mL/hL+0.5*g*hL*hL+SL*(mM-mL);
    }
  if (SM <= 0 and SR >= 0){
      hM = hR*(SR-uuR)/(SR-SM); //sqrt(2.0/g*(0.5*g*hR*hR+hR*(SR-uuR)*(SM-uuR)));
      mM = SM*hM;
      F[0] = mR+SR*(hM-hR);
      F[1] = mR*mR/hR+0.5*g*hR*hR+SR*(mM-mR);
    }
  for (j = 2; j < element::systemsize1; j++){
    F[j]=0;
  }
  return ws;
}  

double swehllcwidth::eval(double* F, double* uL, double* uR, double* qL, double* qR){
  //HLLC Flux for the shallow water equations 
  // including upwind flux for Grass bed updating equation. 
  int j; 
  double hL,mL,hR,mR,uuL,uuR,SL,SR,SM,hM,mM,hbL,hbR,uu,wL,wR,wM;
  double g = gravity_const; // scaled = H_0 / (U_0)^2 g;
  for (j = 0; j < element::systemsize1; j++){  
    F[j]=0.0;
  }
  hL=uL[0];   
  mL=uL[1];   
  uuL=mL/hL;
  hR=uR[0];
  mR=uR[1];
  uuR=mR/hR;
  if (uuL-( sqrt(g*hL))<uuR-(sqrt(g*hR)))
    {
      SL=uuL-sqrt(g*hL);
    }
  else 
    {
      SL=uuR-sqrt(g*hR);
    }
  if (uuL+sqrt(g*hL)>uuR+sqrt(g*hR))
    {
      SR=uuL+sqrt(g*hL);
    }
  else 
    {
      SR=uuR+sqrt(g*hR);
    }   
  double ws = abs(SL); 
  if (abs(SR) > ws){
    ws = abs(SR);
  }
  SM=(0.5*g*hL*hL-0.5*g*hR*hR-uuL*hL*(SL-uuL)+uuR*hR*(SR-uuR))/(hR*(SR-uuR)-hL*(SL-uuL));
  if (SL > 0){
    F[0] = mL;
    F[1] = mL*mL/hL+0.5*g*hL*hL;
  }
  if (SR < 0){
      F[0] = mR;
      F[1] = mR*mR/hR+0.5*g*hR*hR;
    }
  if (SL <= 0 and SM > 0){
      hM = hL*(SL-uuL)/(SL-SM);         //sqrt(2/g*(0.5*g*hL*hL+hL*(SL-uuL)*(SM-uuL)));
      mM = SM*hM;
      F[0] = mL+SL*(hM-hL);
      F[1] = mL*mL/hL+0.5*g*hL*hL+SL*(mM-mL);
    }
  if (SM <= 0 and SR >= 0){
      hM = hR*(SR-uuR)/(SR-SM); //sqrt(2.0/g*(0.5*g*hR*hR+hR*(SR-uuR)*(SM-uuR)));
      mM = SM*hM;
      F[0] = mR+SR*(hM-hR);
      F[1] = mR*mR/hR+0.5*g*hR*hR+SR*(mM-mR);
    }
  for (j = 2; j < element::systemsize1; j++){
    F[j]=0;
  }
  return ws;
}

double swehllctopography::eval(double* F, double* uL, double* uR, double* qL, double* qR){
  //HLLC Flux for the shallow water equations 
  int j; 
  double hL,mL,hR,mR,uuL,uuR,SL,SR,SM,hM,mM,hbL,hbR,uu;
  double g = gravity_const; // scaled = H_0 / (U_0)^2 g;
  for (j = 0; j < element::systemsize1; j++){  
    F[j]=0.0;
  }
  hL=uL[0];   
  mL=uL[1];   
  uuL=mL/hL;
  hR=uR[0];
  mR=uR[1];
  uuR=mR/hR;
  hbL=uL[2];  
  hbR=uR[2];
  //Solution for discontinuous bottom.
  if (hbL>hbR){
    hR -= (hbL-hbR);
    mR = hR*uuR;
  }
  else{
    hL -= (hbR-hbL); 
    mL = hL*uuL;
  }
  if (uuL-( sqrt(g*hL))<uuR-(sqrt(g*hR)))
    {
      SL=uuL-sqrt(g*hL);
    }
  else 
    {
      SL=uuR-sqrt(g*hR);
    }
  if (uuL+sqrt(g*hL)>uuR+sqrt(g*hR))
    {
      SR=uuL+sqrt(g*hL);
    }
  else 
    {
      SR=uuR+sqrt(g*hR);
    }   
  double ws = abs(SL); 
  if (abs(SR) > ws){
    ws = abs(SR);
  }
  SM=(0.5*g*hL*hL-0.5*g*hR*hR-uuL*hL*(SL-uuL)+uuR*hR*(SR-uuR))/(hR*(SR-uuR)-hL*(SL-uuL));
  if (SL > 0){
      F[0] = mL;
      F[1] = mL*mL/hL+0.5*g*hL*hL;
    }
  if (SR < 0){
      F[0] = mR;
      F[1] = mR*mR/hR+0.5*g*hR*hR;
    }
  if (SL <= 0 and SM > 0){
      hM = hL*(SL-uuL)/(SL-SM);  //sqrt(2/g*(0.5*g*hL*hL+hL*(SL-uuL)*(SM-uuL)));
      mM = SM*hM;
      F[0] = mL+SL*(hM-hL);
      F[1] = mL*mL/hL+0.5*g*hL*hL+SL*(mM-mL);
    }
  if (SM <= 0 and SR >= 0){
      hM = hL*(SL-uuL)/(SL-SM);  //sqrt(2.0/g*(0.5*g*hR*hR+hR*(SR-uuR)*(SM-uuR)));
      mM = SM*hM;
      F[0] = mR+SR*(hM-hR);
      F[1] = mR*mR/hR+0.5*g*hR*hR+SR*(mM-mR);
    }
  for (j = 2; j < element::systemsize1; j++){
    F[j]=0;
  }
  return ws;
}  

void swehllctopography::eval_l(double* Fl, double* uL, double* uR, double* qL, double* qR){
  int j; 
  double hL,mL,hR,mR,uuL,uuR,SL,SR,SM,hM,mM,hbL,hbR,h,hb;
  double g = gravity_const; // scaled = H_0 / (U_0)^2 g;
  for (j = 0; j < element::systemsize1; j++){  
    Fl[j]=0;
  }
  hL=uL[0];
  mL=uL[1];
  uuL=mL/hL;
  hR=uR[0];
  mR=uR[1];
  uuR=mR/hR;
  hbL=uL[2];  
  hbR=uR[2];
  if (hbL>hbR){
    h = hR - (hbL-hbR);
    hb = hbL - hbR;
  }
  else {
    h = hR;
    hb = 0;
  }
  Fl[1]+=0.5*(g*hR*hR - g*h*h); //+ g*hR*hbR; 
}

void swehllctopography::eval_r(double* Fr, double * uL, double * uR, double* qL, double* qR){
  int j;
  double hL,mL,hR,mR,uuL,uuR,SL,SR,SM,hM,mM,hbL,hbR,h,hb;
  double g = gravity_const; // scaled = H_0 / (U_0)^2 g;
  for (j = 0; j < element::systemsize1; j++){  
    Fr[j]=0;
  }
  hL=uL[0];
  mL=uL[1];
  uuL=mL/hL;
  hR=uR[0];
  mR=uR[1];
  uuR=mR/hR;
  hbL=uL[2];  
  hbR=uR[2];
  //Bouchut
  if (hbR>hbL){
    h = hL - (hbR-hbL);
    hb = (hbR-hbL);
  }
  else{
    h = hL; 
    hb = 0;
  }
  Fr[1]+=0.5*(g*hL*hL - g*h*h); // + g*hL*hbL; 
}

double sedimenttransport::eval(double* F, double* uL, double* uR, double* qL, double* qR){
  //HLLC Flux for the shallow water equations 
  // including upwind flux for Grass bed updating equation. 
  int j; 
  double hL,mL,hR,mR,uuL,uuR,SL,SR,SM,hM,mM,hbL,hbR,hbM,uu;
  double g = gravity_const; // scaled = H_0 / (U_0)^2 g;
  double A = Grass_const;
  double B = Grass_exponent;
  double ws = 0.0;
  for (j = 0; j < element::systemsize1; j++){  
    F[j]=0.0;
  }
  hL=uL[0];   
  mL=uL[1];   
  uuL=mL/hL;
  hbL=uL[2];  
  hR=uR[0];
  mR=uR[1];
  uuR=mR/hR;
  hbR=uR[2];

//Solution for discontinuous bottom.
   if (hbL>hbR){
     hR -= (hbL-hbR);
     mR = hR*uuR;
   }
   else{
     hL -= (hbR-hbL); 
     mL = hL*uuL;
   }

  double pi = 4.0*atan(1.0);
  int i;
  double dL;
  double dR;
  double aL[4];
  double aR[4];
  double eig[6];
  //characteristic polynomial used for calculating eigenvalues of left face.   
  dL=B*A*abs(pow(uuL,B-1));
  aL[0]=1;
  aL[1]=-2*uuL;
  aL[2]=-g*hL-g*dL+pow(uuL,2);
  aL[3]=g*uuL*dL;
  double QL = (3*aL[2]-pow(aL[1],2))/9;
  double RL = (9*aL[2]*aL[1]-27*aL[3]-2*pow(aL[1],3))/54;
  double theta = acos(RL/sqrt(pow(-QL,3)));
  eig[0] = 2*sqrt(-QL)*cos((theta)/3)-aL[1]/3;
  eig[1] = 2*sqrt(-QL)*cos((theta+2*pi)/3)-aL[1]/3;
  eig[2] = 2*sqrt(-QL)*cos((theta+4*pi)/3)-aL[1]/3;
  //characteristic polynomial used for calculating eigenvalues of right face.   
  dR=B*A*abs(pow(uuR,B-1));
  aR[0]=1;
  aR[1]=-2*uuR;
  aR[2]=-g*hR-g*dR+pow(uuR,2);
  aR[3]=g*uuR*dR;
  double QR = (3*aR[2]-pow(aR[1],2))/9;
  double RR = (9*aR[2]*aR[1]-27*aR[3]-2*pow(aR[1],3))/54;
  theta = acos(RR/sqrt(pow(-QR,3)));
  eig[3] = 2*sqrt(-QR)*cos((theta)/3)-aR[1]/3;
  eig[4] = 2*sqrt(-QR)*cos((theta+2*pi)/3)-aR[1]/3;
  eig[5] = 2*sqrt(-QR)*cos((theta+4*pi)/3)-aR[1]/3;
  SL = eig[0];
  SR = eig[0];
  for (i=1;i<6;i++){
    if (eig[i]<SL){
      SL = eig[i];
    }
    if (eig[i]>SR){
      SR = eig[i];
    }
    if (abs(eig[i])>ws){
      ws = abs(eig[i]);
    }
  }
  SM=(0.5*g*hL*hL-0.5*g*hR*hR-uuL*hL*(SL-uuL)+uuR*hR*(SR-uuR))/(hR*(SR-uuR)-hL*(SL-uuL));
  SM = (SR*uuR-SL*uuL+(uuL*uuL/2+g*(hL+hbL))-(uuR*uuR/2+g*(hR+hbR)))/(SR-SL);
  if (SL > 0){
    F[0] = mL;
    F[1] = mL*mL/hL+0.5*g*hL*hL;
    F[2] = A * uuL*abs(pow(uuL,B-1));
  }
  if (SR < 0){
    F[0] = mR;
    F[1] = mR*mR/hR+0.5*g*hR*hR;
    F[2] = A * uuR*abs(pow(uuR,B-1));
  }
  if (SL <= 0 and SM > 0){
    hM = hL*(SL-uuL)/(SL-SM);  //sqrt(2/g*(0.5*g*hL*hL+hL*(SL-uuL)*(SM-uuL)));
    mM = SM*hM;
    hbM = (SL*hbL - A*uuL*abs(pow(uuL,B-1)) + A*SM*abs(pow(SM,B-1)))/SL;
    //    hbM = (SR*hbR-SL*hbL+A*uuL*abs(pow(uuL,B-1)) - A*uuR*abs(pow(uuR,B-1)))/(SR-SL);
    //    cout << "hbM = " << hbM << "\n";
    //       if (hbL >= hbR){
    // 	hbM = hbR; 
    //       }
    //       else{
    // 	hbM = hbL;
    //       }
    F[0] = mL+SL*(hM-hL);
    F[1] = mL*mL/hL+0.5*g*hL*hL+SL*(mM-mL);
    F[2] = A*uuL*abs(pow(uuL,B-1))+SL*(hbM-hbL);
    //    F[2] = (SR*A*uuL*abs(pow(uuL,B-1))-SL*A*uuR*abs(pow(uuR,B-1))+SL*SR*(hbR-hbL))/(SR-SL);
    //      F[2] = A * pow((uuL+SM)/2,3); 
    //      F[2] = A * pow(SM,3); 
    //      F[2] = A * pow(mL/hL,3); 
    //      F[2] = A * pow(-SL/(SM-SL)*uuL + SM/(SM-SL)*SM,3);
  }
  if (SM <= 0 and SR >= 0){
    hM = hL*(SL-uuL)/(SL-SM);  //sqrt(2.0/g*(0.5*g*hR*hR+hR*(SR-uuR)*(SM-uuR)));
    mM = SM*hM;
    hbM = (SR*hbR - A*uuR*abs(pow(uuR,B-1)) + A*SM*abs(pow(SM,B-1)))/SR;
    //    hbM = (SR*hbR-SL*hbL+A*uuL*abs(pow(uuL,B-1)) - A*uuR*abs(pow(uuR,B-1)))/(SR-SL);
    //       if (hbL >= hbR){
    // 	hbM = hbR; 
    //       }
    //       else{
    // 	hbM = hbL;
    //       }
    F[0] = mR+SR*(hM-hR);
    F[1] = mR*mR/hR+0.5*g*hR*hR+SR*(mM-mR);
    F[2] = A*uuR*abs(pow(uuR,B-1))+SR*(hbM-hbR);
    //    F[2] = (SR*A*uuL*abs(pow(uuL,B-1))-SL*A*uuR*abs(pow(uuR,B-1))+SL*SR*(hbR-hbL))/(SR-SL);
    //      F[2] = A * pow((SM+uuR)/2,3); 
    //      F[2] = A * pow(SM,3); 
    //      F[2] = A * pow(mR/hR,3); 
    //      F[2] = A * pow(-SM/(SR-SM)*SM + SR/(SR-SM)*uuR,3);
  }
  //   if (uuL+uuR > 0) {
  //       F[2] = A * pow(uuL,1);
  //   }
  //   if (uuL+uuR < 0) {
  //     F[2] = A * pow(uuR,1);
  //   }
  for (j = 3; j < element::systemsize1; j++){
    F[j]=0;
  }
  return ws;
}  

void sedimenttransport::eval_l(double* Fl, double* uL, double* uR, double* qL, double* qR){
  int j; 
  double hL,mL,hR,mR,uuL,uuR,SL,SR,SM,hM,mM,hbL,hbR,h,hb;
  double g = gravity_const; // scaled = H_0 / (U_0)^2 g;
  for (j = 0; j < element::systemsize1; j++){  
    Fl[j]=0;
  }
  hL=uL[0];
  mL=uL[1];
  uuL=mL/hL;
  hR=uR[0];
  mR=uR[1];
  uuR=mR/hR;
  hbL=uL[2];  
  hbR=uR[2];
  if (hbL>hbR){
    h = hR - (hbL-hbR);
    hb = hbL - hbR;
  }
  else {
    h = hR;
    hb = 0;
  }
  Fl[1]+=0.5*(g*hR*hR - g*h*h); //+ g*hR*hbR; 
}

void sedimenttransport::eval_r(double* Fr, double * uL, double * uR, double* qL, double* qR){
  int j;
  double hL,mL,hR,mR,uuL,uuR,SL,SR,SM,hM,mM,hbL,hbR,h,hb;
  double g = gravity_const; // scaled = H_0 / (U_0)^2 g;
  for (j = 0; j < element::systemsize1; j++){  
    Fr[j]=0;
  }
  hL=uL[0];
  mL=uL[1];
  uuL=mL/hL;
  hR=uR[0];
  mR=uR[1];
  uuR=mR/hR;
  hbL=uL[2];  
  hbR=uR[2];
  //Bouchut
  if (hbR>hbL){
    h = hL - (hbR-hbL);
    hb = (hbR-hbL);
  }
  else{
    h = hL; 
    hb = 0;
  }
  Fr[1]+=0.5*(g*hL*hL - g*h*h); // + g*hL*hbL; 
}

double hllgrass::eval(double* F, double* uL, double* uR, double* qL, double* qR){
  //HLL flux for grass bed updating equation with SWE in not shock preserving form.
  double g = gravity_const; // scaled = H_0 / (U_0)^2 g;
  double pi = 4.0*atan(1.0);
  int i;
  double A = Grass_const;
  double B = Grass_exponent;
  double dL;
  double dR;
  double aL[4];
  double aR[4];
  double eig[6];
  double hL=uL[0];   
  double hR=uR[0];
  double uuL=uL[1];   
  double uuR=uR[1];
  double hbL=uL[2];  
  double hbR=uR[2];

  //characteristic polynomial used for calculating eigenvalues of left face.   
  dL=B*A*abs(pow(uuL,B-1));
  aL[0]=1;
  aL[1]=-2*uuL;
  aL[2]=-g*hL-g*dL+pow(uuL,2);
  aL[3]=g*uuL*dL;
  double QL = (3*aL[2]-pow(aL[1],2))/9;
  double RL = (9*aL[2]*aL[1]-27*aL[3]-2*pow(aL[1],3))/54;
  double theta = acos(RL/sqrt(pow(-QL,3)));
  eig[0] = 2*sqrt(-QL)*cos((theta)/3)-aL[1]/3;
  eig[1] = 2*sqrt(-QL)*cos((theta+2*pi)/3)-aL[1]/3;
  eig[2] = 2*sqrt(-QL)*cos((theta+4*pi)/3)-aL[1]/3;

  //characteristic polynomial used for calculating eigenvalues of right face.   
  dR=B*A*abs(pow(uuR,B-1));
  aR[0]=1;
  aR[1]=-2*uuR;
  aR[2]=-g*hR-g*dR+pow(uuR,2);
  aR[3]=g*uuR*dR;
  double QR = (3*aR[2]-pow(aR[1],2))/9;
  double RR = (9*aR[2]*aR[1]-27*aR[3]-2*pow(aR[1],3))/54;
  theta = acos(RR/sqrt(pow(-QR,3)));
  //  cout << QR << RR << theta << "\n";
  eig[3] = 2*sqrt(-QR)*cos((theta)/3)-aR[1]/3;
  eig[4] = 2*sqrt(-QR)*cos((theta+2*pi)/3)-aR[1]/3;
  eig[5] = 2*sqrt(-QR)*cos((theta+4*pi)/3)-aR[1]/3;
  double SL = eig[0];
  double SR = eig[0];
  double ws = abs(eig[0]);
  for (i=1;i<6;i++){
    if (eig[i]<SL){
      SL = eig[i];
    }
    if (eig[i]>SR){
      SR = eig[i];
    }
    if (abs(eig[i])>ws){
      ws = abs(eig[i]);
    }
  }
  if (SL>0){
    F[0]=hL*uuL;
    F[1]=pow(uuL,2)/2+g*(hL+hbL);
    F[2]=A*uuL*abs(pow(uuL,B-1));
  }
  if (SR<0){
    F[0]=hR*uuR;
    F[1]=pow(uuR,2)/2+g*(hR+hbR);
    F[2]=A*uuR*abs(pow(uuR,B-1));
  }
  if (SL<=0 and SR>=0){
    F[0]=(SR*hL*uuL-SL*hR*uuR+SL*SR*(hR-hL))/(SR-SL);
    F[1]=(SR*(pow(uuL,2)/2+g*(hL+hbL))-SL*(pow(uuR,2)/2+g*(hR+hbR))+SL*SR*(uuR-uuL))/(SR-SL);
    F[2]=(SR*A*uuL*abs(pow(uuL,B-1))-SL*A*uuR*abs(pow(uuR,B-1))+SL*SR*(hbR-hbL))/(SR-SL);
  }
  return ws;
}

double lfgrass::eval(double* F, double* uL, double* uR, double* qL, double* qR){
  //HLL flux for grass bed updating equation with SWE in not shock preserving form.
  double g = gravity_const; // scaled = H_0 / (U_0)^2 g;
  double pi = 4.0*atan(1.0);
  int i;
  double A = Grass_const;
  double B = Grass_exponent;
  double dL;
  double dR;
  double aL[4];
  double aR[4];
  double eig[6];
  double hL=uL[0];   
  double hR=uR[0];
  double uuL=uL[1];   
  double uuR=uR[1];
  double hbL=uL[2];  
  double hbR=uR[2];
  //characteristic polynomial used for calculating eigenvalues of left face.   
  dL=B*A*abs(pow(uuL,B-1));
  aL[0]=1;
  aL[1]=-2*uuL;
  aL[2]=-g*hL-g*dL+pow(uuL,2);
  aL[3]=g*uuL*dL;
  double QL = (3*aL[2]-pow(aL[1],2))/9;
  double RL = (9*aL[2]*aL[1]-27*aL[3]-2*pow(aL[1],3))/54;
  double theta = acos(RL/sqrt(pow(-QL,3)));
  eig[0] = 2*sqrt(-QL)*cos((theta)/3)-aL[1]/3;
  eig[1] = 2*sqrt(-QL)*cos((theta+2*pi)/3)-aL[1]/3;
  eig[2] = 2*sqrt(-QL)*cos((theta+4*pi)/3)-aL[1]/3;

  //characteristic polynomial used for calculating eigenvalues of right face.   
  dR=B*A*abs(pow(uuR,B-1));
  aR[0]=1;
  aR[1]=-2*uuR;
  aR[2]=-g*hR-g*dR+pow(uuR,2);
  aR[3]=g*uuR*dR;
  double QR = (3*aR[2]-pow(aR[1],2))/9;
  double RR = (9*aR[2]*aR[1]-27*aR[3]-2*pow(aR[1],3))/54;
  theta = acos(RR/sqrt(pow(-QR,3)));
  //  cout << QR << RR << theta << "\n";
  eig[3] = 2*sqrt(-QR)*cos((theta)/3)-aR[1]/3;
  eig[4] = 2*sqrt(-QR)*cos((theta+2*pi)/3)-aR[1]/3;
  eig[5] = 2*sqrt(-QR)*cos((theta+4*pi)/3)-aR[1]/3;
  double SR = eig[0];
  for (i=1;i<6;i++){
    if (abs(eig[i])>abs(SR)){
      SR = eig[i];
    }
  }
  F[0]=0.5*(hL*uuL+hR*uuR)-SR*(hR-hL)/2;
  F[1]=0.5*(pow(uuL,2)/2+g*(hL+hbL)+pow(uuR,2)/2+g*(hR+hbR))-SR*(uuR-uuL)/2;
  F[2]=0.5*(A*uuL*abs(pow(uuL,B-1))+A*uuR*abs(pow(uuR,B-1)))-SR*(hbR-hbL)/2;
  return abs(SR);
}

double lfgrassmomentum::eval(double* F, double* uL, double* uR, double* qL, double* qR){
  //HLL flux for grass bed updating equation with SWE in shock preserving form.
  double g = gravity_const; // scaled = H_0 / (U_0)^2 g;
  double pi = 4.0*atan(1.0);
  int i;
  double A = Grass_const;
  double B = Grass_exponent;
  double dL;
  double dR;
  double aL[4];
  double aR[4];
  double eig[6];
  double hL=uL[0];   
  double hR=uR[0];
  double mL=uL[1];   
  double mR=uR[1];
  double uuL = mL/hL;
  double uuR = mR/hR;
  double hbL=uL[2];  
  double hbR=uR[2];

  //Solution for discontinuous topography.
  if (hbL>hbR){
    hR -= (hbL-hbR);
    mR = hR*uuR;
  }
  else{
    hL -= (hbR-hbL); 
    mL = hL*uuL;
  }

  //characteristic polynomial used for calculating eigenvalues of left face.   
  dL=B*A*abs(pow(uuL,B-1));
  aL[0]=1;
  aL[1]=-2*uuL;
  aL[2]=-g*hL-g*dL+pow(uuL,2);
  aL[3]=g*uuL*dL;
  double QL = (3*aL[2]-pow(aL[1],2))/9;
  double RL = (9*aL[2]*aL[1]-27*aL[3]-2*pow(aL[1],3))/54;
  double theta = acos(RL/sqrt(pow(-QL,3)));
  eig[0] = 2*sqrt(-QL)*cos((theta)/3)-aL[1]/3;
  eig[1] = 2*sqrt(-QL)*cos((theta+2*pi)/3)-aL[1]/3;
  eig[2] = 2*sqrt(-QL)*cos((theta+4*pi)/3)-aL[1]/3;

  //characteristic polynomial used for calculating eigenvalues of right face.   
  dR=B*A*abs(pow(uuR,B-1));
  aR[0]=1;
  aR[1]=-2*uuR;
  aR[2]=-g*hR-g*dR+pow(uuR,2);
  aR[3]=g*uuR*dR;
  double QR = (3*aR[2]-pow(aR[1],2))/9;
  double RR = (9*aR[2]*aR[1]-27*aR[3]-2*pow(aR[1],3))/54;
  theta = acos(RR/sqrt(pow(-QR,3)));
  //  cout << QR << RR << theta << "\n";
  eig[3] = 2*sqrt(-QR)*cos((theta)/3)-aR[1]/3;
  eig[4] = 2*sqrt(-QR)*cos((theta+2*pi)/3)-aR[1]/3;
  eig[5] = 2*sqrt(-QR)*cos((theta+4*pi)/3)-aR[1]/3;
  double SR = eig[0];
  for (i=1;i<6;i++){
    if (abs(eig[i])>abs(SR)){
      SR = eig[i];
    }
  }
  F[0]=0.5*(hL*uuL+hR*uuR)-SR*(hR-hL)/2;
  F[1]=0.5*(hL*pow(uuL,2)+g*pow(hL,2)/2+hR*pow(uuR,2)+g*pow(hR,2)/2) -SR*(mR-mL)/2;
  F[2]=0.5*(A*uuL*abs(pow(uuL,B-1))+A*uuR*abs(pow(uuR,B-1)))-SR*(hbR-hbL)/2;
  return abs(SR);
}

void lfgrassmomentum::eval_l(double* Fl, double* uL, double* uR, double* qL, double* qR){
  int j; 
  double hL,mL,hR,mR,uuL,uuR,SL,SR,SM,hM,mM,hbL,hbR,h,hb;
  double g = gravity_const; // scaled = H_0 / (U_0)^2 g;
  for (j = 0; j < element::systemsize1; j++){  
    Fl[j]=0;
  }
  hL=uL[0];
  mL=uL[1];
  uuL=mL/hL;
  hR=uR[0];
  mR=uR[1];
  uuR=mR/hR;
  hbL=uL[2];  
  hbR=uR[2];
  if (hbL>hbR){
    h = hR - (hbL-hbR);
    hb = hbL - hbR;
  }
  else {
    h = hR;
    hb = 0;
  }
  Fl[1]+=0.5*(g*hR*hR - g*h*h); //+ g*hR*hbR; 
}

void lfgrassmomentum::eval_r(double* Fr, double * uL, double * uR, double* qL, double* qR){
  int j;
  double hL,mL,hR,mR,uuL,uuR,SL,SR,SM,hM,mM,hbL,hbR,h,hb;
  double g = gravity_const; // scaled = H_0 / (U_0)^2 g;
  for (j = 0; j < element::systemsize1; j++){  
    Fr[j]=0;
  }
  hL=uL[0];
  mL=uL[1];
  uuL=mL/hL;
  hR=uR[0];
  mR=uR[1];
  uuR=mR/hR;
  hbL=uL[2];  
  hbR=uR[2];
  //Bouchut
  if (hbR>hbL){
    h = hL - (hbR-hbL);
    hb = (hbR-hbL);
  }
  else{
    h = hL; 
    hb = 0;
  }
  Fr[1]+=0.5*(g*hL*hL - g*h*h); // + g*hL*hbL; 
}

double grassburgers::eval(double* F, double* uL, double* uR, double* qL, double* qR){
  //Enquist Osher Flux for Burger's equation.
  int j;
  double uuL=uL[0];
  double hbL=uL[1];
  double uuR=uR[0];
  double hbR=uR[1];
  double Q = 0.4754;
  //  double R = 
  double H = 1;
  double A = Grass_const;
  double B = Grass_exponent;
  double ws;
//   for (j=0;j<1;j++){
//     if (uuL<0 & uuR<0){
//       F[j]=0.5*uuR*uuR;
//       ws = abs(uuR);
//     }
//     if (uuL>0 & uuR>0){
//       F[j]=0.5*uuL*uuL;
//       ws = abs(uuR);
//     }
//     if (uuL<0 & uuR>0){
//       F[j]=0;
//       ws = 0;
//     }
//     if (uuL>0 & uuR<0){
//       F[j]=0.25*(uuL*uuL+uuR*uuR);
//       ws = abs(0.5*(uuL+uuR));
//     }
//   }
  F[0]=0;
  for (j=1;j<element::systemsize1;j++){
    if (Q > 0){
      //  F[j]=-(3*pow(uuL,3))*log(3*pow(uuL,2)+H-hbL);
      //      ws = (3*pow(uuL,3))/(3*pow(uuL,2)+H-hbL);
      //old approximation
      F[j]=A*Q*abs(pow(Q,B-1))/(pow((H-uL[j]),B));
      ws = A*B*Q*abs(pow(Q,B-1))/(pow((H-uL[j]),B+1));
    }
    if (Q <= 0){
      //      F[j]=-(3*pow(uuR,3))*log(3*pow(uuR,2)+H-hbR);
      //      ws = (3*pow(uuR,3))/(3*pow(uuR,2)+H-hbR);
      //old approximation
      F[j]=A*Q*abs(pow(Q,B-1))/(pow((H-uR[j]),B));
      ws = A*B*Q*abs(pow(Q,B-1))/(pow((H-uR[j]),B+1));
    }
  }
  return ws;
}


double suspendedsediment::eval(double* F, double* uL, double* uR, double* qL, double* qR){
  //HLLC Flux for the shallow water equations 
  int j; 
  double hL,mL,hR,mR,uuL,uuR,SL,SR,SM,hM,mM,chiL,chiR,chiM,hbL,hbR,uu;
  double g = gravity_const; // scaled = H_0 / (U_0)^2 g;
  for (j = 0; j < element::systemsize1; j++){  
    F[j]=0.0;
  }
  hL=uL[0];   
  mL=uL[1];   
  uuL=mL/hL;
  hR=uR[0];
  mR=uR[1];
  uuR=mR/hR;
  hbL=uL[2];  
  hbR=uR[2];
  chiL = uL[3]/hL;
  chiR = uR[3]/hR;
  //Solution for discontinuous bottom.
  if (hbL>hbR){
    hR -= (hbL-hbR);
    mR = hR*uuR;
  }
  else{
    hL -= (hbR-hbL); 
    mL = hL*uuL;
  }
  if (uuL-( sqrt(g*hL))<uuR-(sqrt(g*hR)))
    {
      SL=uuL-sqrt(g*hL);
    }
  else 
    {
      SL=uuR-sqrt(g*hR);
    }
  if (uuL+sqrt(g*hL)>uuR+sqrt(g*hR))
    {
      SR=uuL+sqrt(g*hL);
    }
  else 
    {
      SR=uuR+sqrt(g*hR);
    }   
  double ws = abs(SL); 
  if (abs(SR) > ws){
    ws = abs(SR);
  }
  SM=(0.5*g*hL*hL-0.5*g*hR*hR-uuL*hL*(SL-uuL)+uuR*hR*(SR-uuR))/(hR*(SR-uuR)-hL*(SL-uuL));
  if (SL > 0){
      F[0] = mL;
      F[1] = mL*mL/hL+0.5*g*hL*hL;
      F[2] = 0;
      F[3] = mL*chiL;
    }
  if (SR < 0){
      F[0] = mR;
      F[1] = mR*mR/hR+0.5*g*hR*hR;
      F[2] = 0;
      F[3] = mR*chiR;
    }
  if (SL <= 0 and SM > 0){
      hM = hL*(SL-uuL)/(SL-SM);  //sqrt(2/g*(0.5*g*hL*hL+hL*(SL-uuL)*(SM-uuL)));
      mM = SM*hM;
      chiM = hL*chiL*(SL-uuL)/(hM*(SL-SM));
      F[0] = mL+SL*(hM-hL);
      F[1] = mL*mL/hL+0.5*g*hL*hL+SL*(mM-mL);
      F[2] = 0;
      F[3] = mL*chiL+SL*(chiM*hM-chiL*hL);
  }
  if (SM <= 0 and SR >= 0){
      hM = hL*(SL-uuL)/(SL-SM);  //sqrt(2.0/g*(0.5*g*hR*hR+hR*(SR-uuR)*(SM-uuR)));
      mM = SM*hM;
      chiM = hR*chiR*(SR-uuR)/(hM*(SR-SM));
      F[0] = mR+SR*(hM-hR);
      F[1] = mR*mR/hR+0.5*g*hR*hR+SR*(mM-mR);
      F[2] = 0;
      F[3] = mR*chiR+SR*(chiM*hM-chiR*hR);
  }
  //  F[2] = 0;
  //  for (j = 2; j < element::systemsize1; j++){
  //    F[j]=0;
  //  }
  return ws;
}  

void suspendedsediment::eval_l(double* Fl, double* uL, double* uR, double* qL, double* qR){
  int j; 
  double hL,mL,hR,mR,uuL,uuR,SL,SR,SM,hM,mM,hbL,hbR,h,hb;
  double g = gravity_const; // scaled = H_0 / (U_0)^2 g;
  for (j = 0; j < element::systemsize1; j++){  
    Fl[j]=0;
  }
  hL=uL[0];
  mL=uL[1];
  uuL=mL/hL;
  hR=uR[0];
  mR=uR[1];
  uuR=mR/hR;
  hbL=uL[2];  
  hbR=uR[2];
  if (hbL>hbR){
    h = hR - (hbL-hbR);
    hb = hbL - hbR;
  }
  else {
    h = hR;
    hb = 0;
  }
  Fl[1]+=0.5*(g*hR*hR - g*h*h); //+ g*hR*hbR; 
}

void suspendedsediment::eval_r(double* Fr, double * uL, double * uR, double* qL, double* qR){
  int j;
  double hL,mL,hR,mR,uuL,uuR,SL,SR,SM,hM,mM,hbL,hbR,h,hb;
  double g = gravity_const; // scaled = H_0 / (U_0)^2 g;
  for (j = 0; j < element::systemsize1; j++){  
    Fr[j]=0;
  }
  hL=uL[0];
  mL=uL[1];
  uuL=mL/hL;
  hR=uR[0];
  mR=uR[1];
  uuR=mR/hR;
  hbL=uL[2];  
  hbR=uR[2];
  //Bouchut
  if (hbR>hbL){
    h = hL - (hbR-hbL);
    hb = (hbR-hbL);
  }
  else{
    h = hL; 
    hb = 0;
  }
  Fr[1]+=0.5*(g*hL*hL - g*h*h); // + g*hL*hbL; 
}

double lfgrassmomentumdiffusion::eval(double* F, double* uL, double* uR, double* qL, double* qR){
  //HLL flux for grass bed updating equation with SWE in shock preserving form.
  double g = gravity_const; // scaled = H_0 / (U_0)^2 g;
  double pi = 4.0*atan(1.0);
  int i;
  double A = Grass_const;
  double B = Grass_exponent;
  double dL;
  double dR;
  double aL[4];
  double aR[4];
  double eig[6];
  double hL=uL[0];   
  double hR=uR[0];
  double mL=uL[1];   
  double mR=uR[1];
  double uuL = mL/hL;
  double uuR = mR/hR;
  double hbL=uL[2];  
  double hbR=uR[2];

  //Solution for discontinuous topography.
  if (hbL>hbR){
    hR -= (hbL-hbR);
    mR = hR*uuR;
  }
  else{
    hL -= (hbR-hbL); 
    mL = hL*uuL;
  }

  //characteristic polynomial used for calculating eigenvalues of left face.   
  dL=B*A*abs(pow(uuL,B-1));
  aL[0]=1;
  aL[1]=-2*uuL;
  aL[2]=-g*hL-g*dL+pow(uuL,2);
  aL[3]=g*uuL*dL;
  double QL = (3*aL[2]-pow(aL[1],2))/9;
  double RL = (9*aL[2]*aL[1]-27*aL[3]-2*pow(aL[1],3))/54;
  double theta = acos(RL/sqrt(pow(-QL,3)));
  eig[0] = 2*sqrt(-QL)*cos((theta)/3)-aL[1]/3;
  eig[1] = 2*sqrt(-QL)*cos((theta+2*pi)/3)-aL[1]/3;
  eig[2] = 2*sqrt(-QL)*cos((theta+4*pi)/3)-aL[1]/3;

  //characteristic polynomial used for calculating eigenvalues of right face.   
  dR=B*A*abs(pow(uuR,B-1));
  aR[0]=1;
  aR[1]=-2*uuR;
  aR[2]=-g*hR-g*dR+pow(uuR,2);
  aR[3]=g*uuR*dR;
  double QR = (3*aR[2]-pow(aR[1],2))/9;
  double RR = (9*aR[2]*aR[1]-27*aR[3]-2*pow(aR[1],3))/54;
  theta = acos(RR/sqrt(pow(-QR,3)));
  //  cout << QR << RR << theta << "\n";
  eig[3] = 2*sqrt(-QR)*cos((theta)/3)-aR[1]/3;
  eig[4] = 2*sqrt(-QR)*cos((theta+2*pi)/3)-aR[1]/3;
  eig[5] = 2*sqrt(-QR)*cos((theta+4*pi)/3)-aR[1]/3;
  double SR = eig[0];
  for (i=1;i<6;i++){
    if (abs(eig[i])>abs(SR)){
      SR = eig[i];
    }
  }
  F[0]=0.5*(hL*uuL+hR*uuR)-SR*(hR-hL)/2;
  F[1]=0.5*(hL*pow(uuL,2)+g*pow(hL,2)/2+hR*pow(uuR,2)+g*pow(hR,2)/2) -SR*(mR-mL)/2;
  F[2]=-SR*(hbR-hbL)/2;
  if (uuL >= 0){
    F[2] += 0.5*(A*(1-0.1*qL[0])*abs(pow(uuL,B)));
  }
  else {
    F[2] += 0.5*(A*(-1-0.1*qL[0])*abs(pow(uuL,B)));
  }
  if (uuR >= 0){
    F[2] += 0.5*(A*(1-0.1*qL[0])*abs(pow(uuR,B)));
  }
  else {
    F[2] += 0.5*(A*(-1-0.1*qL[0])*abs(pow(uuR,B)));
  }
  return (abs(SR)+abs(qL[0]));
}

void lfgrassmomentumdiffusion::eval_l(double* Fl, double* uL, double* uR, double* qL, double* qR){
  int j; 
  double hL,mL,hR,mR,uuL,uuR,SL,SR,SM,hM,mM,hbL,hbR,h,hb;
  double g = gravity_const; // scaled = H_0 / (U_0)^2 g;
  for (j = 0; j < element::systemsize1; j++){  
    Fl[j]=0;
  }
  hL=uL[0];
  mL=uL[1];
  uuL=mL/hL;
  hR=uR[0];
  mR=uR[1];
  uuR=mR/hR;
  hbL=uL[2];  
  hbR=uR[2];
  if (hbL>hbR){
    h = hR - (hbL-hbR);
    hb = hbL - hbR;
  }
  else {
    h = hR;
    hb = 0;
  }
  Fl[1]+=0.5*(g*hR*hR - g*h*h); //+ g*hR*hbR; 
}

void lfgrassmomentumdiffusion::eval_r(double* Fr, double * uL, double * uR, double* qL, double* qR){
  int j;
  double hL,mL,hR,mR,uuL,uuR,SL,SR,SM,hM,mM,hbL,hbR,h,hb;
  double g = gravity_const; // scaled = H_0 / (U_0)^2 g;
  for (j = 0; j < element::systemsize1; j++){  
    Fr[j]=0;
  }
  hL=uL[0];
  mL=uL[1];
  uuL=mL/hL;
  hR=uR[0];
  mR=uR[1];
  uuR=mR/hR;
  hbL=uL[2];  
  hbR=uR[2];
  //Bouchut
  if (hbR>hbL){
    h = hL - (hbR-hbL);
    hb = (hbR-hbL);
  }
  else{
    h = hL; 
    hb = 0;
  }
  Fr[1]+=0.5*(g*hL*hL - g*h*h); // + g*hL*hbL; 
}
