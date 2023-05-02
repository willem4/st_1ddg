#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "node.h"
#include "state.h"
#include "element.h"
#include "flux.h"
#include "source.h"
#include "flux2.h"

#include "enumtypes.h"
#include "readtest.h"
#include "readmesh.h"
#include "physpara.h"

#include "initcon.h"
#include "output.h"

#include "flux2proc.h"
#include "eulerforw.h"
#include "rk3.h"
#include "cranknicolsonwp.h"
#include "updatenodes.h"

using namespace std;

string itoa(int i)
{
  int rem = 0;
  string ret;
  while(0 < i)
  {
    rem = i%10;
    i = i/10;
    ret = (char)(rem + '0') + ret;
  }
  return ret;
}

double getdt(vector <element> Elements, testtypedescr testtype)
{ 
  int i;
  double minlength = 123456; //arbitrary initial value.
  double maxspeed = 0;
  double u;
  double h;
  double m;
  double q;
  double g = gravity_const;
  double A = Grass_const;
  switch ( testtype )
    {
    case Advection: 
      for (i = 0; i < Elements.size(); i++){
	if (Elements[i].length < minlength){
	  minlength = Elements[i].length;
	}
      }
      maxspeed = advection_speed;
      break;
    case Burgers:
      for (i = 0; i < Elements.size(); i++){
	if (Elements[i].length < minlength){
	  minlength = Elements[i].length;
	}
      }
      for (i = 0; i < Elements.size(); i++){
	u = Elements[i].U[0].get(-1.0);
	if (u*u*0.5>maxspeed) {
	  maxspeed = u*u*0.5;
	}
	u = Elements[i].U[0].get(1.0);
	if (u*u*0.5>maxspeed) {
	  maxspeed = u*u*0.5;
	}
      }
      break;
    case SWEBurgers: 
    case SWELinearWaveSolution: 
    case SWERiemannProblemLeftRarefactionRightShock:
    case SWERiemannProblemLeftShockRightRarefaction: 
    case SWERiemannProblemLeftShockRightShock: 
    case SWERiemannProblemLeftRarefactionRightRarefaction:
    case SWEContinuousTopography: 
    case SWEDiscontinuousTopography: 
    case SWEFlowOverIsolatedContinuousRidgeI: 
    case SWEFlowOverIsolatedDiscontinuousRidgeI: 
    case SWEFlowOverIsolatedContinuousRidgeIV:
    case SWEFlowOverIsolatedDiscontinuousRidgeIV:
    case SWEFlowOverSinyBedI:
      for (i = 0; i < Elements.size(); i++){
	if (Elements[i].length < minlength){
	  minlength = Elements[i].length;
	}
      }
      for (i = 0; i < Elements.size(); i++){
	h = Elements[i].U[0].get(-1.0);
	m = Elements[i].U[1].get(-1.0);
	u = m/h;
	if (sqrt(u*u+g*h)>maxspeed) {
	  maxspeed = sqrt(u*u+g*h);
	}
	h = Elements[i].U[0].get(1.0);
	m = Elements[i].U[1].get(1.0);
	u = m/h;
	if (sqrt(u*u+g*h)>maxspeed) {
	  maxspeed = sqrt(u*u+g*h);
	}
      }
      break;
    case BurgersDiffusive:
      for (i = 0; i < Elements.size(); i++){
	if (Elements[i].length < minlength){
	  minlength = Elements[i].length;
	}
      }
      for (i = 0; i < Elements.size(); i++){
	u = Elements[i].U[0].get(-1.0);
	if (u*u*0.5+sqrt(u*u)/Elements[i].length>maxspeed) {
	  maxspeed = u*u*0.5+abs(u)/Elements[i].length;
	}
	u = Elements[i].U[0].get(1.0);
	if (u*u*0.5+sqrt(u*u)/Elements[i].length>maxspeed) {
	  maxspeed = u*u*0.5+abs(u)/Elements[i].length;
	}
      }
      break;
    case SedimentTransport:
    case HLLGrass: //Better prediction possible. Pass from flux.
    case LFGrass: //Better prediction possible. Pass from flux.
      for (i = 0; i < Elements.size(); i++){
	if (Elements[i].length < minlength){
	  minlength = Elements[i].length;
	}
      }
      for (i = 0; i < Elements.size(); i++){
	h = Elements[i].U[0].get(-1.0);
	m = Elements[i].U[1].get(-1.0);
	u = m/h;
	if (sqrt(u*u+g*h)>maxspeed) {
	  maxspeed = sqrt(u*u+g*h);
	}
// 	if (abs(A*u)>maxspeed) {
// 	  maxspeed = abs(A*u);
// 	}
	h = Elements[i].U[0].get(1.0);
	m = Elements[i].U[1].get(1.0);
	u = m/h;
	if (sqrt(u*u+g*h)>maxspeed) {
	  maxspeed = sqrt(u*u+g*h);
	}
// 	if (abs(A*u)>maxspeed) {
// 	  maxspeed = abs(A*u);
// 	}
      }
      break;
    }
  return minlength/maxspeed;
}

int main(int argc, char * argv[])
{
  if(argc!=2)
  {
    cout<<"Correct input is: st_1ddg testfile\n";
    return 0;
  }
  /*  if(access(argv[1], 00)) //access returns 0 if the file can be accessed
  {                //under the specified method (00)
    cout<<"Testfile does not exist\n";  //because it checks file existence
    return 0;
    }*/

  timeintmethod timeint;
  testtypedescr testtype;
  string meshname;
  double it;          //initial time.
  double dt;          //time step.
  double et;          //final time.
  int deg;            //degree of functions 1 or 2.
  string ext;         //file extension. Used for Error Testing. 
  int bcl, bcr;       //Left and Right Boundary Conditions.
  int count = 0;
  int countmax = 100;
  double CFL; 

  readtest2(testtype, timeint, meshname, CFL, it, et, deg, argv[1]);

  //Test for reading and showing contents of file.
  //  ifstream testfile;  //ifstream is used for file input
  //  testfile.open(argv[1]); //argv[1] is the second argument passed in
                         //presumably the file name
  //char x;
  //the_file.get(x);
  //while(!the_file.eof())  //eof is defined as the end of the file
  //{
  //  cout<<x;
  //  the_file.get(x);//Notice we always let the loop check x for the end of
  //}  //file to avoid bad output when it is reached
  //the_file.close();  //Always clean up
  //return 0;
  // argc:  Number of items in argv 
  // argv:  Address of an array of C-strings 
  //  int Nn; //Number of Nodes
  //  int Ne = 10; //Number of Elements 
  
  vector <element> Elements;
  vector <node> Nodes;

  element::CFL = CFL;
  state::statesize = deg;      //degree of elements (1 or 2). 

  //Determine the systemsizes 
  switch ( testtype )
    {
    case Advection: 
    case Burgers:
      element::systemsize1 = 1;  //systems of the form u_t + f(u)_x = s
      element::systemsize2 = 0;  //systems of the form u + g(u)_x = s 
      break;
    case SWEBurgers: 
    case SWELinearWaveSolution: 
    case SWERiemannProblemLeftRarefactionRightShock:
    case SWERiemannProblemLeftShockRightRarefaction: 
    case SWERiemannProblemLeftShockRightShock: 
    case SWERiemannProblemLeftRarefactionRightRarefaction:
    case DecoupledGrassSedimentOnly:
      element::systemsize1 = 2;  //systems of the form u_t + f(u)_x = s
      element::systemsize2 = 0;  //systems of the form u + g(u)_x = s 
      break;
    case SWEContinuousTopography: 
    case SWEDiscontinuousTopography: 
    case SWEFlowOverIsolatedContinuousRidgeI: 
    case SWEFlowOverIsolatedDiscontinuousRidgeI: 
    case SWEFlowOverIsolatedContinuousRidgeIV:
    case SWEFlowOverIsolatedDiscontinuousRidgeIV: 
    case SWEFlowOverSinyBedI:
    case SedimentTransport:
    case HLLGrass:
    case LFGrass:
    case LFGrassMomentum:
      element::systemsize1 = 3;  //systems of the form u_t + f(u)_x = s
      element::systemsize2 = 0;  //systems of the form u + g(u)_x = s 
      break;
    case BurgersDiffusive:
      element::systemsize1 = 1;  //systems of the form u_t + f(u)_x = s
      element::systemsize2 = 1;  //systems of the form u + g(u)_x = s 
      break;
    case SuspendedSediment:
      element::systemsize1 = 4;  //systems of the form u_t + f(u)_x = s
      element::systemsize2 = 0;  //systems of the form u + g(u)_x = s 
      break;
    case LFGrassMomentumDiffusion:
      element::systemsize1 = 3;  //systems of the form u_t + f(u)_x = s
      element::systemsize2 = 1;  //systems of the form u + g(u)_x = s 
      break;      
    }

  //Choose the flux, source and flux2 to be used.
  fluxbase * fluxtype;
  sourcebase * sourcetype;
  flux2base * flux2type;
  switch ( testtype )
    {
    case Advection: 
      cout << "Program for Advection equation using upwind flux\n";
      fluxtype = new advection;
      sourcetype = new nosource;
      flux2type = new noflux2;
      break;
    case Burgers:
      cout << "Program for Burgers equation using Engquist-Osher flux\n";
      fluxtype = new burgers;
      sourcetype = new nosource;
      flux2type = new noflux2;
      break;
    case SWEBurgers: 
    case SWELinearWaveSolution: 
    case SWERiemannProblemLeftRarefactionRightShock:
    case SWERiemannProblemLeftShockRightRarefaction: 
    case SWERiemannProblemLeftShockRightShock: 
    case SWERiemannProblemLeftRarefactionRightRarefaction:
      cout << "Program for Shallow Water Equations - HLLC flux\n";
      fluxtype = new swehllc;
      sourcetype = new nosource;
      flux2type = new noflux2;
      break;
    case SWEContinuousTopography: 
    case SWEDiscontinuousTopography: 
    case SWEFlowOverIsolatedContinuousRidgeI: 
    case SWEFlowOverIsolatedDiscontinuousRidgeI: 
    case SWEFlowOverIsolatedContinuousRidgeIV:
    case SWEFlowOverIsolatedDiscontinuousRidgeIV: 
    case SWEFlowOverSinyBedI:
      fluxtype = new swehllctopography;
      sourcetype = new topography;
      flux2type = new noflux2;
      cout << "Program for Shallow Water Equations - HLLC flux \n";
      cout << "Including Source Term\n";
      break;
    case SedimentTransport:
      fluxtype = new sedimenttransport;
      sourcetype = new topography;
      flux2type = new noflux2;
      cout << "Program for Sediment Transport ";
      cout << "Shallow Water Equations including source term - HLLC flux & \n";
      cout << "Upwind treatment of sediment movement. \n";
      break;
    case BurgersDiffusive:
      fluxtype = new burgersdiffusion;
      sourcetype = new nosource;
      flux2type = new diffusion;
      cout << "Program for Burger's equation using Engquist Osher flux \n";
      cout << "(Including diffusive term)\n";
      break;
    case HLLGrass:
      fluxtype = new hllgrass;
      sourcetype = new nosource;      
      flux2type = new noflux2;
      cout << "Program for Grass bed updating equation with \nSWE in non shock preserving form\n using HLL flux.\n";
      break;
    case LFGrass:
      fluxtype = new lfgrass;
      sourcetype = new nosource;      
      flux2type = new noflux2;
      cout << "Program for Grass bed updating equation with \nSWE in non shock preserving form\n using Lax-Friedrichs flux.\n";
      break;
    case LFGrassMomentum:
      fluxtype = new lfgrassmomentum;
      sourcetype = new topography;
      flux2type = new noflux2;
      cout << "Program for Grass bed updating equation with \nSWE in shock preserving form\n using Lax-Friedrichs flux.\n";
      break;
    case DecoupledGrassSedimentOnly:
      fluxtype = new grassburgers;
      sourcetype = new nosource;      
      flux2type = new noflux2;
      cout << "Program for Grass bed updating equation - Sediment Only in the assymtotic case \n";
      break;
    case SuspendedSediment:
      fluxtype = new suspendedsediment;
      sourcetype = new exchange;      
      flux2type = new noflux2;
      cout << "Program for Suspended Sediment \n";
      break;
    case LFGrassMomentumDiffusion:
      fluxtype = new lfgrassmomentumdiffusion;
      sourcetype = new topography;
      flux2type = new sedimentdiffusion;
      cout << "Program for Grass bed updating equation with \nSWE in shock preserving form\n using Lax-Friedrichs flux including LDG diffusion part.\n";
      break;
    }
  
  readmesh(meshname.c_str(), Elements, Nodes, bcl, bcr);

  cout << "Succesfully read meshfile: " << meshname.c_str() << "\n";
  cout << "Number of Elements:        " << Elements.size() << "\n";
  cout << "Number of Nodes:           " << Nodes.size() << "\n";
  cout << "CFL Condition              " << CFL << "\n";
  cout << "Final time:                " << et << "\n";
  cout << "Degree of freedom:         " << state::statesize << "\n";
  cout << "Boundary conditions:       \n";
  cout << "---------------------------" << "\n";
  cout << "Initial condition at:  t = " << it << "\n"; 

  initcon(Nodes, Elements, testtype,it); //Fill Elements with initial conditions 
  ext = itoa(Elements.size());
  output(Nodes, Elements, "init." + ext);

  switch ( timeint ) 
    { 
    case EulerForward:
      cout << "Euler forward time integration: \n";
      while (it<et){ 
	dt = CFL*getdt(Elements,testtype);
	count++;
	if (it+dt > et){
	  dt = et - it;
	}
	it += dt;
	if (count == countmax) {
	  cout << "Evaluating at time: " << it << "\n";
	  count = 0;
	}
	eulerforw(dt, Nodes, Elements, bcl, bcr, testtype, fluxtype, sourcetype, flux2type, it-dt);
      }
      break;
    case RungeKutta3:
      cout << "Runge Kutta third order time integration: \n";
      while (it<et){ 
	dt = CFL*getdt(Elements,testtype);
	count++;
	if (it+dt > et){
	  dt = et - it;
	}
	it += dt;
	if (count == countmax) {
	  cout << "Evaluating at time: " << it << "\n";
	  count = 0;
	}
       	rk3(dt, Nodes, Elements, bcl, bcr, testtype, fluxtype, sourcetype, flux2type, it-dt);
      }
      break;
    case CrankNicolsonWithPredictor:
      cout << "Crank Nicolson scheme with predictor step time integration: \n";
      while (it<et){ 
	dt = CFL*getdt(Elements,testtype);
	count++;
	if (it+dt > et){
	  dt = et - it;
	}
	it += dt;
	if (count == countmax) {
	  cout << "Evaluating at time: " << it << "\n";
	  count = 0;
	}
       	cranknicolsonwp(dt, Nodes, Elements, bcl, bcr, testtype, fluxtype, sourcetype, flux2type, it-dt);
      }
      break;
    } 

  output(Nodes, Elements, "out." + ext);
  cout << "Final result  at:  t = " << it << "\n"; 

//   switch ( timeint ) 
//     { 
//     case EulerForward:
//       cout << "Euler forward time integration: \n";
//       while (it<et){ 
// 	dt = et - it; //CFL*getdt(Elements,testtype);
// 	count++;
// 	//	if (it+dt > et){
// 	//	  dt = et - it;
// 	//	}
// 	eulerforw(dt, Nodes, Elements, bcl, bcr, testtype, fluxtype, sourcetype, flux2type, it-dt);
// 	it += dt;
	
// 	if (count == countmax) {
// 	  cout << "Evaluating at time: " << it << "\n";
// 	  count = 0;
// 	}

//       }
//       break;
//     case RungeKutta3:
//       cout << "Runge Kutta third order time integration: \n";
//       while (it<et){ 
// 	dt = et-it; //CFL*getdt(Elements,testtype);
// 	count++;
// 	//	if (it+dt > et){
// 	//	  dt = et - it;
// 	//	}
//        	rk3(dt, Nodes, Elements, bcl, bcr, testtype, fluxtype, sourcetype, flux2type, it);
// 	it += dt;
	
// 	if (count == countmax) {
// 	  cout << "Evaluating at time: " << it << "\n";
// 	  count = 0;
// 	}

//       }
//       break;
//     case CrankNicolsonWithPredictor:
//       cout << "Crank Nicolson scheme with predictor step time integration: \n";
//       while (it<et){ 
// 	dt = et-it;//CFL*getdt(Elements,testtype);
// 	count++;
// 	//	if (it+dt > et){
// 	//	  dt = et - it;
// 	//	}
//        	cranknicolsonwp(dt, Nodes, Elements, bcl, bcr, testtype, fluxtype, sourcetype, flux2type, it-dt);
// 	it += dt;

// 	if (count == countmax) {
// 	  cout << "Evaluating at time: " << it << "\n";
// 	  count = 0;
// 	}
// 	//	cout << dt; 
//       }
//       break;
//     } 

}
