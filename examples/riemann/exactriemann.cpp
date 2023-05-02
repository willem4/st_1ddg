/*-------------------------------------------------------------------------------------*/
/* Solution of Riemann Problem for Shallow Water Equations                             */
/*                  --> Vijaya Raghav Ambati                                           */
/*-------------------------------------------------------------------------------------*/

#include<iostream>
#include<fstream>
#include<vector>
#include<math.h>
#include<unistd.h>
using namespace std;

static double g=9.81;
//      Functions to calculate h^star.
double Function1(double hk, double h, int n);
double Function2(double hk, double h, int n);

double exactriemann(double hl, double ul, double hr, double ur, double xs, double x, double t, int j){
//double xl,xr,xs,
//double x,t;
//double hl,hr,ul,ur,
  double hs,us,fl,fr,dfl,dfr,new_hs,h,u;
//   string file ="riemann_exact.dat";
//   ofstream output(file.c_str());
//   system("tput clear");
//   cout<<"Enter xl : ";
//   cin>>xl; // left end of domain
//   cout<<"Enter xr : ";
//   cin>>xr; // right end of domain
//   cout<<"Enter xs : ";
//   cin>>xs;//point of discontinuity
//   cout<<"Enter hl : ";
//   cin>>hl; // hl
//   cout<<"\nEnter hr : ";
//   cin>>hr; // hr
//   cout<<"\nEnter ul : ";
//   cin>>ul; // ul
//   cout<<"\nEnter ur : ";
//   cin>>ur; // ur.
//   //      Assume h star = hleft
  new_hs=hl;
  
  // Iterative Procedure to obtain (h^star) using
  do{
    hs=new_hs;
    if(hs>hl){
      if(hs>hr){
	// hs > hl so use shock wave function
	fl =Function2(hl,hs,0);
	dfl=Function2(hl,hs,1);
	// hs > hr so use shock wave function
	fr =Function2(hr,hs,0);
	dfr=Function2(hr,hs,1);
	
      }
      else{
	// hs > hl so use shock wave function
	fl =Function2(hl,hs,0);
	dfl=Function2(hl,hs,1);
	// hs < hr so use rarefaction wave function
	fr =Function1(hr,hs,0);
	dfr=Function1(hr,hs,1);
      }
    }
    else{
      if(hs>hr){
	// hs < hl so use rarefaction wave function
	fl =Function1(hl,hs,0);
	dfl=Function1(hl,hs,1);
	// hs > hr so use shock wave function
	fr =Function2(hr,hs,0);
	dfr=Function2(hr,hs,1);
	
      }
      else{
	// hs < hl so use rarefaction wave function
	fl =Function1(hl,hs,0);
	dfl=Function1(hl,hs,1);
// hs < hr so use rarefaction wave functon;
	fr =Function1(hr,hs,0);
	dfr=Function1(hr,hs,1);
      }
    }
    new_hs=hs-((fl+fr+ur-ul)/(dfl+dfr));
  }while(fabs(new_hs-hs)>1e-16);
  
// Calculation of h star is finished.
  
// Calculate u star
//   cout << endl;
//   cout <<"|--------------------------|\n";
//   cout <<"|  Solution of star state  |\n";
//   cout <<"|--------------------------|\n";
//   cout <<" h star = "<<hs<<endl;
  us=(ul+ur)/2+(fr-fl)/2;
  //        cout <<" u star = "<<us<<endl;
  //        cout <<"|--------------------------|\n";
  
       // Enter x and t:
  //        x = 0;
//        t = 0.1;
  // Exact solution
  double al=sqrt(g*hl);
  double ar=sqrt(g*hr);
  double as=sqrt(g*hs);
  // Left wave:
  if(hs > hl){
    double sl = (hs * us - hl * ul) / (hs - hl);
    // a) Shock wave
    if(x < xs + sl * t){
      h = hl;
      u = ul;
    }
    else{
      h = hs;
      u = us;
    }
  }
  else{
    // b) Rarefaction wave
    double shl = ul - al;
    double stl = us - as;
    if(x < xs + shl * t){
      h = hl;
      u = ul;
    }
    else{
      if(x < xs + stl * t){
	h = (1/(9.0 * g)) * pow(ul + 2 * al - (x-xs)/t,2);
	u = ul + 2.0 * al - 2.0 * sqrt(g*h);
      }
      else{
	h = hs;
	u = us;
      }
    }
  }
  // Right Wave:
  if(hs > hr){
    // a) Shock wave
    double sr = (hs * us - hr * ur) / (hs - hr);
    if(x > xs + sr * t){
      h = hr;
      u = ur;
    }
  }
  else{
    // b) Rarefaction wave
    double shr = us + as;
    double str = ur + ar;
    if(x > xs + shr * t){
      if( x < xs + str * t){
	h = (1/(9.0 * g)) * pow(2 * ar - ur + (x - xs)/t,2);
	u = ur + 2.0 * sqrt(g * h) - 2.0 * ar;
      }
      else{
	h = hr;
	u = ur;
      }
    }
  }
  switch(j) {
  case 0:
    return h;
    break;
  case 1:
    return h*u;
    break;
  }
  //       cout << "Height = " << h << endl;
  //       cout << "velocity = " << u << endl;
}

// Funtion for Rarefaction Wave
double Function1(double hk, double h, int n){
  double f,df;
  if(n==0){
    f=2.0*sqrt(g*h)-2.0*sqrt(g*hk);
    return(f);
  }
  else{
    df=sqrt(g/h);
    return(df);
  }
}
// Function for Shock Wave
double Function2(double hk, double h, int n){
  double f,df;
  if(n==0){
    f=sqrt((g*pow((h-hk),2)*(h+hk))/(2.0*h*hk));
    return(f);
  }
  else{
    df=sqrt(g/(2*hk))*((2+(hk/h)+pow((hk/h),2))/(2*sqrt(1+(hk/h))));
    return(df);
  }
}

// int main(){
//   double hl = 2.0; 
//   double ul = 0.5;
//   double hr = 1.0;
//   double ur = 0.0;
//   double xs = 0.5;
//   double x = 0.93;
//   double t = 0.1;
//   int j = 1;
//   cout << exactriemann(hl, ul, hr, ur, xs, x, t, j) << "\n";
// }
