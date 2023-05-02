#ifndef FLUX2_H
#define FLUX2_H

using namespace std;

//! Local Dicontinuous Galerkin Flux Class
class flux2base {
public:
  //! Virtual Evaluation of the LDG flux.
  virtual double eval(double *, double *, double *) = 0;
};

//! LDG noflux class Derived from flux2base. 
class noflux2: public flux2base {
public:
  //! Evaluation which does nothing.
  double eval(double *, double *, double *){}
};

//! LDG diffusion flux class used in comination with Burgers' equation.
class diffusion: public flux2base {
public:
  //! Evaluates the flux for the LDG variables for the difussion.
  double eval(double *, double *, double *);
};

//! LDG diffusion flux class used in comination with Grass Bed Updating equation.
class sedimentdiffusion: public flux2base {
 public:
  //! Evaluates the flux for the LDG variables for the sedimentdifussion.
  double eval(double *, double *, double *);
};

#endif
