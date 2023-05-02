#ifndef FLUX_H
#define FLUX_H

using namespace std;

/* extern int advection(double *,double *, double *); */
/* extern int burgers(double *, double *, double *); */
/* extern int flux(double *, double *, double *); */
/* extern int fluxl(double *, double *, double *); */
/* extern int fluxr(double *, double *, double *); */
/* extern int flux(double *, double *); */

//! Class Fluxbase 
class fluxbase {
public:
  //! Virtual evaluation of the flux.
  virtual double eval(double *, double *, double *, double *, double *) = 0;
  //! Virtual non conservative evaluation of the flux on the left of the element.
  virtual void eval_l(double *, double *, double *, double *, double *) = 0;
  //! Virtual non conservative evaluation of the flux on the right of the element.
  virtual void eval_r(double *, double *, double *, double *, double *) = 0;
};

//! Derived fluxbase class for Advection.
class advection: public fluxbase {
public:
  //! Evaluation of the advection flux.
  double eval(double *, double *, double *, double *, double *);
  void eval_l(double *, double *, double *, double *, double *){}
  void eval_r(double *, double *, double *, double *, double *){}
};

//! Derived fluxbase class for Burgers' equation.
class burgers: public fluxbase {
public:
  //! Evaluation of the flux for Burgers' equation.
  double eval(double *, double *, double *, double *, double *);
  void eval_l(double *, double *, double *, double *, double *){}
  void eval_r(double *, double *, double *, double *, double *){}
};

//! Derived fluxbase class for Burgers' equation.
class burgersdiffusion: public fluxbase {
public:
  //! Evaluation of the flux for Burgers' equation including diffusion.
  double eval(double *, double *, double *, double *, double *);
  void eval_l(double *, double *, double *, double *, double *){}
  void eval_r(double *, double *, double *, double *, double *){}
};

//! Derived fluxbase class for Shallow Water equations using HLLC with no topography.
class swehllc: public fluxbase {
public:
  //! Evaluation of the HLLC flux.
  double eval(double *, double *, double *, double *, double *);
  void eval_l(double *, double *, double *, double *, double *){}
  void eval_r(double *, double *, double *, double *, double *){}
};

//! Derived fluxbase class for Shallow Water equations using HLLC with discontinuous topography.
class swehllctopography: public fluxbase {
public:
  //! Evaluation of the HLLC flux with topography correction.
  double eval(double *, double *, double *, double *, double *);
  //! Non conservative flux as suggested by Bouchut et al. (2004).
  void eval_l(double *, double *, double *, double *, double *);
  //! Non conservative flux as suggested by Bouchut et al. (2004).
  void eval_r(double *, double *, double *, double *, double *);
};


//! Derived fluxbase class for sediment transport.
class sedimenttransport: public fluxbase {
public:
  //! Evaluation of the ... flux with topography correction.
  double eval(double *, double *, double *, double *, double *);
  //! Non conservative flux as suggested by Bouchut et al. (2004).
  void eval_l(double *, double *, double *, double *, double *);
  //! Non conservative flux as suggested by Bouchut et al. (2004).
  void eval_r(double *, double *, double *, double *, double *);
};

//! Derived fluxbase class for assymtotic formulation of the Shallow Water equations including Grass' Bed Updating Equation using HLL.
class hllgrass: public fluxbase {
public:
  //! Evaluation of the HLL flux.
  double eval(double *, double *, double *, double *, double *);
  void eval_l(double *, double *, double *, double *, double *){}
  void eval_r(double *, double *, double *, double *, double *){}
};

//! Derived fluxbase class for assymtotic formulation of the Shallow Water equations including Grass' Bed Updating Equation using LF.
class lfgrass: public fluxbase {
public:
  //! Evaluation of the Lax Friedrichs flux.
  double eval(double *, double *, double *, double *, double *);
  void eval_l(double *, double *, double *, double *, double *){}
  void eval_r(double *, double *, double *, double *, double *){}
};

//! Derived fluxbase class for the Shallow Water equations including Grass' Bed Updating Equation using LF.
class lfgrassmomentum: public fluxbase {
public:
  //! Evaluation of the Lax Friedrichs flux with topography correction.
  double eval(double *, double *, double *, double *, double *);
  //! Non conservative flux as suggested by Bouchut et al. (2004).
  void eval_l(double *, double *, double *, double *, double *);
  //! Non conservative flux as suggested by Bouchut et al. (2004).
  void eval_r(double *, double *, double *, double *, double *);
};

class grassburgers: public fluxbase {
public:
  //! Evaluation of the assymtotic Grass bed updating equation as a sort of burgers equation.
  double eval(double *, double *, double *, double *, double *);
  void eval_l(double *, double *, double *, double *, double *){}
  void eval_r(double *, double *, double *, double *, double *){}
};

class suspendedsediment: public fluxbase {
public:
  //! Evaluation of the HLLC flux with topography correction.
  double eval(double *, double *, double *, double *, double *);
  //! Non conservative flux as suggested by Bouchut et al. (2004).
  void eval_l(double *, double *, double *, double *, double *);
  //! Non conservative flux as suggested by Bouchut et al. (2004).
  void eval_r(double *, double *, double *, double *, double *);
};

//! Derived fluxbase class for the Shallow Water equations including Grass' Bed Updating Equation using LF and LDG for diffusion term.
class lfgrassmomentumdiffusion: public fluxbase {
public:
  //! Evaluation of the Lax Friedrichs flux with topography correction.
  double eval(double *, double *, double *, double *, double *);
  //! Non conservative flux as suggested by Bouchut et al. (2004).
  void eval_l(double *, double *, double *, double *, double *);
  //! Non conservative flux as suggested by Bouchut et al. (2004).
  void eval_r(double *, double *, double *, double *, double *);
};

#endif
