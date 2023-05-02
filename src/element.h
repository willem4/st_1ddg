#ifndef ELEMENT_H
#define ELEMENT_H

#include "state.h"

//! The element class.
/*!
  The element class is the building block of how the information per element is stored.
*/
class element{
public:
  //! The dimension of the Discontinuous Galerkin variables.
  static int systemsize1;  
  //! The dimension of the Local Discontinuous Galerkin variables. 
  static int systemsize2;
  //! Interpolation point used for two point Gaussian quadrature
  static double ep;
  //! The Courant coefficient number.
  static double CFL;
  //! The node numbers which bound the element.
  int Node1, Node2;
  //! The length of the element.
  double length;
  //! The stored Discontinuous Galerkin variables.
  state * U;
  //! The right hand side used while calculating one time step.
  state * RHSofU;
  //! A backup of U needed for some time integration schemes.
  state * U_old;
  //! A backup of RHSofU needed for some time integration schemes.
  state * RHSofU_old;
  //! The Local Discontinuous Galerkin variables. 
  state * Q;
  //! The right hand side used for Local Discontinuous Galerkin calculations.
  state * RHSofQ;
  //! The value of the stabilisation operator.
  state StabOp;
  //! A backup of the stabilisation operator used in some time integration schemes.
  state StabOp_old;
  //! Discontinuity detector variable from Krividonova et al., 2004.
  double Kriv;
  //! Count of how many faces have information which passes into the cell \f$\| \partial\Omega^- \|\f$.
  int KrivCount;
  //! The maximum wavespeed sorted per element used to determine the time step.
  double ws;
  //! The maximum d_xx wavespeed sorted per element used to determine the time step.
  double dds;
  //! The values of the DG variables in -ep.
  double * u1;
  //! The values of the DG variables in ep.
  double * u2;
  //! The values of slopes of the DG variables in -ep.
  double * u1slope;
  //! The values of slopes of the DG variables in ep.
  double * u2slope;
  //! The values of the LDG variables in -ep.
  double * q1;
  //! The values of the LDG variables in -ep.
  double * q2;
  //! Default constructor.
  element();
  //! Overloaded constructor.
  element(int, int, double);
  //! Updates u1, u2, u1slope and u2slope in each time step.
  void updateu();
  //! Updates q1 and q2 in every time step. 
  void updateq();
  //! Default destructor.
  ~element();
};

#endif

