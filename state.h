#ifndef STATE_H
#define STATE_H

//! The state class. 
/*! 
  The state describes how each variable in the element class is stored.
*/
class state{
public:
  //! Size of the state vector
  static int statesize;
  //! Temporary state vector used to prevent memory overflow.
  static state g;
  //! Expansion coefficients of state vector
  double* coefficient; 
  //! Default constructor: creates an empty state vector
  state();	     
  //! Asigns one state vector to another state vector.
  state operator=(state vector);
  //! Function used to get the value of the state at \f$\zeta \in [-1,1]\f$.
  double get(double);
  //! Function used to get the value of the first derivative of the state at \f$\zeta \in [-1,1]\f$.
  double getslope(double);
  //! Clears the values in the state vector.
  state clear();
  //! Adds the values of one state vector to another
  state operator+=(state vector);
  //! Multiplies the values in the state vector by \f$ \alpha \f$
  state operator*(double alpha);
  //! Divides the values in the state vector by \f$ \alpha \f$
  state operator/(double alpha);
};

#endif
