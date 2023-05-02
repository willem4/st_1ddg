#ifndef SOURCE_H
#define SOURCE_H

using namespace std;

//! Source base class.
class sourcebase {
public:
  //! Virtual evaluation of the source. 
  virtual void eval(double *, double *, double *, double *) = 0;
};

//! Derived NoSource Class from SourceBase
class nosource: public sourcebase {
public:
  //!Evaluation of the source which does nothing.
  void eval(double *, double *, double *, double *){}
};

//! Derived Topography Class from SourceBase.
class topography: public sourcebase {
public:
  //! Evaluation of the Source in the Shallow Water Equations.
  void eval(double *, double *, double *, double *);
};

//! Derived Exchange Class from SourceBase.
class exchange: public sourcebase {
public:
  //! Evaluation of the Source in the Shallow Water Equations and Bed Sediment and Suspended Sediment Exchange.
  void eval(double *, double *, double *, double *);
};
#endif
