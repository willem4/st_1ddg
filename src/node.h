#ifndef NODE_H
#define NODE_H

//! Node Class (= Face)
class node{
public:
  //! x-value of the node.
  double x;  
  //! Default constructor: creates x;
  node();    
  //! Overloading constructor assigning x = value;
  node(double);
  //! Value of system U on the left of the face.    
  double * uL; 
  //! Value of system U on the right of the face.    
  double * uR; 
  //! Value of system Q on the left of the face.    
  double * qL; 
  //! Value of system Q on the right of the face.    
  double * qR; 
};

#endif
