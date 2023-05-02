#include "node.h"
#include "element.h"

using namespace std;

// Static data members
// Construtor : To construct the node numbers of the state class.
node::node(){
  double x;
}

// Overloading Constructor.
node::node(double value){
    x=value; // Assigns value to x.
    uL = new double [element::systemsize1]; //value of system U left of face.
    uR = new double [element::systemsize1]; //value of system U right of face.
    qL = new double [element::systemsize2]; //value of system Q left of face.
    qR = new double [element::systemsize2]; //value of system Q right of face.
}
