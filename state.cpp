#include"state.h"

using namespace std;

// Static data members
int state::statesize;
state state::g;

// Construtor : To construct the node numbers of the state class.
state::state(){
  coefficient=new double [state::statesize];
  for (int i=0; i<state::statesize; i++){
    coefficient[i]=0;
  }
}

//Destructor
// state::~state(){
//    delete [] coefficient; 
// }

state state::operator=(state vector){
  for(int i=0; i<state::statesize; i++){
    coefficient[i]=vector.coefficient[i];
  }
  return *this;
}

double state::get(double xi){
  double f = coefficient[0];
  if (state::statesize ==2){ 
    f+=xi*coefficient[1];
  }
  return f;
}

double state::getslope(double xi){
  double f = 0;
  if (state::statesize ==2){ 
    f+=coefficient[1];
  }
  return f;
}
  
state state::clear(){
  for(int i=0; i<state::statesize; i++){
    coefficient[i]=0;
  }
  return *this;
}

state state::operator+=(state vector){
  for(int i=0; i<state::statesize; i++){
    coefficient[i]+=vector.coefficient[i];
  }
  return *this;
}

state state::operator*(double alpha){
  for(int i=0; i<state::statesize; i++){
    g.coefficient[i]=alpha*coefficient[i];
  }
  return g;
}

state state::operator/(double alpha){
  state g;
  for(int i=0; i<state::statesize; i++){
    g.coefficient[i]=coefficient[i]/alpha;
  }
  return g;
}
  
