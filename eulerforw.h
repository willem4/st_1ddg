#include "node.h"
#include "state.h"
#include "element.h"
#include "flux.h"
#include "source.h"
#include "flux2.h"

#include<vector>

using namespace std;

//! Euler Forward Time Integration.
extern int eulerforw(double &dt, vector<node> &Nodes, vector<element> &Element, int bcl, int bcr, testtypedescr testtype, fluxbase * fluxtype, sourcebase * sourcetype, flux2base * flux2type, double t);
