#include "node.h"
#include "state.h"
#include "element.h"
#include "flux.h"
#include "source.h"
#include "flux2.h"

#include<vector>

using namespace std;

//! Runge Kutta Third Order Time Integration.
extern int rk3(double &dt, vector<node> &Nodes, vector<element>& Elements, int bcl, int bcr, testtypedescr testtype, fluxbase * fluxtype, sourcebase * sourcetype, flux2base * flux2type, double t);
