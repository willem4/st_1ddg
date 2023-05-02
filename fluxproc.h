#include <vector>

#include "state.h"
#include "element.h"
#include "node.h"
#include "flux.h"

#include "enumtypes.h"

using namespace std;

extern int fluxproc(vector<node> &Nodes, vector<element>& Elements, int bcl, int bcr, testtypedescr testtype, fluxbase * fluxtype, double t);
