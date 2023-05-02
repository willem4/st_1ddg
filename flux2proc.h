#include <vector>

#include "state.h"
#include "element.h"
#include "node.h"
#include "flux2.h"
#include "enumtypes.h"

using namespace std;

extern int flux2proc(vector<node> &Nodes, vector<element>& Elements, int bcl, int bcr, testtypedescr testtype, flux2base * flux2type, double t);

