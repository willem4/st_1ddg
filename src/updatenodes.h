#include <vector>

#include "state.h"
#include "element.h"
#include "node.h"
#include "flux.h"

#include "enumtypes.h"

using namespace std;

extern int updatenodesu(vector<node> &Nodes, vector<element> &Elements, int bcl, int bcr, testtypedescr testtype, double t);
extern int updatenodesq(vector<node> &Nodes, vector<element> &Elements, int bcl, int bcr, testtypedescr testtype, double t);
