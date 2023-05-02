#include "state.h"
#include "element.h"
#include "node.h"

using namespace std;

//! Function which calculates the source for the complete system.
extern int sourceproc(vector<node> Nodes, vector<element>& Elements, sourcebase * sourcetype);
