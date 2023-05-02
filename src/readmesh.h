#ifndef READMESH_H
#define READMESH_H

#include "element.h"
#include "node.h"

#include <vector>
#include <string>

using namespace std;
//! Function for reading the meshfile.
extern int readmesh(string meshname, vector <element> &Elements, vector <node> &Nodes, int &bcl, int &bcr);

#endif
