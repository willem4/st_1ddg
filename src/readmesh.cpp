#include<iostream>
#include<fstream>
#include<string>
#include<vector>

#include "readmesh.h"
//#include "count.h"
#include "element.h"
#include "node.h"

using namespace std;
int readmesh(string meshname, vector <element> &Elements, vector <node> &Nodes, int &bcl, int &bcr)
{  
  int Nn;
  int Ne; 
  int nodeint1,nodeint2;
  int i,j,val, ElemTemp, Nb, Temp;
  double Version, NodeTemp, length;
  string a;
  ifstream InFile(meshname.c_str(), ios::in);
  InFile >> a;   
  InFile >> Version;
  InFile >> val;
  if (val != 0) //This means the mesh is one dimensional.
    {
      cout<<"The program only works for one dimensional meshes" << "\n";
    }
  //Nodal data 
  InFile >> Nn; //Total number of nodes.
  //	vector<double> Node(Nn+1);
  for(i=1;i<=Nn;i++) //Describe the x-values of the nodes. 
    {   
      InFile >> NodeTemp; 
      Nodes.push_back(node(NodeTemp));// >> Node.at(i);
    }  

  //Element data
  InFile >> Ne; //Total number of elements.	
  for(i=1;i<=Ne;i++) //Loop over elements. 
    {
      InFile >> nodeint1;
      InFile >> nodeint2;
      length = Nodes[nodeint2-1].x - Nodes[nodeint1-1].x; 
      Elements.push_back(element(nodeint1-1,nodeint2-1,length));
 	  //Which nodes make up the Element.
     }
  //Boundary data.      
  InFile >> Nb;
  //Boundary conditions.
  InFile >> Temp; 
  InFile >> bcl;
  InFile >> Temp; 
  InFile >> bcr;
  Temp = 10*bcl+bcr;
  if (Temp != 11 and Temp != 22 and Temp != 23 and Temp != 32 and Temp != 33 and Temp != 44)
    {
      cout << "[readmesh:] Invalid Boundary conditions \n";
    }
  InFile.close();
  return 0;
}
  
  
