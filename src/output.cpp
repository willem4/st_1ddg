#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "element.h"
#include "node.h" 
#include "state.h"

using namespace std;

int output(vector<node> Nodes, vector<element> Elements, string filename)
{
  int i,j;
  ofstream OutFile(filename.c_str(), ios::out);
  OutFile << element::systemsize1 << "\t" << element::systemsize2+1 << "\n";
  OutFile << 2*Elements.size() << "\n";
  for (j=0; j<element::systemsize1; j++){
    for(i=0; i<Elements.size(); i++){
      OutFile << setprecision(15) << Nodes[Elements[i].Node1].x << "\t" << setprecision(15) <<Elements[i].U[j].get(-1) << "\n";
      OutFile << setprecision(15) << Nodes[Elements[i].Node2].x << "\t" << setprecision(15) <<Elements[i].U[j].get(1) << "\n";
    }
  }
  for (j=0; j<element::systemsize2; j++){
    for(i=0; i<Elements.size(); i++){
      OutFile << setprecision(15) << Nodes[Elements[i].Node1].x << "\t" << setprecision(15) << Elements[i].Q[j].get(-1) << "\n";
      OutFile << setprecision(15) << Nodes[Elements[i].Node2].x << "\t" << setprecision(15) << Elements[i].Q[j].get(1) << "\n";
    }
  }
  //Krividonova.
  for (j=0; j<1; j++){
    for(i=0; i<Elements.size(); i++){
      OutFile << setprecision(15) << Nodes[Elements[i].Node1].x << "\t" << setprecision(15) <<Elements[i].Kriv << "\n";
      OutFile << setprecision(15) << Nodes[Elements[i].Node2].x << "\t" << setprecision(15) <<Elements[i].Kriv << "\n";
    }
  }
  OutFile.close();
}
