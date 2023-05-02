#include <string>

#include "enumtypes.h"

//! Function for reading the test. 
extern int readtest(testtypedescr &testtype, timeintmethod &timeint, string &meshname, double &CFL, double &it, double &et, int &deg);
//! Function for writing the test (not used any longer). 
extern int writetest(testtypedescr &testtype, timeintmethod &timeint, string &meshname, double &CFL, double &it, double &et, int &deg);
//! Function for reading the test with filename. 
extern int readtest2(testtypedescr &testtype, timeintmethod &timeint, string &meshname, double &CFL, double &it, double &et, int &deg, string testname);
