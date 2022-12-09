#include "Global.h"

//grid
vector<double> Hstep;
vector<double> Center;
CMatrix<double> Element;

//u0,u1,u2
CMatrix<double> U;                 

//u~,u~~
CMatrix<double> Umod;              

//u+,u-
CMatrix<double> Upm; 