#ifndef _GLOBAL_H_
#define _GLOBAL_H_
#pragma once

#include <vector>
#include "Matrix.h"

using namespace std;

///////////////////////////////////
///////////
///////////////////////////////////

//const
const unsigned int Begin = 1;
const double M2 = 9.0;


//grid:
extern vector<double> Hstep;              //每个单元的空间步长向量
extern vector<double> Center;             //每个单元的中心坐标
extern vector<double> real_solve;
extern CMatrix<double> Element;           //每个单元的左右边界

//u0,u1,u2
extern CMatrix<double> U;                 //3阶DG方法U自由度u0,u1,u2

//u~,u~~
extern CMatrix<double> Umod;              //每个单元区间的u~,u~~

//u+,u-
extern CMatrix<double> Upm;               //每个单元区间的u+,u-

#endif