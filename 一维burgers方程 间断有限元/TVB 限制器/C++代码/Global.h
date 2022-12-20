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
extern vector<double> Hstep;              //ÿ����Ԫ�Ŀռ䲽������
extern vector<double> Center;             //ÿ����Ԫ����������
extern vector<double> real_solve;
extern CMatrix<double> Element;           //ÿ����Ԫ�����ұ߽�

//u0,u1,u2
extern CMatrix<double> U;                 //3��DG����U���ɶ�u0,u1,u2

//u~,u~~
extern CMatrix<double> Umod;              //ÿ����Ԫ�����u~,u~~

//u+,u-
extern CMatrix<double> Upm;               //ÿ����Ԫ�����u+,u-

#endif