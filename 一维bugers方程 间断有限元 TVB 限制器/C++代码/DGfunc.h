#ifndef _DGFUNC_H_
#define _DGFUNC_H_
#pragma once

#include <cmath>
#include <cassert>
#include "Global.h"
#include "MyMath.h"
using namespace std;

//基函数
inline double basefun0(int i,double x)
{
	return 1;
}
inline double basefun1(int i,double x)
{
   return (x - Center[i]);
}
inline double basefun2(int i,double x)
{
   return ( (x - Center[i]) * (x - Center[i]) - Hstep[i] * Hstep[i] / 12.0 );
}



//函数 u(x,t)
inline double u(int i,double x)
{
	double a0 = 1;
	double a1 = 12.0 / Hstep[i];
	double a2 = 180.0 / (Hstep[i] * Hstep[i]);

	double u(0);

	u = a0 * U(i,0) + a1 * U(i,1) * (x - Center[i]) + a2 * U(i,2) * ( (x - Center[i]) * (x - Center[i]) - Hstep[i] * Hstep[i] / 12.0 );

	return u;
}

//初值函数
inline double u0(double x)
{
	return sin( - Pi * x);
}
inline double fu(double u)
{
	return ( 0.5 * u * u );
}

//通量函数
inline double LLFflux(double a,double b)
{
	double beta(0);

	beta = max(fabs(a),fabs(b));

	return ( 0.5 * ( fu(a) + fu(b) - beta * (b - a) ) );
}
#endif