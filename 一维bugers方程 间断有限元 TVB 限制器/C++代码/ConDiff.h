#ifndef _CONDIFF_H_
#define _CONDIFF_H_

#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <iomanip>
#include "Global.h"
#include "Matrix.h"
#include "DGfunc.h"
#include "MyMath.h"

using namespace std;

////////////////////////////////////////
//   Burgers 方程的DG方法 
////////////////////////////////////////

class CDGWENO
{
public:
	CDGWENO();
	~CDGWENO();
	CDGWENO(double left,double right,double time,unsigned int num);
public:
	void InitialDG();
	void SetBoundary();
public:
	void GetUmod();
	void GetUpm();
public:
	void GetM();
	double TVBLimiter(double a,double b,double c,double M,double h);
	void UseLimiter();
	void UpdateVariables();
public:
	void Getflux();
public:
	double Intergal0(int i);
	double Intergal1(int i);
	double Intergal2(int i);
	double Lh0(int i);
	double Lh1(int i);
	double Lh2(int i);
public:
	void Solvestep();
	void RungeKutta();
	void MainSolver();
	void Output();



private:
	double m_a;
	double m_b;
	double m_time;
	unsigned int m_num;
	double m_t;

	unsigned int End;
	unsigned int Num;

	vector<double> m_M;
	vector<double> m_Hflux;
	
};
#endif