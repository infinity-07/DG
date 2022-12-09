#include "ConDiff.h"

///////////////////////////////////
/////////
///////////////////////////////////
//public:
CDGWENO::CDGWENO()
:m_a(0),m_b(0),m_num(0),m_time(0)
{}
CDGWENO::~CDGWENO()
{}
CDGWENO::CDGWENO(double left, double right, double time, unsigned int num)
{
	m_a = left;
	m_b = right;
	m_time = time;
	m_num = num;
	//CFL = 0.1
	m_t = 0.1 * (m_b - m_a) / (double) num;

	End = Begin + num;
	Num = Begin + End;

	Hstep.resize(Num);
	Center.resize(Num);
	Element.Resize(Num,2);

	U.Resize(Num,3);
	Umod.Resize(Num,2);
	Upm.Resize(Num,2);
	//Hstep
	for(unsigned int i = Begin; i != End; i++)
	{
		Hstep[i] = (m_b - m_a) / (double) num;
	}
	//Element
    Element(Begin,0) = m_a;
	Element(Begin,1) = m_a + Hstep[Begin];
	for(unsigned int i = Begin+1; i != End; i++)
	{
		Element(i,0) = Element(i-1,1);
		Element(i,1) = Element(i,0) + Hstep[i];	
	}
	//center
	for(unsigned int i = Begin; i != End; i++)
	{
		Center[i] = 0.5 * (Element(i,0) + Element(i,1));
	}
	//M
	m_M.resize(Num);
	m_Hflux.resize(Num+1);
}
//public:
void CDGWENO::InitialDG()
{
	double x0,x1,x2;
	double w0,w1,w2;

	w0 = 0.555555555555556;
	w1 = 0.888888888888889;
	w2 = 0.555555555555556;

	for(unsigned int i = Begin; i != End; i++)
	{
		x0 = -0.774596669241483 * 0.5 * (Element(i,1) - Element(i,0)) + Center[i];
		x1 = Center[i];
		x2 = 0.774596669241483 * 0.5 * (Element(i,1) - Element(i,0)) + Center[i];

		U(i,0) = w0 * u0(x0) + w1 * u0(x1) + w2 * u0(x2);
		U(i,0) = 0.5 * (Element(i,1) - Element(i,0)) * U(i,0) / Hstep[i];

		U(i,1) = w0 * u0(x0) * basefun1(i,x0) + w1 * u0(x1) * basefun1(i,x1) + w2 * u0(x2) * basefun1(i,x2);
		U(i,1) = 0.5 * (Element(i,1) - Element(i,0)) * U(i,1) / (Hstep[i] * Hstep[i]);

		U(i,2) = w0 * u0(x0) * basefun2(i,x0) + w1 * u0(x1) * basefun2(i,x1) + w2 * u0(x2) * basefun2(i,x2);
		U(i,2) = 0.5 * (Element(i,1) - Element(i,0)) * U(i,2) / (Hstep[i] * Hstep[i] * Hstep[i]);
	}
}

void CDGWENO::SetBoundary()
{
	//周期边界
	U(Begin-1,0) = U(End-1,0);
	U(Begin-1,1) = U(End-1,1);
	U(Begin-1,2) = U(End-1,2);

	U(End,0) = U(Begin,0);
	U(End,1) = U(Begin,1);
	U(End,2) = U(Begin,2);
}

//public:
void CDGWENO::GetUmod()
{
	for(unsigned int i = Begin; i != End; i++)
	{
		//u~
		Umod(i,0) = 6.0 * U(i,1) + 30.0 * U(i,2);
		//u~~
		Umod(i,1) = 6.0 * U(i,1) - 30.0 * U(i,2);
	}
}
void CDGWENO::GetUpm()
{
	for(unsigned int i = Begin; i != End; i++)
	{
		//u+
		Upm(i,0) = U(i,0) - Umod(i,1);
		//u-
		Upm(i,1) = U(i,0) + Umod(i,0);
	}
}
//public:
void CDGWENO::GetM()
{
	double delta_u1(0),delta_u2(0);

	for(unsigned int i = Begin; i != End; i++)
	{
		delta_u1 = fabs( U(i+1,0) - U(i,0) );
		delta_u2 = fabs( U(i,0) - U(i-1,0) );
		m_M[i] = 2.0/9.0 * (3.0 + 10.0 * M2) * M2 * Hstep[i] * Hstep[i] / (Hstep[i] * Hstep[i] + delta_u1 + delta_u2);
	}
}
double CDGWENO::TVBLimiter(double a, double b, double c, double M, double h)
{
	if(fabs(a) <= M * h * h)
	{
		return a;
	}
	else
		return (Minmod(a,b,c));
}
void CDGWENO::UseLimiter()
{
	double delta_u1(0),delta_u2(0);

	for(unsigned int i = Begin; i != End; i++)
	{
		delta_u1 = U(i+1,0) - U(i,0);
		delta_u2 = U(i,0) - U(i-1,0);
		Umod(i,0) = TVBLimiter(Umod(i,0),delta_u1,delta_u2,m_M[i],Hstep[i]);
		Umod(i,1) = TVBLimiter(Umod(i,1),delta_u1,delta_u2,m_M[i],Hstep[i]);
	}
}
void CDGWENO::UpdateVariables()
{
	for(unsigned int i = Begin; i != End; i++)
	{
		U(i,1) = 1.0/12.0 * (Umod(i,0) + Umod(i,1));
		U(i,2) = 1.0/60.0 * (Umod(i,0) - Umod(i,1));
	}
	for(unsigned int i = Begin; i != End; i++)
	{
		Upm(i,0) = U(i,0) - Umod(i,1);
		Upm(i,1) = U(i,0) + Umod(i,0);
	}
    SetBoundary();
}
//public:
void CDGWENO::Getflux()
{
	double l(0),r(0);

	l = U(Begin-1,0) + 6.0 * U(Begin-1,1) + 30.0 * U(Begin-1,2);
	m_Hflux[Begin] = LLFflux(l,Upm(Begin,0));

	for(unsigned int i = Begin+1; i != End; i++)
	{
		m_Hflux[i] = LLFflux(Upm(i-1,1),Upm(i,0));
	}

	r = U(End,0) - 6.0 * U(End,1) + 30.0 * U(End,2);
	m_Hflux[End] = LLFflux(Upm(End-1,1),r);
}
//public:    
double CDGWENO::Intergal0(int i)
{
	return 0;
}
double CDGWENO::Intergal1(int i)
{
	double x1,x2;
	x1 = -0.447213595 * 0.5 * Hstep[i] + Center[i];
	x2 =  0.447213595 * 0.5 * Hstep[i] + Center[i];

	double r1,r2;
	r1 = u(i,x1);
	r2 = u(i,x2);

	double s;

	s = 0.5 * Hstep[i] * (m_Hflux[i] + 5.0 * fu(r1) + 5.0 * fu(r2) + m_Hflux[i+1]) / 6.0;

	return s;
}
double CDGWENO::Intergal2(int i)
{
	double x1,x2;
	x1 = -0.447213595 * 0.5 * Hstep[i] + Center[i];
	x2 =  0.447213595 * 0.5 * Hstep[i] + Center[i];

	double r1,r2;
	r1 = u(i,x1);
	r2 = u(i,x2);

	double s;

	s = 0.5 * Hstep[i] * ( m_Hflux[i] * (Element(i,0) - Center[i]) + 5.0 * fu(r1) * (x1 - Center[i])
		                   + 5.0 * fu(r2) * (x2 - Center[i]) + m_Hflux[i+1] * (Element(i,1) - Center[i]) ) / 6.0;

	return s;
}
double CDGWENO::Lh0(int i)
{
	return ( -(m_Hflux[i+1] - m_Hflux[i]) / Hstep[i] );
}
double CDGWENO::Lh1(int i)
{
	double s;

	s = - (m_Hflux[i+1] + m_Hflux[i]) / (2.0 * Hstep[i]) + Intergal1(i) / (Hstep[i] * Hstep[i]);

	return s;
}
double CDGWENO::Lh2(int i)
{
	double s;

	s = - (m_Hflux[i+1] - m_Hflux[i]) / (6.0 * Hstep[i]) + 2.0 * Intergal2(i) / (Hstep[i] * Hstep[i] * Hstep[i]);

	return s;
}
//public:
void CDGWENO::Solvestep()
{
	SetBoundary();
	GetUmod();
	GetM();
	UseLimiter();
	UpdateVariables();

	Getflux();
}
void CDGWENO::RungeKutta()
{
	CMatrix<double> U0,U1,U2,Un;
	U0.Resize(Num,3);
	U1.Resize(Num,3);
	U2.Resize(Num,3);
	Un.Resize(Num,3);


	for(unsigned int i = Begin; i != End; i++)
	{
		U0(i,0) = U(i,0);
		U0(i,1) = U(i,1);
		U0(i,2) = U(i,2);
	}

    Solvestep();
	for(unsigned int i = Begin; i != End; i++)
	{
		U1(i,0) = U0(i,0) + m_t * Lh0(i);
		U1(i,1) = U0(i,1) + m_t * Lh1(i);
		U1(i,2) = U0(i,2) + m_t * Lh2(i);
	}
	for(unsigned int i = Begin; i != End; i++)
	{
		U(i,0) = U1(i,0);
		U(i,1) = U1(i,1);
		U(i,2) = U1(i,2);
	}

    Solvestep();
	for(unsigned int i = Begin; i != End; i++)
	{
		U2(i,0) = 0.75 * U0(i,0) + 0.25 * U1(i,0) + 0.25 * m_t * Lh0(i);
		U2(i,1) = 0.75 * U0(i,1) + 0.25 * U1(i,1) + 0.25 * m_t * Lh1(i);
		U2(i,2) = 0.75 * U0(i,2) + 0.25 * U1(i,2) + 0.25 * m_t * Lh2(i);
	}
	for(unsigned int i = Begin; i != End; i++)
	{
		U(i,0) = U2(i,0);
		U(i,1) = U2(i,1);
		U(i,2) = U2(i,2);
	}

    Solvestep();
	for(unsigned int i = Begin; i != End; i++)
	{
		Un(i,0) = 1.0/3.0 * U0(i,0) + 2.0/3.0 * U2(i,0) + 2.0/3.0 * m_t * Lh0(i);
		Un(i,1) = 1.0/3.0 * U0(i,1) + 2.0/3.0 * U2(i,1) + 2.0/3.0 * m_t * Lh1(i);
		Un(i,2) = 1.0/3.0 * U0(i,2) + 2.0/3.0 * U2(i,2) + 2.0/3.0 * m_t * Lh2(i);
	}
	for(unsigned int i = Begin; i != End; i++)
	{
		U(i,0) = Un(i,0);
		U(i,1) = Un(i,1);
		U(i,2) = Un(i,2);
	}
    //rungekutta 输出之前的限制器
	SetBoundary();
	GetUmod();
	GetM();
	UseLimiter();
	UpdateVariables();
}

void CDGWENO::MainSolver()
{
	double now(0);
	InitialDG();
	while (now <= m_time)
	{
		cout<<"time:"<<now<<endl;
		RungeKutta();
		now = now + m_t;
	}
}
void CDGWENO::Output()
{
	ofstream out_U("scalar DG U.txt",ios::out);
	out_U.setf( ios::fixed | ios::showpoint );
	out_U << "X, U0, U1,U2"<<endl;
	for(unsigned int i = Begin; i != End; i++)
	{
		out_U << Center[i] <<" "<< U(i,0)<<" "<< U(i,1)<<" "<<U(i,2) <<endl;
	}
	out_U.close();

	ofstream out_test("scalar DG result.txt",ios::out);
    out_test.setf( ios::fixed | ios::showpoint );
	out_test.precision(12);
	out_test << "X, Y"<<endl;
	for(unsigned int i = Begin; i != End; i++)
	{
		double temp(0);
		temp = u(i,Center[i]);
		out_test << setw(12) << Center[i] <<" "
			     << setw(12) << temp <<endl;
	}
	out_test.close();

}