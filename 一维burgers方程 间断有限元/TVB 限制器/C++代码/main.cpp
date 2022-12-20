#include <iostream>
#include "ConDiff.h"
using namespace std;

int main()
{
	CDGWENO A(-1.0,1.0,0.5/Pi,20);
	A.MainSolver();
	A.Output();
	A.get_real_solve();
	
	return 0;
}