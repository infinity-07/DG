#include <iostream>
#include "ConDiff.h"
using namespace std;

int main()
{
	CDGWENO A(-1.0,1.0,0.5/Pi,160);
	A.MainSolver();
	A.Output();
	return 0;
}