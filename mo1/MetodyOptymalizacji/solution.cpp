//Do not edit the code below (unless you know what you are doing)

#include"solution.h"

int solution::f_calls = 0;
int solution::g_calls = 0;
int solution::H_calls = 0;

solution::solution(double L)
{
	x = matrix(L);
	g = NAN;
	H = NAN;
	y = NAN;
}

solution::solution(double *A, int n)
{
	x = matrix(A, n);
	g = NAN;
	H = NAN;
	y = NAN;
}

//You can edit the following code

void solution::fit_fun()
{
	matrix Y0(new double[3]{ 5,1,10 }, 3);
	matrix *Y = solve_ode(0, 1, 1000, Y0, x);
	int *w = get_size(Y[1]);
	double max = Y[1](0, 2);
	for (int i = 1; i < w[0]; ++i) 
		if (max < Y[1](i, 2)) max = Y[1](i, 2);
	y = abs(max - 50);
	++f_calls;
}

void solution::grad()
{
	g = NAN;
	++g_calls;
}

void solution::hess()
{
	H = NAN;
	++H_calls;
}

