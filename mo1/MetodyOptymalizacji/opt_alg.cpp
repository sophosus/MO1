#include"opt_alg.h"
#include<fstream>
double *expansion(double x0, double d, double alfa, double Nmax)
{
	double *p = new double[2];
	solution X0(x0), X1(x0 + d);
	X0.fit_fun();
	X1.fit_fun();

	if (X0.y == X1.y) {
		p[0] = X0.x(0);
		p[1] = X1.x(0);
		return p;
	}
	if (X0.y < X1.y) {
		d *= -1;
		X1.x = X0.x + d;
		X1.fit_fun();
		if (X0.y <= X1.y) {
			p[0] = X1.x(0);
			p[1] = X0.x(0) - d;
		}
	}

	solution X2;
	int i = 1;
	while (true) {
		X2.x = x0 + pow(alfa, i)*d;
		X2.fit_fun();
		if (X2.y > X1.y || solution::f_calls > Nmax) break;
		X0 = X1;
		X1 = X2;
		++i;
	}
	d > 0 ? p[0] = X0.x(0), p[1] = X2.x(0) : (p[0] = X2.x(0), p[1] = X0.x(0));
	return p;
}

solution fib(double a, double b, double epsilon) {
	int n = static_cast<int>(ceil(log2(sqrt(5)*(b - a) / epsilon) / log2((1 + sqrt(5)) / 2)));
	int *F = new int[n] {1, 1};
	for (int i = 2; i < n; ++i) {
		F[i] = F[i - 2] + F[i - 1];
	}

	solution A(a), B(b), C, D;
	
	C.x = B.x - 1.0*F[n - 2] / F[n - 1] * (B.x - A.x);
	D.x = A.x + B.x - C.x;
	C.fit_fun();
	D.fit_fun();

	for (int i = 0; i <= n - 3; ++i) {
		//cout << endl << B.x - A.x;
		if (C.y < D.y)
			B = D;
		else
		{
			A = C;
		}
		C.x = B.x - 1.0*F[n - i - 2] / F[n - i - 1] * (B.x - A.x);
		D.x = A.x + B.x - C.x;
		C.fit_fun();
		D.fit_fun();
		
	}
	
	//ofstream plik("out2.csv");
	//plik << (B.x - A.x) << ',';
	
	return C;
}

solution lag(double a, double b, double epsilon, int Nmax ) {
	solution A(a), B(b), C, D;
	C.x = 0.5*(A.x + B.x);
	A.fit_fun();
	B.fit_fun();
	C.fit_fun();

	int i = 0;
	while (true)
	{
		matrix As(new double*[3]{ new double[3]{A.x(0)*A.x(0),A.x(0),1},
		new double[3]{B.x(0)*B.x(0),B.x(0),1},
		new double[3]{C.x(0)*C.x(0),C.x(0),1} }, 3, 3);

		matrix bs(new double[3]{ A.y(0),B.y(0),C.y(0) }, 3);
		matrix xs(3, 1);

		try {
			xs = inv(As)*bs;
		}
		catch (char*) {
			C = NAN;
			return C;
		}
		if (xs(0) <= 0) {
			C = NAN;
			return C;
		}
		D.x = -xs(1) / (2 * xs(0));
		if (D.x<A.x || D.x>B.x) {
			C = NAN;
			return C;
		}
		D.fit_fun();
		if (D.x < C.x)
			if (D.y < C.y)
			{
				B = C;
				C = D;
			}
			else
				A = D;
		else
			if (D.y < C.y)
			{
				A = C;
				C = D;
			}
			else
				B = D;
		//ofstream plik("out1.csv");
		//plik << B.x - A.x << ',';
		//cout << B.x - A.x << endl;
		if (((B.x - A.x) < epsilon) || solution::f_calls > Nmax) return C;

	}
}