#include<random>
#include<fstream>
#include <string>
#include"opt_alg.h"
#include"ode_solver.h"

int main()
{
	try
	{
		double x0, a = 1e-4, b = 1e-2, epsilon = 1e-6, d = 1e-5, alfa = 1.5;
		int Nmax = 1000;
		random_device rd;
		string name ;
		double *mac = new double[2]{ 1,100 };
		alfa = 1.15;

		for (int i = 0; i < 3; i++)
		{
			name= "op_1_"+to_string(i);
			ofstream out(name+".csv");
			out << "lp,x0,liczba wywolan,a,b,x*,y*,liczba wywolan,x*,y*,liczba wywolan \n";
			alfa += 0.23*i;
			for (int i = 0; i < 5; i++)
				{
				x0 = (b - a)*rd() / rd.max() + a;
				solution::f_calls = 0;
				double *p = expansion(x0, d, alfa, Nmax);
				out <<i+1<<","<< x0 << ","  << solution::f_calls << "," << p[0] << "," << p[1] ;
				
				solution::f_calls = 0;
				solution opt_f = fib(p[0], p[1], epsilon);
				out << "," << opt_f.x << "," << opt_f.y << "," << solution::f_calls ;

				solution::f_calls = 0;
				solution opt_l= lag(p[0], p[1], epsilon, Nmax);
				out << "," << opt_l.x << "," << opt_l.y << "," << solution::f_calls << "," << "\n";

			}
			out << "przyklad bez ekspansji: , , ,1,100,";
			solution::f_calls = 0;
			solution opt_f = fib(mac[0], mac[1], epsilon);
			out <<opt_f.x<<","<<opt_f.y<<","<< solution::f_calls << ",";
			solution::f_calls = 0;
			solution opt_l = lag(mac[0], mac[1], epsilon, Nmax);
			out << opt_l.x << "," << opt_l.y << "," << solution::f_calls << ",";
			cout << endl << "finished" << endl<<name<<endl;
			out.close();
		}
		
		//solution::f_calls = 0;
		//solution opt_f = fib(a, b, epsilon);
		//cout<<endl << endl << "," << opt_f.x << "," << opt_f.y << "," << solution::f_calls << endl;

		//solution::f_calls = 0;
		//solution opt_l = lag(a, b, epsilon, Nmax);
		//cout << endl << "," << opt_l.x << "," << opt_l.y << "," << solution::f_calls << "," << "\n" << endl;
	}
	catch (char * EX_INFO)
	{
		cout << EX_INFO << endl;
	}
	system("pause");
	return 0;
}
