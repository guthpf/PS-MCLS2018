//Copyright © 2017-2018 Gustavo Thebit Pfeiffer / Y. Sato laboratory 
/*
    This file is part of PS-MCLS2018.

    PS-MCLS2018 is free software: you can redistribute it and/or modify it
    under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    PS-MCLS2018 is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
    License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with PS-MCLS2018. If not, see <http://www.gnu.org/licenses/>.
*/

#include "methods.h"
#include "problems.h"

using namespace std;

int main() {

	srand(time(NULL));
	//srand(0);

	const unsigned int n = 2;
	const unsigned int m = 3;
	Probl1 problem;
	Vec<n> x0, xstar;
	xstar = problem.xstar();
	x0[0] = 2;
	x0[1] = 1;

	DISP(x0.transpose())
	
	Method<n> * f[] = {
		new AdamSqrt<n,m>(&problem, "Adam, N=2", 2, 0, .1),
		new AdamSqrt<n,m>(&problem, "Adam, N=10", 10, 0, .1),
		new AdamSqrt<n,m>(&problem, "Adam, N=30", 30, 0, .1),
		new AdamSqrtHybrid<n,m>(&problem, "IP-Adam Hybrid", forget_lin<2,1>, 0, .1),
	};

	const int nf = sizeof(f)/sizeof(f[0]);

	compareSeveral<n>(f, nf, 1000, x0, xstar);

	for(int i = 0; i < nf; i++)
		delete f[i];

	return 0;
}


