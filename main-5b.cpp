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

	//srand(time(NULL));
	srand(0);

	const unsigned int n = 10;
	const unsigned int m = 25;
	Probl2<n,m> problem;
	Vec<n> x0, xstar = problem.xstar();
	for(int i = 0; i < n; i++)
		x0[i] = xstar[i] + rgauss()/std::sqrt(n);
	problem.repair(x0);

	DISP(x0.transpose())
		
	Method<n> * f[] = {
		new SGN<n,m>(&problem, "SGN", 1., 0., forget_pow<2,1,1,2>),
		new AvgSGN<n,m>(&problem, "aSGN", 1., .75, forget_lin<2,1>),
		new Bordes<n,m>(&problem, "Bordes2009", 20, .1, 1.),
		new Wang<n,m>(&problem, "Wang2015", 20, 0, 1.),
		new AdaGradHybrid<n,m>(&problem, "IP-AdaGrad Hybrid", .1, forget_lin<2,1>),
		new AdamSqrtHybrid<n,m>(&problem, "IP-Adam Hybrid", forget_lin<2,1>, 0, .1),
		new AdamHybridAvg<n,m>(&problem, "Avg IP-Adam Hybrid", forget_lin<2,1>, .75, 1.),
	};


	const int nf = sizeof(f)/sizeof(f[0]);

	compareSeveral<n>(f, nf, 1000, x0, xstar);

	for(int i = 0; i < nf; i++)
		delete f[i];

	return 0;
}


