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

	const unsigned int n = 4;
	const unsigned int m = 5;
	Probl2<n,m> problem;
	Vec<n> x0, xstar = problem.xstar();
	for(int i = 0; i < n; i++)
		x0[i] = xstar[i] + rgauss()/std::sqrt(n);
	problem.repair(x0);

	DISP(x0.transpose())

	Mat<n, n> H = problem.JtJ_block(xstar, 100000).inverse();
	Method<n> * f[] = {
		new SGD<n,m>(&problem, "SGD (N=10)", H, 0., 10),
		new AvgIP<n,m>(&problem, "aSGD (α=3/4, N=10)", H, 1., .75, 10., 0.), //AvgIP with q=0 is equivalent to SGD
		new IP<n,m>(&problem, "IP (N_{k}=2k)", H, 1, 0, 2),
		new AvgIP<n,m>(&problem, "aIP (α=3/4, N_{k}=2k)", H, 1., .75, 2., 1.),
		new Hybrid<n,m>(&problem, "Hybrid (ζ_{k}~√\\{3/k\\})", 1., 0, H, forget_pow<2,1,1,2>),
		new AvgHybrid<n,m>(&problem, "aHybrid (α=3/4, ζ_{k}~3/k)", 1., .75, H, forget_lin<2,1>),
	};

	const int nf = sizeof(f)/sizeof(f[0]);

	compareSeveral<n>(f, nf, 1000, x0, xstar);

	for(int i = 0; i < nf; i++)
		delete f[i];

	return 0;
}

