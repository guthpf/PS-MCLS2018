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

template <unsigned n>
void solution_precision(Method<n> & method, const Vec<n> x0, int tmax, int nsol, int secsize, double init_rand_rad) {
	Vec<n>* x = new Vec<n>[nsol];
	Method<n> ** m = new Method<n>*[nsol];
	m[0] = &method;
	for(int i = 1; i < nsol; i++)
		m[i] = m[0]->clone();
	for(int i = 0; i < nsol; i++) {
		m[i]->reset(x0);
		for(int j = 0; j < n; j++)
			m[i]->x_[j] += init_rand_rad*rgauss();
	}

	while(method.t_ < tmax) {
		DISP(method.t_)
		for(int i = 0; i < nsol; i++) {
			for(int k = 0; k < secsize; k++)
				m[i]->iterate();
			x[i] = m[i]->x_;
		}

		cout << "x0: " << x0.transpose() << endl;
		Vec<n> xsum;
		Mat<n, n> xsum2;
		xsum.setZero();
		xsum2.setZero();
		for(int i = 0; i < nsol; i++) {
			cout << "x[" << i << "] = "  << x[i].transpose() << endl;
			xsum += x[i];
			xsum2 += x[i] * x[i].transpose();
		}
		xsum /= nsol;
		xsum2 /= nsol;
		Mat<n, n> err2 = (xsum2 - xsum*xsum.transpose())/(nsol-1.);
		cout << endl << "x* = " << endl;
		for(int i = 0; i < n; i++) {
			cout << xsum[i] << " ± " << 3*sqrt(err2(i,i)) << endl;
		}
		cout << endl;
		cout << "variance matrix = " << endl << err2 << endl << endl;
	}

	for(int i = 1; i < nsol; i++)
		delete m[i];
	delete[] m;

	delete[] x;
}

int main() {

	srand(time(NULL));

	const unsigned int n = 2;
	const unsigned int m = 3;
	Probl1 problem;
	Vec<n> x0;
	x0.setZero();
	x0[0] = 1;
	x0[1] = 1;
	Mat<n, n> H = problem.JtJ_block(x0, 1000).inverse();
	AvgIP<n,m> method(&problem, "aSGD", H, 1., .66, 10., 0.);
	solution_precision<n>(method, x0, 1e9, 8, 1000, 0);

	return 0;
}

/*

x0: 1 1
x[0] = 0.660888  2.28403
x[1] = 0.660434  2.28789
x[2] = 0.659873  2.28704
x[3] = 0.660924  2.28522
x[4] = 0.661463  2.28439
x[5] = 0.661087   2.2867
x[6] = 0.660774  2.28517
x[7] = 0.66157 2.28341

x* = 
0.660877 ± 0.000578703
2.28548 ± 0.00167635

variance matrix = 
3.72107e-08  -7.943e-08
 -7.943e-08 3.12237e-07

method.t_: 40000000

*/

