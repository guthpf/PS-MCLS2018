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

#pragma once
#include "core.h"

//PROBLEM #0:
//Q(x1, x2) = [.5*x1-.1, .5*x2-.3, .5*x1+.5*x2-.9]
//^Q = [U1*x1,U2*x2,U3*x1+U4*x2]
//
//sol: (.5333, .9333)
//
class Probl0 : public MCILS<2,3> {
public:

	QJ jacobian(const Vec<2> & x, std::list<double> * data, bool store) const {
		double u[] = {
			3*(datastore(data, store, unif)-.5)+.5,
			3*(datastore(data, store, unif)-.5)+.5,
			3*(datastore(data, store, unif)-.5)+.5,
			3*(datastore(data, store, unif)-.5)+.5
		};
		QJ qj;
		qj.Q[0] = u[0] * x[0]  -.1;
		qj.Q[1] = u[1] * x[1] -.3;
		qj.Q[2] = u[2] * x[0] + u[3] * x[1] - .9;
		qj.J(0,0) = u[0];
		qj.J(0,1) = 0;
		qj.J(1,0) = 0;
		qj.J(1,1) = u[1];
		qj.J(2,0) = u[2];
		qj.J(2,1) = u[3];
		return qj;
	}
	
	//Minimum x* (theoretical)
	Vec<2> xstar() const {
		Vec<2> x;
		x[0] = 1.6 / 3.;
		x[1] = 2.8 / 3.;
		return x;
	}

	//(JtJ)^{-1} at x* (theoretical)
	Mat<2, 2> Hstar() const {
		Mat<2, 2> H;
		H(0,0) = H(1,1) = 8./3;
		H(0,1) = H(1,0) = -4./3;
		return H;
	}

};


//PROBLEM #1:
// x = [μ, σ]
//
//Q(μ, σ) = 
//[ int_y g(y,μ,σ) q1(y) dy - .5]
//[ int_y g(y,μ,σ) q2(y) dy - .5]
//[ int_y g(y,μ,σ) q3(y) dy - .2]
//
// Derivative:
// dg(y,μ,σ) / [d{μ,σ}.g(y,μ,σ)] =
// d{-(y-μ)²/2σ²} / d{μ,σ} = ((y-μ)/σ², (y-μ)²/σ³ - 1/σ)
//
// q1 = line(y)
// q2 = line(y-1)
// q3 = line(2y-3)
//
// line(y) = (|y|-|y-1|+1)/2
//
// Solution was computed numerically (see xstar()):

double line(double y) {
	return (std::abs(y) - std::abs(y-1) + 1) / 2;
}

double dline(double y) {
	return (y >= 0. && y <= 1.)?1.:0.;
}

class Probl1 : public MCILS<2,3> {
protected:

	void repair(Vec<2>& x) const {
		if(x[0] < 0)
			x[0] = 0;
		if(x[0] > 2)
			x[0] = 2;
		if(x[1] < 1)
			x[1] = 1;
		if(x[1] > 3)
			x[1] = 3;
	}

	QJ jacobian(const Vec<2> & x, std::list<double> * data, bool store) const {
		double mu = x[0];
		double sigma = x[1];
		QJ qj;
		double gvar = datastore(data, store, rgauss);
		double y = mu + sigma * gvar;
		double qf[] = {
			line(y),
			line(y-1),
			line(2*y-3)
		};
		double dqf[] = {
			dline(y),
			dline(y-1),
			dline(2*y-3)*2
		};
		for(int c = 0; c < 3; c++) {
			qj.Q[c] = qf[c];
			qj.J(c,0) = dqf[c];
			qj.J(c,1) = dqf[c] * gvar;
		}
		qj.Q[0] -= .5;
		qj.Q[1] -= .5;
		qj.Q[2] -= .2;
		return qj;
	}

public:
	
	//Minimum x* (numerically estimated)
	Vec<2> xstar() const {
		Vec<2> x;
		x[0] = 0.660877;
		x[1] = 2.28548;
		return x;
	}

};

//PROBLEM #2:
template<unsigned n, unsigned m>
class Probl2 : public MCILS<n,m> {

	Vec<n> xstar_, lb_, ub_;
	Mat<m, n> A_;
	Vec<m> y_;

public:

	static double arcsinh(double x) {
		return std::log(x + std::sqrt(x*x+1));
	}

	Probl2() {
		for(int i = 0; i < m; i++) {
			for(int j = 0; j < n; j++)
				A_(i, j) = std::exp(-.5*(n*m*std::pow(((i+1)/(double)m-(j+1)/(double)n),2)));
			y_[i] = (i+1)*(i+1);
		}
		Vec<n> u = 2*(A_.transpose()*A_).inverse() * A_.transpose() * y_;
		double r = 1.;
		static bool printv = true;
		for(int i = 0; i < n; i++) {
			xstar_[i] = arcsinh(u[i]);
			lb_[i] = arcsinh(u[i]-r);
			ub_[i] = arcsinh(u[i]+r);
			if(printv)
				std::cout << "ub-lb:" << ub_[i] - lb_[i] << std::endl;
		}
		printv = false;
	}

	typename MCILS<n,m>::QJ jacobian(const Vec<n>& x, std::list<double> * data, bool store) const {
		typename MCILS<n,m>::QJ qj;
		Vec<n> Bu, Bdu;
		for(int i = 0; i < n; i++) {
			double B = .5 + .7*(datastore(data, store, unif) - .5);
			double ex = std::exp(x[i]);
			double sinhx = (ex - 1/ex)/2;
			double coshx = (ex + 1/ex)/2;
			Bu[i] = B * sinhx;
			Bdu[i] = B * coshx;
		}
		qj.Q = A_*Bu-y_;
		qj.J = A_*Bdu.asDiagonal();
		return qj;
	}
	
	//Minimum x* (theoretical)
	Vec<n> xstar() {
		return xstar_;
	}
	
	void repair(Vec<n>& x) const {
		double r = 1.;
		for(int i = 0; i < n; i++) {
			if(x[i] < lb_[i])
				x[i] = lb_[i];
			if(x[i] > ub_[i])
				x[i] = ub_[i];
		}
	}
};

