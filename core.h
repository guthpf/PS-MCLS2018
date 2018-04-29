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
#include <eigen3/Eigen/Eigen>
#include <list>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include "plotter/plotter.h"

#define ASSERT(x) if(!(x)){std::cout << "ASSERTION ERROR AT LINE " << __LINE__ << " FILE " << __FILE__ << std::endl; exit(1);}

#define DISP(x) std::cout << #x << ": " << x << std::endl;

template <typename T>
inline T sqr(const T & x) {
	return x*x;
}

//Generates a uniform random variable in the interval [0,1)
double unif()
{
	double u;
	u = rand()/(double)RAND_MAX;
	return u;
}

//Generates a Gaussian random variable of zero mean and unit variance
double rgauss()
{
	static bool state = false;
	static double v;
	state = !state;
	if(state)
	{
		double r = sqrt(-2 * log(unif()));
		double t = 2 * M_PI * unif();
		v = r * cos(t);
		return r * sin(t);
	}
	else
	{
		return v;
	}
}

//Input:
//  nsamples: sum(1)
//  sum1: sum(x)
//  sum2: sum(x²)
//Output:
//  string with format "<mean> ± <error>", where
//    <mean> is the average value of x
//    <error> is 3 times the standard deviation of the mean
std::string printpm(double sum1, double sum2, int nsamples)
{
	double e = sum1/nsamples;
	double v = 3*sqrt((sum2/nsamples - e*e)/(nsamples-1.));
	std::ostringstream o;
	o << e << " ± " << v;
	return o.str();
}

template <unsigned n>
using Vec = Eigen::Matrix<double, n, 1>;

template <unsigned n, unsigned m>
using Mat = Eigen::Matrix<double, n, m>;

//Stochastic Objective Function
//requires defining functions to estimate objective function and stochastic gradient
//n: dimensionality
template <unsigned n>
class StochasticOF {
public:
	virtual double objfun(const Vec<n>& x) const = 0;

	virtual Vec<n> gradient(const Vec<n>& x) const = 0;
	
	virtual ~StochasticOF(){}
	
	//If x violates a problem constraint, project it back to the closest feasible value.
	//Does nothing by default.
	virtual void repair(Vec<n>& x) const {
	}

};

//Monte-Carlo Integration Least-Squares type objective function: f(x) = ||Q(x)||²/2
//requires defining function to compute (Q,J) pairs (residual and Jacobian)
//n: dimensionality
//m: height of the Q vectors
template <unsigned n, unsigned m>
class MCILS : public StochasticOF<n> {
public:
	struct QJ {
		Vec<m> Q; //residual
		Mat<m, n> J; //Jacobian
	};

	virtual ~MCILS(){}

	//Estimates (Q, J) pair at x
	//'data' and 'store' may be used by stochastic quasi-Newton type methods, which reuse values of random variables. In this case, use store=true when calling the first time, and store=false when calling for the second. Do not call more than two times for each x.
	virtual QJ jacobian(const Vec<n>& x, std::list<double> * data = NULL, bool store = false) const = 0;

	//Unbiasedly estimates objective function from two (Q,J) pairs
	double objfun(const Vec<n>& x) const {
		QJ qj = jacobian(x);
		QJ qj2 = jacobian(x);
		return .5 * (qj2.Q.transpose() * qj.Q)[0];
	}

	//Unbiasedly estimates objective function from <size> (Q,J) pairs
	//Use size >= 2
	double objfun_block(const Vec<n>& x, int size) const {
		Vec<m> Qsum;
		double bias = 0;
		Qsum.setZero();
		for(int i = 0; i < size; i++) {
			QJ qj = jacobian(x);
			Qsum += qj.Q;
			bias += qj.Q.transpose() * qj.Q;
		}
		return .5 * (Qsum.transpose() * Qsum - bias) / (size * (size - 1.));
	}

	//Unbiasedly estimates JtJ from <size> (Q,J) pairs
	Mat<n, n> JtJ_block(const Vec<n>& x, int size) const {
		Mat<m, n> Jsum;
		Mat<n, n> bias;
		Jsum.setZero();
		bias.setZero();
		for(int i = 0; i < size; i++) {
			QJ qj = jacobian(x);
			Jsum += qj.J;
			bias += qj.J.transpose() * qj.J;
		}
		return (Jsum.transpose() * Jsum - bias) / (size * (size - 1.));
	}

	//Same as JtJ_block(...), but biased (sums bias term to the output instead of subtracting)
	Mat<n, n> JtJ_block_regularized(const Vec<n>& x, int size) const {
		Mat<m, n> Jsum;
		Mat<n, n> bias;
		Jsum.setZero();
		bias.setZero();
		for(int i = 0; i < size; i++) {
			QJ qj = jacobian(x);
			Jsum += qj.J;
			bias += qj.J.transpose() * qj.J;
		}
		return (Jsum.transpose() * Jsum + bias) / (size * (size + 1.));
	}
	
	void objfun_compare(const Vec<n>& x1, const Vec<n>& x2, int block_size, int nsamples) {
		double fsum1 = 0;
		double fsum12 = 0;
		double fsum2 = 0;
		double fsum22 = 0;
		int count = 0;
		std::cout << x1.transpose() << "\t" << x2.transpose();
		for(int i = 0; i < nsamples; i++) {
			double f1 = objfun_block(x1, block_size);
			double f2 = objfun_block(x2, block_size);
			std::cout << i << "\t";
			std::cout.flush();
			fsum1 += f1;
			fsum12 += f1*f1;
			fsum2 += f2;
			fsum22 += f2*f2;
			count++;
			if(count > 1)
				std::cout << printpm(fsum1, fsum12, count) << ",    " << printpm(fsum2, fsum22, count) << ",    " << std::abs((fsum1-fsum2)/count);
			std::cout << std::endl;
		}
	}

	//Unbiasedly estimates gradient from two (Q,J) pairs
	Vec<n> gradient(const Vec<n>& x) const {
		QJ qj = jacobian(x);
		QJ qj2 = jacobian(x);
		return .5 * (qj.J.transpose() * qj2.Q + qj2.J.transpose() * qj.Q);
	}

	//Unbiasedly estimates gradient from <size> (Q,J) pairs
	Vec<n> gradient_block(const Vec<n>& x, int size, std::list<double> * data = NULL, bool store = false) const {
		Vec<m> Qsum;
		Mat<m, n> Jsum;
		Vec<n> bias;
		Qsum.setZero();
		Jsum.setZero();
		bias.setZero();
		for(int i = 0; i < size; i++) {
			QJ qj = jacobian(x, data, store);
			Qsum += qj.Q;
			Jsum += qj.J;
			bias += qj.J.transpose() * qj.Q;
		}
		return (Jsum.transpose() * Qsum - bias) / (size * (size - 1.));
	}

	//tests correctness of the Jacobian estimator comparing to SPSA
	//SPSA estimates J ~ (Q(x+cd) - Q(x-cd)) * d'/2c  (E[dd'] = I)
	//Uses c = c0 * .5^k, and nsamples = 8^k
	//Outputs error between two estimators, for 0<=k<8.
	//Output should converge to zero as k->+∞
	void test_jacobian(const Vec<n> & x, double c0) const {

		//J ~ (Q(x+cd) - Q(x-cd)) * d'/2c  (E[dd'] = I)
		//
		// √2 noise / 2c + (dd'-I) + O(c²)
		// noise² / 2c² + (n-1)I + O(c²)
		//
		//err² ∝ 1/itnum_k.c_k²
		//
		DISP(x.transpose())
		std::cout << std::endl;
		for(int k = 0; k < 8; k++) {
			double c = c0 * std::pow(.5,k);
			double itnum = std::pow(2,3*k);
			Vec<n> xplus, xminus, d;
			Mat<m,n> J, Jcomp;
			J.setZero();
			Jcomp.setZero();
			for(int it = 0; it < itnum; it++) {
				for(int j = 0; j < n; j++) {
					d[j] = (unif() < .5)?1:-1;
					xplus[j] = x[j] + c*d[j];
					xminus[j] = x[j] - c*d[j];
				}
				Jcomp += (jacobian(xplus).Q - jacobian(xminus).Q) * d.transpose();
				J += jacobian(x).J;
			}
			Jcomp /= itnum * 2. * c;
			J /= itnum;
			std::cout << (Jcomp - J).squaredNorm() << std::endl;
		}
		std::cout << std::endl << std::endl;
	}

};

//Input:
//  fun: random number generator
//  data: when NULL, the function simply calls 'fun' and returns its value. When non-NULL, the function will read and pop the values in 'data' if 'store' is true, or generate from 'fun' and write them to 'data' when 'store' is false.
double datastore(std::list<double> * data, bool store, double (*fun)()) {
	if(data) {
		if(store) {
			double u = fun();
			data->push_back(u);
			return u;
		} else {
			double u = data->front();
			data->pop_front();
			return u;
		}
	} else {
		return fun();
	}
}

template <unsigned n>
class Method {
public:
	int k_; //Current iteration
	int t_; //Current time (sum_iterations Nsamples)
	Vec<n> x_; //Current parameters
	std::string caption_; //For plotting
	
	virtual ~Method(){}
	
	virtual void reset(const Vec<n> & x0) {
		x_ = x0;
		t_ = k_ = 0;
	}
	
	//Updates k_, t_, and x_
	virtual void iterate()=0;
	
	//Allocates a clone of the method, with the same parameters. Must delete after using.
	virtual Method<n> * clone()=0;

};

//Clones and compares the methods allocated to f_[*] and shows progress in Plotter
//Input:
//  f_ is an nf-sized array of Method<n>*
//  nf: number of methods
//  fnsamples: number of independent runs per method
//  names: an array of strings of size nf with the names of the methods (for plotting)
//  x0 and xstar are the initial and optimal parameters, respectively
template <unsigned n>
void compareSeveral(Method<n> ** f_, int nf, int fnsamples, Vec<n> x0, Vec<n> xstar) {

	//Init methods
	Method<n> *** f = new Method<n>**[nf];
	for(int i = 0; i < nf; i++) {
		f[i] = new Method<n>*[fnsamples];
		f[i][0] = f_[i];
		f[i][0]->reset(x0);
		for(int j = 1; j < fnsamples; j++) {
			f[i][j] = f[i][0]->clone();
			f[i][j]->reset(x0);
		}
	}

	//init Plotter
	Plotter plotter;
	plotter.xaxisname = "log_{10}(t)";
	plotter.yaxisname = "log_{10}(E[||x_{k}-x*||²])";
	{
		std::ostringstream oss;
		oss << "x_{1} = (" << x0[0];
		for(int i = 1; i < n; i++)
			oss << ", " << x0[i];
		oss << ")";
		plotter.title = oss.str();
	}
	int plotter_nmax = 15000000;
	Plotter::Curve ** c = new Plotter::Curve * [nf];
	for(int i = 0; i < nf; i++) {
		c[i] = plotter.addCurve();
		c[i]->name = f[i][0]->caption_;
		c[i]->x = new double[plotter_nmax];
		c[i]->y = new double[plotter_nmax];
		c[i]->y_u = new double[plotter_nmax];
		c[i]->y_d = new double[plotter_nmax];
		
		//generate curve color automatically
		switch(i) {
		case 0:
			c[i]->col = 0;
			break;
		case 1:
			c[i]->col = 0xFF;
			break;
		case 2:
			c[i]->col = 0xFF0000;
			break;
		case 3:
			c[i]->col = 0x8000;
			break;
		case 4:
			c[i]->col = 0xFF00FF;
			break;
		case 5:
			c[i]->col = 0xC0C0;
			break;
		case 6:
			c[i]->col = 0xF08000;
			break;
		case 7:
			c[i]->col = 0xC0C000;
			break;
		default:
			c[i]->col = rand()%0xFFFFFF;
			break;
		}
		
		//col_u and col_d are the average between col and white
		c[i]->col_u = c[i]->col_d =
			((c[i]->col/0x10000) % 0x100 + 0xFF)/2 * 0x10000 +
			((c[i]->col/0x100) % 0x100 + 0xFF)/2 * 0x100 +
			(c[i]->col % 0x100 + 0xFF)/2;
	}

	//Initial curve values (at t=0)
	for(int i = 0; i < nf; i++) {
		c[i]->x[c[i]->n] = std::log(f[i][0]->t_)/std::log(10);
		double y = 0;
		double err = 0;
		for(int j = 0; j < fnsamples; j++) {
			double pc = (f[i][j]->x_ - xstar).squaredNorm();
			y += pc;
			err += pc*pc;
		}
		y /= fnsamples;
		err = 3*std::sqrt((err/fnsamples - y*y)/(fnsamples - 1));
		if(!(err==err)) //err may be NaN when computing e.g. 3*sqrt(-1e-5)
			err = 0;
		c[i]->y[c[i]->n] = std::log(y)/std::log(10);
		c[i]->y_d[c[i]->n] = std::log(y-err)/std::log(10);
		c[i]->y_u[c[i]->n] = std::log(y+err)/std::log(10);
		c[i]->n++;
	}

	//Main loop
	int ltime;
	while(!plotter.quit) {
		{
			//Select the curve #i that is the most behind in time
			int tmin = f[0][0]->t_;
			int imin = 0;
			for(int i = 1; i < nf; i++)
				if(f[i][0]->t_ < tmin) {
					tmin = f[i][0]->t_;
					imin = i;
				}
			int i = imin;

			//Iterate for this curve #i
			for(int j = 0; j < fnsamples; j++)
				f[i][j]->iterate();

			//Update x, y, etc.
			c[i]->x[c[i]->n] = std::log(f[i][0]->t_)/std::log(10);
			double y = 0;
			double err = 0;
			for(int j = 0; j < fnsamples; j++) {
				double pc = (f[i][j]->x_ - xstar).squaredNorm();
				y += pc;
				err += pc*pc;
			}
			y /= fnsamples;
			err = 3*std::sqrt((err/fnsamples - y*y)/(fnsamples - 1));
			if(!(err==err)) //err may be NaN when computing e.g. 3*sqrt(-1e-5)
				err = 0;
			c[i]->y[c[i]->n] = std::log(y)/std::log(10);
			c[i]->y_d[c[i]->n] = std::log(y-err)/std::log(10);
			c[i]->y_u[c[i]->n] = std::log(y+err)/std::log(10);
			c[i]->n++;

			//avoid overflow
			if(c[i]->n >= plotter_nmax || f[i][0]->t_ > 10500000) {
				plotter.exportSVG();
				plotter.quit = true;
			}
		}
		
		if(plotter.quit)
			break;

		//Draw graph (once per second)
		int ctime = time(NULL);
		if(ctime > ltime)
			plotter.refresh();
		ltime = ctime;
	}
	
	//cleanup	
	for(int i = 0; i < nf; i++) {
		delete[] c[i]->x;
		delete[] c[i]->y;
		delete[] c[i]->y_u;
		delete[] c[i]->y_d;
		for(int j = 1; j < fnsamples; j++)
			delete f[i][j];
		delete[] f[i];
	}
	delete[] f;
	delete[] c;
}

