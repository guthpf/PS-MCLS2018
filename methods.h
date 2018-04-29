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

//A workaround for the fact that one may not pass a double as a template.
//The value is passed as a pair of integers (nummerator and denominator) instead.
#define DOUBLE(c) (c##_num/(double)c##_den)

//Asymptotically constant forget rate
template <int c_num, int c_den>
inline double forget_const(int k) {
	const double c = DOUBLE(c);
	return c/(1-pow(1-c,k));
}

template <int p_num, int p_den>
inline double forget_lin(int k) {
	double p = DOUBLE(p);
	return (1.+p)/(k+p);
}

template <int p_num, int p_den, int fp_num, int fp_den>
inline double forget_pow(int k) {
	double p = DOUBLE(p);
	return std::pow((1.+p)/(k+p), DOUBLE(fp));
}

#undef DOUBLE

template <unsigned n, unsigned m>
class SGD : public Method<n> {
public:
	MCILS<n,m> * problem_;
	Mat<n,n> H_;
	double c_;
	int size_;
	
	SGD(MCILS<n,m> * problem, std::string caption, const Mat<n,n> & H, double c, int size) {
		problem_ = problem;
		this->caption_ = caption;
		H_ = H;
		c_ = c;
		size_ = size;
	};

	Method<n> * clone() {
		return new SGD(problem_, this->caption_, H_, c_, size_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
	}
	
	void iterate() {
		this->k_++;
		this->t_ += size_;
		double a = 1./(this->k_ + c_);
		Vec<n> x0 = this->x_;
		this->x_ -= a*H_*problem_->gradient_block(this->x_, size_);
		if(!(this->x_ == this->x_))
			this->x_ = x0;
		problem_->repair(this->x_);
	}
};

template <unsigned n, unsigned m>
class IP : public Method<n> {
public:
	MCILS<n,m> * problem_;
	Mat<n,n> H_;
	double p_;
	double c_;
	int size0_;
	
	IP(MCILS<n,m> * problem, std::string caption, const Mat<n,n> & H, double p, double c, int size0) {
		problem_ = problem;
		this->caption_ = caption;
		H_ = H;
		p_ = p;
		c_ = c;
		size0_ = size0;
	};

	Method<n> * clone() {
		return new IP(problem_, this->caption_, H_, p_, c_, size0_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
	}
	
	void iterate() {
		this->k_++;
		int size = size0_ * std::pow(this->k_, p_);
		this->t_ += size;
		double stepsize = (1+p_)/(this->k_ + p_ + c_);
		Vec<n> x0 = this->x_;
		this->x_ -= stepsize * H_*problem_->gradient_block(this->x_, size);
		if(!(this->x_ == this->x_))
			this->x_ = x0;
		problem_->repair(this->x_);
	}
};

template <unsigned n, unsigned m>
class AvgIP : public Method<n> {
public:
	MCILS<n,m> * problem_;
	Vec<n> x2_;
	double nsamples_lin_;
	double nsamples_pow_;
	double stepsize_lin_;
	double stepsize_pow_;
	Mat<n, n> H_;
	
	AvgIP(
		MCILS<n,m> * problem, std::string caption,
		Mat<n, n> H,
		double stepsize_lin,
		double stepsize_pow,
		double nsamples_lin,
		double nsamples_pow) {
		
		ASSERT(nsamples_pow >= 0)
		ASSERT(stepsize_pow > 0)
		ASSERT(stepsize_pow < 1)
		ASSERT(stepsize_pow > (1-nsamples_pow)/2.)

		problem_ = problem;
		this->caption_ = caption;
		H_ = H;
		nsamples_lin_ = nsamples_lin;
		nsamples_pow_ = nsamples_pow;
		stepsize_lin_ = stepsize_lin; 
		stepsize_pow_ = stepsize_pow; 

	}

	Method<n> * clone() {
		return new AvgIP(problem_, this->caption_, H_, stepsize_lin_, stepsize_pow_, nsamples_lin_, nsamples_pow_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
		x2_ = x0;
	}
	
	void iterate() {
		this->k_++;
		int nsamples = nsamples_lin_ * std::pow(this->k_, nsamples_pow_);
		double stepsize = stepsize_lin_ * std::pow(this->k_, -stepsize_pow_);
	
		Vec<n> step = stepsize * H_ * problem_->gradient_block(x2_, nsamples);
		Vec<n> x0 = x2_;
		x2_ -= step;
		problem_->repair(x2_);
		if(!(x2_ == x2_))
			x2_ = x0;
		this->x_ += (nsamples / (double)(this->t_ + nsamples)) * (x2_ - this->x_);

		this->t_ += nsamples;
	}
};

template <unsigned n, unsigned m>
class Hybrid : public Method<n> {
public:
	MCILS<n,m> * problem_;
	double s_;
	double c_;
	double (*forget_fun_)(int);
	Mat<n, n> H_;

	Mat<m, n> Jmean;
	Vec<m> Qmean;
	
	Hybrid(MCILS<n,m> * problem, std::string caption, double s, double c, const Mat<n, n> & H, double (*forget_fun)(int)) {
		problem_ = problem;
		this->caption_ = caption;
		s_ = s;
		c_ = c;
		H_ = H;
		forget_fun_ = forget_fun;
	};

	Method<n> * clone() {
		return new Hybrid(problem_, this->caption_, s_, c_, H_, forget_fun_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
		Jmean.setZero();
		Qmean.setZero();
	}
	
	void iterate() {
		this->k_++;
		this->t_++;
		double a = 1./(s_*this->k_ + c_);
		typename MCILS<n,m>::QJ qj = problem_->jacobian(this->x_);
		Vec<n> x0 = this->x_;
		if(this->k_ > 1)
			this->x_ -= H_*((a*.5)*(Jmean.transpose() * qj.Q + qj.J.transpose() * Qmean));
		problem_->repair(this->x_);
		if(!(this->x_ == this->x_))
			this->x_ = x0;
		double forget = forget_fun_(this->k_);
		Jmean = (1-forget) * Jmean + forget * qj.J;
		Qmean = (1-forget) * Qmean + forget * qj.Q;
	}
};

template <unsigned n, unsigned m>
class AvgHybrid : public Method<n> {
public:

	MCILS<n,m> * problem_;
	double lin_;
	double alpha_;
	double (*forget_fun_)(int);

	Mat<n, n> H_;
	Mat<m, n> Jmean;
	Vec<m> Qmean;
	Vec<n> x2;
	
	AvgHybrid(MCILS<n,m> * problem, std::string caption, double lin, double alpha, const Mat<n, n> & H, double (*forget_fun)(int)) {
		problem_ = problem;
		this->caption_ = caption;
		lin_ = lin;
		alpha_ = alpha;
		H_ = H;
		ASSERT(alpha > .5)
		ASSERT(alpha < 1.)
		forget_fun_ = forget_fun;
	};
	
	Method<n> *  clone() {
		return new AvgHybrid(problem_, this->caption_, lin_, alpha_, H_, forget_fun_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
		x2 = x0;
		Jmean.setZero();
		Qmean.setZero();
	}
	
	void iterate() {
		this->k_++;
		this->t_++;
		double a = lin_*std::pow(this->k_, - alpha_);
		typename MCILS<n,m>::QJ qj = problem_->jacobian(x2);
		Vec<n> x0 = x2;
		if(this->k_ > 1)
			x2 -= H_*((a*.5)*(Jmean.transpose() * qj.Q + qj.J.transpose() * Qmean));
		problem_->repair(x2);
		if(!(x2 == x2))
			x2 = x0;
		problem_->repair(x2);
		this->x_ += (x2 - this->x_)/this->k_;
		double forget = forget_fun_(this->k_);
		Jmean += forget * (qj.J - Jmean);
		Qmean += forget * (qj.Q - Qmean);
	}
};

template <unsigned n, unsigned m>
class SGN : public Method<n> {
public:

	MCILS<n,m> * problem_;
	double s_;
	double c_;
	double (*forget_fun_)(int);

	Mat<n, n> H_, R_;
	Mat<m, n> Jmean, Jhess;
	Vec<m> Qmean;
	
	SGN(MCILS<n,m> * problem, std::string caption, double s, double c, double (*forget_fun)(int)) {
		problem_ = problem;
		this->caption_ = caption;
		s_ = s;
		c_ = c;
		forget_fun_ = forget_fun;
	};
	
	Method<n> *  clone() {
		return new SGN(problem_, this->caption_, s_, c_, forget_fun_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
		Jmean.setZero();
		Qmean.setZero();
		H_.setZero();
		R_.setZero();
		Jhess.setZero();
	}
	
	void iterate() {
		this->k_++;
		this->t_++;
		double a = 1./(s_*this->k_ + c_);
		typename MCILS<n,m>::QJ qj = problem_->jacobian(this->x_);
		Vec<n> step = H_*((a*.5)*(Jmean.transpose() * qj.Q + qj.J.transpose() * Qmean));
		Vec<n> xold = this->x_;
		this->x_ -= step;
		problem_->repair(this->x_);
		if(!(this->x_ == this->x_))
			this->x_ = xold;
		double forget = forget_fun_(this->k_);
		Jmean += forget * (qj.J - Jmean);
		Qmean += forget * (qj.Q - Qmean);
		Jhess += qj.J;
		R_ += qj.J.transpose() *  qj.J;
		Mat<m, n> Jtemp = (Jhess - Jmean / forget);
		H_ = ((Jtemp.transpose() * Jtemp + R_)/(sqr(this->k_ - 1./forget) + this->k_)).inverse();
	}
};


template <unsigned n, unsigned m>
class AvgSGN : public Method<n> {
public:

	MCILS<n,m> * problem_;
	double lin_;
	double alpha_;
	double (*forget_fun_)(int);

	Mat<n, n> H_, R_;
	Mat<m, n> Jmean, Jhess;
	Vec<m> Qmean;
	Vec<n> x2;
	
	AvgSGN(MCILS<n,m> * problem, std::string caption, double lin, double alpha, double (*forget_fun)(int)) {
		problem_ = problem;
		this->caption_ = caption;
		lin_ = lin;
		alpha_ = alpha;
		ASSERT(alpha > .5)
		ASSERT(alpha < 1.)
		forget_fun_ = forget_fun;
	};
	
	Method<n> *  clone() {
		return new AvgSGN(problem_, this->caption_, lin_, alpha_, forget_fun_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
		x2 = x0;
		Jmean.setZero();
		Qmean.setZero();
		H_.setZero();
		R_.setZero();
		Jhess.setZero();
	}
	
	void iterate() {
		this->k_++;
		this->t_++;
		double a = lin_*std::pow(this->k_, - alpha_);
		typename MCILS<n,m>::QJ qj = problem_->jacobian(x2);
		Vec<n> step = (.5*a) * (H_ * (Jmean.transpose() * qj.Q + qj.J.transpose() * Qmean));
		if(!(step.transpose() * step < 1./0.))//To avoid the case when H_ is inf or nan (may happen in the initial iterations)
			step.setZero();
		Vec<n> xold = x2;
		x2 -= step;
		problem_->repair(x2);
		this->x_ += (x2 - this->x_)/this->k_;
		double forget = forget_fun_(this->k_);
		Jmean += forget * (qj.J - Jmean);
		Qmean += forget * (qj.Q - Qmean);
		Jhess += qj.J;
		R_ += qj.J.transpose() *  qj.J;
		Mat<m, n> Jtemp = (Jhess - Jmean / forget);
		H_ = ((Jtemp.transpose() * Jtemp + R_)/(sqr(this->k_ - 1./forget) + this->k_)).inverse();
	}
};

//Implements a modified version of averaged SGN that employs two forget rates:
//q_k = k^forget1 for the gradient terms, and q_k = k^forget2 for the Gauss-Newton terms.
template <unsigned n, unsigned m>
class AvgSGNForget : public Method<n> {
public:

	MCILS<n,m> * problem_;
	double lin_;
	double alpha_;
	double forget1_;
	double forget2_;

	Mat<n, n> H_, R_;
	Mat<m, n> qJsum, wJsum, wqJsum;
	Vec<m> qQsum;
	Vec<n> x2;
	double qsum, w2sum, wsum, wqsum;
	
	AvgSGNForget(MCILS<n,m> * problem, std::string caption, double lin, double alpha, double forget1, double forget2) {
		problem_ = problem;
		this->caption_ = caption;
		lin_ = lin;
		alpha_ = alpha;
		ASSERT(alpha > .5)
		ASSERT(alpha < 1.)
		forget1_ = forget1;
		forget2_ = forget2;
		ASSERT(forget1 >= 0)
		ASSERT(forget2 >= 0)
	};
	
	Method<n> *  clone() {
		return new AvgSGNForget(problem_, this->caption_, lin_, alpha_, forget1_, forget2_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
		x2 = x0;
		qJsum.setZero();
		qQsum.setZero();
		H_.setZero();
		R_.setZero();
		qJsum.setZero();
		wJsum.setZero();
		wqJsum.setZero();
		qsum = w2sum = wsum = wqsum = 0;
	}
	
	void iterate() {
		this->k_++;
		this->t_++;
		double a = lin_*std::pow(this->k_, - alpha_);
		typename MCILS<n,m>::QJ qj = problem_->jacobian(x2);
		Vec<n> step;
		if(this->k_ == 1)
			step.setZero();
		else
			step = (.5*a) * (H_ * (qJsum.transpose() * qj.Q + qj.J.transpose() * qQsum)/qsum);
		if(!(step.transpose() * step < 1./0.))//To avoid the case when H_ is inf or nan (may happen in the initial iterations)
			step.setZero();
		Vec<n> xold = x2;
		x2 -= step;
		problem_->repair(x2);
		this->x_ += (x2 - this->x_)/this->k_;
		double q = std::pow(this->k_,forget1_);
		double w = std::pow(this->k_,forget2_);
		qQsum += q * qj.Q;
		qJsum += q * qj.J;
		wqJsum += w * q * qj.J;
		wJsum += w * qj.J;
		R_ += w * w * qj.J.transpose() *  qj.J;
		qsum += q;
		wsum += w;
		wqsum += w*q;
		w2sum += w*w;
		Mat<m, n> Jtemp = wJsum - wqJsum / q;
		H_ = ((Jtemp.transpose() * Jtemp + R_)/(sqr(wsum - wqsum/q) + w2sum)).inverse();
	}
};

//Implements Schraudolph et al.'s method (Schraudolph et al., "A stochastic quasi-Newton method for online convex optimization", Artificial Intelligence and Statistics, 2007)
template <unsigned n, unsigned m>
class Schraudolph : public Method<n> {
public:
	MCILS<n,m> * problem_;
	int size_;
	double lambda_;
	double eps_;
	double s_;
	double c_;

	Mat<n, n> H;
	
	Schraudolph(MCILS<n,m> * problem, std::string caption, int size, double lambda, double c, double s = .1, double eps = 1e-10) {
		problem_ = problem;
		this->caption_ = caption;
		size_ = size;
		lambda_ = lambda;
		eps_ = eps;
		c_ = c;
		s_ = s;
	};

	Method<n> * clone() {
		return new Schraudolph(problem_, this->caption_, size_, lambda_, c_, s_, eps_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
		H.setZero();
		for(int i = 0; i < n; i++)
			H(i, i) = eps_;
	}
	
	void iterate() {
		this->k_++;
		this->t_ += size_;
		std::list<double> q;
		Vec<n> g = problem_->gradient_block(this->x_, size_, &q, true);
		Vec<n> p = -H * g;
		Vec<n> x0 = this->x_;
		this->x_ += p / (s_ * this->k_ + c_);
		problem_->repair(this->x_);
		if(!(this->x_ == this->x_)) {
			this->x_ = x0;
			return;
		}
		Vec<n> svec = this->x_ - x0;
		Vec<n> y = problem_->gradient_block(this->x_, size_, &q, false) - g + lambda_ * svec;
		double rho = 1./ (svec.transpose() * y);
		if(this->k_ == 1) {
			double val = 1. / (rho * y.transpose() * y);
			for(int i = 0; i < n; i++)
				H(i, i) = val;
		}
		H -= rho * (H * y) * svec.transpose();
		H -= rho * svec * (y.transpose() * H);
		H += s_ * rho * svec * svec.transpose();
		if(!(H==H)) {
			double val = 1. / (rho * y.transpose() * y);
			for(int i = 0; i < n; i++)
				H(i, i) = val;
		}
	}
};

template <unsigned n, unsigned m>
class Wang : public Method<n> {
public:
	MCILS<n,m> * problem_;
	int size_;
	double eps_;
	double s_;
	double c_;

	Mat<n, n> H;
	
	Wang(MCILS<n,m> * problem, std::string caption, int size, double c, double s = 1., double eps = 1e-10) {
		problem_ = problem;
		this->caption_ = caption;
		size_ = size;
		eps_ = eps;
		c_ = c;
		s_ = s;
	};

	Method<n> * clone() {
		return new Wang(problem_, this->caption_, size_, c_, s_, eps_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
		H.setZero();
		for(int i = 0; i < n; i++)
			H(i, i) = eps_;
	}
	
	void iterate() {
		this->k_++;
		this->t_ += size_;
		
		//Update x
		std::list<double> q;
		Vec<n> g = problem_->gradient_block(this->x_, size_, &q, true);
		Vec<n> p = -H * g;
		Vec<n> x0 = this->x_;
		this->x_ += p / (s_ * this->k_ + c_);
		problem_->repair(this->x_);
		if(!(this->x_ == this->x_)) {
			this->x_ = x0;
			return;
		}
		
		//Calculate s and y
		Vec<n> svec = this->x_ - x0;
		Vec<n> y = problem_->gradient_block(this->x_, size_, &q, false) - g;
		double rho = 1./ (svec.transpose() * y);
		if(this->k_ == 1) {
			double val = 1. / (rho * y.transpose() * y);
			double invdelta_ = 1./0.;
			if(val < 0 || val > invdelta_)
				val = invdelta_;
			for(int i = 0; i < n; i++)
				H(i, i) = val;
		}
		
		//Regularize
		Vec<n> Bs = H.inverse()*svec;
		double theta1 = svec.transpose()*Bs;
		double theta2 = 1/rho;
		double theta = ((theta2 < .25 * theta1)?(.75 * theta1/(theta1 - theta2)):1);
		Vec<n> y_reg = theta * y + (1-theta) * Bs;
		double rho_reg = 1./ (svec.transpose() * y_reg);

		//Update H
		Mat<n,n> H0 = H;
		H -= rho_reg * (H * y_reg) * svec.transpose();
		H -= rho_reg * svec * (y_reg.transpose() * H);
		H += s_ * rho_reg * svec * svec.transpose();
		if(!(H==H)) {
			double val = 1. / (rho_reg * y_reg.transpose() * y_reg);
			for(int i = 0; i < n; i++)
				H(i, i) = val;
			if(!(H==H))
				H = H0;
		}
	}
};

template <unsigned n, unsigned m>
class Bordes : public Method<n> {
public:
	MCILS<n,m> * problem_;
	int size_;
	double lambda_, c_;
	bool extracond_;

	Mat<n, n> H;
	
	Bordes(MCILS<n,m> * problem, std::string caption, int size, double lambda, double c, bool extracond = true) {
		problem_ = problem;
		this->caption_ = caption;
		size_ = size;
		lambda_ = lambda;
		c_ = c;
		extracond_ = extracond;
	};

	Method<n> * clone() {
		return new Bordes(problem_, this->caption_, size_, lambda_, c_, extracond_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
		H.setZero();
		for(int i = 0; i < n; i++)
			H(i, i) = 1./lambda_;
	}
	
	void iterate() {
		this->k_++;
		this->t_ += size_;
		std::list<double> q;
		Vec<n> g = problem_->gradient_block(this->x_, size_, &q, true);
		Vec<n> x0 = this->x_;
		this->x_ -= (H * g) / (this->k_ + c_);
		problem_->repair(this->x_);
		if(!(this->x_ == this->x_)) {
			this->x_ = x0;
			return;
		}
		Vec<n> svec = this->x_ - x0;
		Vec<n> y = problem_->gradient_block(this->x_, size_, &q, false) - g;
		for(int i = 0; i < n; i++) {
			double parc = svec[i]/y[i];

			double & Hii = H(i, i);
			Hii += 2./(this->k_ + 1) * (parc - Hii);
			if(extracond_ && !(Hii <= 100./lambda_))
				Hii = 100./lambda_;
			if(!(Hii >= .01/lambda_))
				Hii = .01/lambda_;

		}
	}
};

//Implements Wei's method (C. Z. Wei, "Multivariate adaptive stochastic approximation", The Annals of Statistics, 1987)
template <unsigned n, unsigned m>
class Wei : public Method<n> {
public:
	MCILS<n,m> * problem_;
	int size_;
	Mat<n, n> Bsum;
	double lambdasum;
	double k0_;
	double c0_;
	
	Wei(MCILS<n,m> * problem, std::string caption, int size, double k0, double c0) {
		problem_ = problem;
		this->caption_ = caption;
		size_ = size;
		k0_ = k0;
		c0_ = c0;
	};

	Method<n> * clone() {
		return new Wei(problem_, this->caption_, size_, k0_, c0_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
		Bsum.setZero();
		lambdasum = 0;
	}
	
	void iterate() {
		this->k_++;
		this->t_ += 2 * size_ * n;
		double stepsize = 1./(this->k_ + k0_);
		double c = c0_/std::sqrt(this->k_*std::sqrt(std::log(this->k_+1)/std::log(2)));
		double lambda = c*c;

		Vec<n> g;
		g.setZero();
		Vec<n> d;
		Mat<n, n> Bcur;
		d.setZero();
		for(int i = 0; i < n; i++) {
			d[i] = c;
			Vec<n> g1 = problem_->gradient_block(this->x_ + d, size_);
			Vec<n> g2 = problem_->gradient_block(this->x_ - d, size_);
			g += g1 + g2;
			for(int j = 0; j < n; j++)
				Bcur(j, i) = g1[j] - g2[j];
			d[i] = 0;
		}
		Bsum += lambda * (Bcur+Bcur.transpose()) / (4*c);
		lambdasum += lambda;
		g /= 2*n;
		Mat<n, n> H = (Bsum/lambdasum).inverse();
		Vec<n> x0 = this->x_;
		this->x_ -= stepsize * (H * g);
		problem_->repair(this->x_);
		if(!(this->x_ == this->x_))
			this->x_ = x0;
	}
};

template <unsigned n, unsigned m>
class AdaGrad : public Method<n> {
	MCILS<n,m> * problem_;
	double a_;
	int size_;

	Mat<n, n> GGt;
	Mat<n, n> H;
	
public:
	AdaGrad(MCILS<n,m> * problem, std::string caption, double a, int size) {
		problem_ = problem;
		this->caption_ = caption;
		a_ = a;
		size_ = size;
	};

	Method<n> * clone() {
		return new AdaGrad(problem_, this->caption_, a_, size_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
		H.setZero();
		GGt.setZero();
	}
	
	void iterate() {
		this->k_++;
		this->t_ += size_;
		Vec<n> g = problem_->gradient_block(this->x_, size_);
		GGt += g*g.transpose();
		for(int j = 0; j < n; j++)
			H(j,j) = 1./std::sqrt(GGt(j,j));
		Vec<n> x0 = this->x_;
		this->x_ -= a_*H*g;
		problem_->repair(this->x_);
		if(!(this->x_ == this->x_))
			this->x_ = x0;
	}
	
};

//Implements AdaGrad-IP hybrid method
template <unsigned n, unsigned m>
class AdaGradHybrid : public Method<n> {
public:
	MCILS<n,m> * problem_;
	double a_;
	double (*forget_fun_)(int);

	Mat<m, n> Jmean;
	Vec<m> Qmean;
	Mat<n, n> GGt;
	Mat<n, n> H;
	
	AdaGradHybrid(MCILS<n,m> * problem, std::string caption, double a, double (*forget_fun)(int)) {
		problem_ = problem;
		this->caption_ = caption;
		a_ = a;
		forget_fun_ = forget_fun;
	};

	Method<n> * clone() {
		return new AdaGradHybrid(problem_, this->caption_, a_, forget_fun_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
		Jmean.setZero();
		Qmean.setZero();
		GGt.setZero();
		H.setZero();
	}
	
	void iterate() {
		this->k_++;
		this->t_++;
		typename MCILS<n,m>::QJ qj = problem_->jacobian(this->x_);
		Vec<n> x0 = this->x_;
		if(this->k_ > 1) {
			Vec<n> g = .5*(Jmean.transpose() * qj.Q + qj.J.transpose() * Qmean);
			GGt += g*g.transpose();
			for(int j = 0; j < n; j++)
				H(j,j) = 1./std::sqrt(GGt(j,j));
			this->x_ -= a_*(H*g);
		}
		problem_->repair(this->x_);
		if(!(this->x_ == this->x_))
			this->x_ = x0;
		double forget = forget_fun_(this->k_);
		Jmean = (1-forget) * Jmean + forget * qj.J;
		Qmean = (1-forget) * Qmean + forget * qj.Q;
	}
};

//Implements Adam with constant learning rate a
template <unsigned n, unsigned m>
class AdamConst : public Method<n> {

	MCILS<n,m> * problem_;
	double a_, b1_, b2_, c_;
	int size_;

	Vec<n> mv;
	Mat<n, n> v, H;
	
public:
	AdamConst(MCILS<n,m> * problem, std::string caption, int size, double a=.001, double b1=.9, double b2=.999, double c=1e-8) {
		problem_ = problem;
		this->caption_ = caption;
		a_ = a;
		b1_ = b1;
		b2_ = b2;
		c_ = c;
		size_ = size;
	};

	Method<n> * clone() {
		return new AdamConst(problem_, this->caption_, size_, a_, b1_, b2_, c_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
		mv.setZero();
		v.setZero();
		H.setZero();
	}
	
	void iterate() {
		this->k_++;
		this->t_ += size_;
		Vec<n> g = problem_->gradient_block(this->x_, size_);
		mv = b1_ * mv + (1 - b1_) * g;
		v = b2_ * v + (1 - b2_) * g*g.transpose();
		Vec<n> mfix = mv/(1-std::pow(b1_,this->k_));
		Mat<n, n> vfix = v/(1-std::pow(b2_,this->k_));
		for(int j = 0; j < n; j++)
			H(j,j) = 1./(std::sqrt(vfix(j,j))+c_);
		Vec<n> x0 = this->x_;
		this->x_ -= a_*H*mfix;
		problem_->repair(this->x_);
		if(!(this->x_ == this->x_))
			this->x_ = x0;
	}
	
};

//Implements IP-Adam hybrid method with constant learning rate a
template <unsigned n, unsigned m>
class AdamConstHybrid : public Method<n> {

	MCILS<n,m> * problem_;
	double a_, b1_, b2_, c_;
	Mat<m, n> Jmean;
	Vec<m> Qmean;
	double (*forget_fun_)(int);

	Vec<n> mv;
	Mat<n, n> v, H;
	
public:
	AdamConstHybrid(MCILS<n,m> * problem, std::string caption, double (*forget_fun)(int), double a=.001, double b1=.9, double b2=.999, double c=1e-8) {
		problem_ = problem;
		this->caption_ = caption;
		forget_fun_ = forget_fun;
		a_ = a;
		b1_ = b1;
		b2_ = b2;
		c_ = c;
	};

	Method<n> * clone() {
		return new AdamConstHybrid(problem_, this->caption_, forget_fun_, a_, b1_, b2_, c_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
		mv.setZero();
		v.setZero();
		H.setZero();
		Jmean.setZero();
		Qmean.setZero();
	}
	
	void iterate() {
		this->k_++;
		this->t_++;
		typename MCILS<n,m>::QJ qj = problem_->jacobian(this->x_);
		Vec<n> g = .5*(Jmean.transpose() * qj.Q + qj.J.transpose() * Qmean);
		if(this->k_ > 1) {
			mv = b1_ * mv + (1 - b1_) * g;
			v = b2_ * v + (1 - b2_) * g*g.transpose();
			Vec<n> mfix = mv/(1-std::pow(b1_,this->k_-1));
			Mat<n, n> vfix = v/(1-std::pow(b2_,this->k_-1));
			for(int j = 0; j < n; j++)
				H(j,j) = 1./(std::sqrt(vfix(j,j))+c_);
			Vec<n> x0 = this->x_;
			this->x_ -= a_*H*mfix;
			problem_->repair(this->x_);
			if(!(this->x_ == this->x_))
				this->x_ = x0;
		}
		double forget = forget_fun_(this->k_);
		Jmean = (1-forget) * Jmean + forget * qj.J;
		Qmean = (1-forget) * Qmean + forget * qj.Q;
	}
	
};

//Implements Adam with learning rate alpha_k = a/(k+ac)^.5
template <unsigned n, unsigned m>
class AdamSqrt : public Method<n> {

	MCILS<n,m> * problem_;
	double a_, ac_, b1_, b2_, c_;
	int size_;

	Vec<n> mv;
	Mat<n, n> v, H;
	
public:
	AdamSqrt(MCILS<n,m> * problem, std::string caption, int size, double ac, double a=.001, double b1=.9, double b2=.999, double c=1e-8) {
		problem_ = problem;
		this->caption_ = caption;
		a_ = a;
		ac_ = ac;
		b1_ = b1;
		b2_ = b2;
		c_ = c;
		size_ = size;
	};

	Method<n> * clone() {
		return new AdamSqrt(problem_, this->caption_, size_, ac_, a_, b1_, b2_, c_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
		mv.setZero();
		v.setZero();
		H.setZero();
	}
	
	void iterate() {
		this->k_++;
		this->t_ += size_;
		Vec<n> g = problem_->gradient_block(this->x_, size_);
		mv = b1_ * mv + (1 - b1_) * g;
		v = b2_ * v + (1 - b2_) * g*g.transpose();
		Vec<n> mfix = mv/(1-std::pow(b1_,this->k_));
		Mat<n, n> vfix = v/(1-std::pow(b2_,this->k_));
		for(int j = 0; j < n; j++)
			H(j,j) = 1./(std::sqrt(vfix(j,j))+c_);
		Vec<n> x0 = this->x_;
		this->x_ -= a_*H*mfix/sqrt(this->k_+ac_);
		problem_->repair(this->x_);
		if(!(this->x_ == this->x_))
			this->x_ = x0;
	}
	
};

//Implements IP-Adam hybrid method with learning rate alpha_k = a/(k+ac)^.5
template <unsigned n, unsigned m>
class AdamSqrtHybrid : public Method<n> {

	MCILS<n,m> * problem_;
	double a_, ac_, b1_, b2_, c_;
	Mat<m, n> Jmean;
	Vec<m> Qmean;
	double (*forget_fun_)(int);

	Vec<n> mv;
	Mat<n, n> v, H;
	
public:
	AdamSqrtHybrid(MCILS<n,m> * problem, std::string caption, double (*forget_fun)(int), double ac, double a=.001, double b1=.9, double b2=.999, double c=1e-8) {
		problem_ = problem;
		this->caption_ = caption;
		forget_fun_ = forget_fun;
		a_ = a;
		b1_ = b1;
		b2_ = b2;
		c_ = c;
		ac_ = ac;
	};

	Method<n> * clone() {
		return new AdamSqrtHybrid(problem_, this->caption_, forget_fun_, ac_, a_, b1_, b2_, c_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
		mv.setZero();
		v.setZero();
		H.setZero();
		Jmean.setZero();
		Qmean.setZero();
	}
	
	void iterate() {
		this->k_++;
		this->t_++;
		typename MCILS<n,m>::QJ qj = problem_->jacobian(this->x_);
		Vec<n> g = .5*(Jmean.transpose() * qj.Q + qj.J.transpose() * Qmean);
		if(this->k_ > 1) {
			mv = b1_ * mv + (1 - b1_) * g;
			v = b2_ * v + (1 - b2_) * g*g.transpose();
			Vec<n> mfix = mv/(1-std::pow(b1_,this->k_-1));
			Mat<n, n> vfix = v/(1-std::pow(b2_,this->k_-1));
			for(int j = 0; j < n; j++)
				H(j,j) = 1./(std::sqrt(vfix(j,j))+c_);
			Vec<n> x0 = this->x_;
			this->x_ -= a_*H*mfix/sqrt(this->k_+ac_);
			problem_->repair(this->x_);
			if(!(this->x_ == this->x_))
				this->x_ = x0;
		}
		double forget = forget_fun_(this->k_);
		Jmean = (1-forget) * Jmean + forget * qj.J;
		Qmean = (1-forget) * Qmean + forget * qj.Q;
	}
	
};

//Averaged IP-Adam hybrid method with learning rate a/k^alpha
template <unsigned n, unsigned m>
class AdamHybridAvg : public Method<n> {

	MCILS<n,m> * problem_;
	double a_, alpha_, b1_, b2_, c_;
	Mat<m, n> Jmean;
	Vec<m> Qmean;
	double (*forget_fun_)(int);

	Vec<n> mv, x2;
	Mat<n, n> v, H;
	
public:
	AdamHybridAvg(MCILS<n,m> * problem, std::string caption, double (*forget_fun)(int), double alpha, double a=.001, double b1=.9, double b2=.999, double c=1e-8) {
		problem_ = problem;
		this->caption_ = caption;
		forget_fun_ = forget_fun;
		a_ = a;
		b1_ = b1;
		b2_ = b2;
		c_ = c;
		alpha_ = alpha;
	};

	Method<n> * clone() {
		return new AdamHybridAvg(problem_, this->caption_, forget_fun_, alpha_, a_, b1_, b2_, c_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
		x2 = x0;
		mv.setZero();
		v.setZero();
		H.setZero();
		Jmean.setZero();
		Qmean.setZero();
	}
	
	void iterate() {
		this->k_++;
		this->t_++;
		typename MCILS<n,m>::QJ qj = problem_->jacobian(x2);
		Vec<n> g = .5*(Jmean.transpose() * qj.Q + qj.J.transpose() * Qmean);
		if(this->k_ > 1) {
			mv = b1_ * mv + (1 - b1_) * g;
			v = b2_ * v + (1 - b2_) * g*g.transpose();
			Vec<n> mfix = mv/(1-std::pow(b1_,this->k_-1));
			Mat<n, n> vfix = v/(1-std::pow(b2_,this->k_-1));
			for(int j = 0; j < n; j++)
				H(j,j) = 1./(std::sqrt(vfix(j,j))+c_);
			Vec<n> x0 = x2;
			x2 -= a_*H*mfix/pow(this->k_,alpha_);
			problem_->repair(x2);
			if(!(x2 == x2))
				x2 = x0;
		}
		double forget = forget_fun_(this->k_);
		Jmean = (1-forget) * Jmean + forget * qj.J;
		Qmean = (1-forget) * Qmean + forget * qj.Q;
		this->x_ += (x2 - this->x_)/this->k_;
	}
	
};

//Implements Adam method with learning rate alpha_k = a/(k+ac)
template <unsigned n, unsigned m>
class AdamLin : public Method<n> {

	MCILS<n,m> * problem_;
	double a_, ac_, b1_, b2_, c_;
	int size_;

	Vec<n> mv;
	Mat<n, n> v, H;
	
public:
	AdamLin(MCILS<n,m> * problem, std::string caption, int size, double ac, double a=.001, double b1=.9, double b2=.999, double c=1e-8) {
		problem_ = problem;
		this->caption_ = caption;
		a_ = a;
		ac_ = ac;
		b1_ = b1;
		b2_ = b2;
		c_ = c;
		size_ = size;
	};

	Method<n> * clone() {
		return new AdamLin(problem_, this->caption_, size_, ac_, a_, b1_, b2_, c_);
	}

	void reset(const Vec<n> & x0) {
		Method<n>::reset(x0);
		mv.setZero();
		v.setZero();
		H.setZero();
	}
	
	void iterate() {
		this->k_++;
		this->t_ += size_;
		Vec<n> g = problem_->gradient_block(this->x_, size_);
		mv = b1_ * mv + (1 - b1_) * g;
		v = b2_ * v + (1 - b2_) * g*g.transpose();
		Vec<n> mfix = mv/(1-std::pow(b1_,this->k_));
		Mat<n, n> vfix = v/(1-std::pow(b2_,this->k_));
		for(int j = 0; j < n; j++)
			H(j,j) = 1./(std::sqrt(vfix(j,j))+c_);
		Vec<n> x0 = this->x_;
		this->x_ -= a_*H*mfix/(this->k_+ac_);
		problem_->repair(this->x_);
		if(!(this->x_ == this->x_))
			this->x_ = x0;
	}
	
};

