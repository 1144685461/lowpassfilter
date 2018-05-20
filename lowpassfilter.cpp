/*
 *  Copyright(C) 2018 I-CON. All rights reserved.
 */

/*
 *  lowpassfilter.cpp
 *
 *  Original Author: zhuyongsheng, 2018/05/10
 *
 *  Description:
 *  Low Pass Filter Input.
 *
 *  FIR:
 *     FIR  is desiged by windows Function.
 *     h(n) = w(n)hd(n)
 *     h(n) = (wc/PI)*sin(wc(n-(N-1)/2))/(wc(n-(N-1)/2)), 0 <= n <= N-1.
 *
 *  IIR:
 *    G(s) = (kn*s + wn)/(s+ wn) or
 *    G(s) = (s^2 + s*wn/(A*Q) + wn^2)/(s^2 + wn/Q + wn^2)
 *
 *    History:
 *
 *    V1.0    zhuyongsheng    2018/05/10
 *            Create
 */

#include "lowpassfilter.h"
#include <glog/logging.h>

LowPassFilter* LowPassFilterFactory::getLowPassFilter(enum LowPassType type, int N, double t, double wc)
{
	switch(type) {
	case LOW_PASS_FILTER_FIR:
		return new FirLowPassFilter(N, wc, t);
	case LOW_PASS_FILTER_IIR:
		return new IIRLowPassFilter(t, wc, N);
	default:
		return new LowPassFilter(N, t, wc);
	}
}

LowPassFilter::LowPassFilter(int N, double t, double wc):N_(N+1), T_(t),wc_(wc)
{
	input_.resize(N_);
}

void LowPassFilter::saveInputData(double &value)
{
	for (int i = 0; i < N_ - 1; ++i) {
		input_[i] = input_[i+1];
	}

	if (N_-1 >= 0) {
		input_[N_-1] = value;
	}
}

// FIR
FirLowPassFilter::FirLowPassFilter(int N, double wc, double t, enum WinType type):LowPassFilter(N, wc, t),type_(type)
{
	wc_ = wc*t;
	coefficient_.resize(N_);
	setparam();
	//print();
}

void FirLowPassFilter::setparam()
{
	double n;

	for (int i = 0; i < (N_)/2; ++i) {
		n = i - (N_-1)/2.0;
		coefficient_[i] = sin(wc_*n)/(PI*n)*window(i);
		coefficient_[N_- 1 - i] = coefficient_[i];
	}

	if ((N_-1)%2 == 0)
		coefficient_[(N_-1)/2] = wc_/PI;
}

double FirLowPassFilter::filter(double &value)
{
	saveInputData(value);

	double output = 0.0;

	for (int i = 0; i < N_; ++i) {
		output += coefficient_[i]*input_[N_-1-i];
	}

	return output;
}

double FirLowPassFilter::window(int index)
{
	int w  = 1.0;

	switch(type_) {
	case FIR_ORTHOGON:
		w = 1.0;
		break;
	case FIR_TRIANGLE:
		w = 1.0 - fabs(1.0-2*index*(N_-1.0));
		break;
	case FIR_HANNING:
		w = 0.5*(1 - cos(2*PI*index/(N_-1)));
		break;
	case FIR_HAMMING:
		w = 0.54 - 0.46*cos(2*PI*index/(N_-1));
		break;
	case FIR_BLOCKMAN:
		w = 0.42 - 0.5*cos(2*PI*index/(N_-1)) + 0.08*cos(4*PI*index/(N_-1));
	default:
		DCHECK(false);
		break;
	}
	return w;
}

void FirLowPassFilter::print()
{
	for (int i = 0; i < N_; ++i) {
		fprintf(stderr, " coefficient_[%d] = %f\n", i, coefficient_[i]);
	}
}

IIRLowPassFilter::IIRLowPassFilter(double t, double wc, int N):LowPassFilter(N, t, wc)
{
	order_ = static_cast<enum Order>(N);

	DCHECK(T_ != 0);
	c_ = 2/T_;

	output_.resize(N_);
	sA_.resize(N_);
	sB_.resize(N_);
	zA_.resize(N_);
	zB_.resize(N_);

	doExtractHSCoef();
	setparam();
	print();
}

void IIRLowPassFilter::doExtractHSCoef()
{
	switch(order_) {
	case FIRST_ORDER:
		sA_[0] = wc_;
		sB_[0] = wc_;
		sB_[1] = 1;
		break;
	case SECOND_ORDER:
		sA_[0] = wc_*wc_;
		sA_[1] = wc_/(A*Q);
		sA_[2] = 1;

		sB_[0] = wc_*wc_;
		sB_[1] = wc_/Q;
		sB_[2] = 1;
		break;
	case THIRD_ORDER:
	default:
		DCHECK(false);
		break;
	}
}

void IIRLowPassFilter::setparam()
{
	switch(order_) {
	case FIRST_ORDER:
		r_     = sB_[0] + sB_[1]*c_;
		zA_[0] = (sA_[0] + sA_[1]*c_)/r_;
		zA_[1] = (sA_[0] - sA_[1]*c_)/r_;
		zB_[0] = 1;
		zB_[1] = (sB_[0] - sB_[1]*c_)/r_;
		break;
	case SECOND_ORDER:
		r_     = sB_[0] + sB_[1]*c_ + sB_[2]*c_*c_;
		zA_[0] = (sA_[0] + sA_[1]*c_ + sA_[2]*c_*c_)/r_;
		zA_[1] = (2*sA_[0] - 2*sA_[2]*c_*c_)/r_;
		zA_[2] = (sA_[0] - sA_[1]*c_ + sA_[2]*c_*c_)/r_;
		zB_[0] = 1;
		zB_[1] = (2*sB_[0] - 2*sB_[2]*c_*c_)/r_;
		zB_[2] = (sB_[0] - sB_[1]*c_ + sB_[2]*c_*c_)/r_;
		break;
	case THIRD_ORDER:
#if 0
		r_     = sB_[0] + sB_[1]*c_ + sB_[2]*c_*c_ + sB_[3]*c_*c_*c_;
		zA_[0] = (sA_[0] + sA_[1]*c_ + sA_[2]*c_*c_ + sA_[3]*c_*c_*c_)/r_;
		zA_[1] = (3*sA_[0] + sA_[1]*c_ - sA_[2]*c_*c_ - 3*sA_[3]*c_*c_*c_)/r_;
		zA_[2] = (3*sA_[0] - sA_[1]*c_ + sA_[2]*c_*c_ + 3*sA_[3]*c_*c_*c_)/r_;
		zA_[3] = (sA_[0] - sA_[1]*c_ + sA_[2]*c_*c_ - sA_[3]*c_*c_*c_)/r_;

		zB_[0] = 1;
		zB_[1] = (3*sB_[0] + sB_[1]*c_ - sB_[2]*c_*c_ - 3*sB_[3]*c_*c_*c_)/r_;
		zB_[2] = (3*sB_[0] - sB_[1]*c_ - sB_[2]*c_*c_ + 3*sB_[3]*c_*c_*c_)/r_;
		zB_[3] = (sB_[0] - sB_[1]*c_ + sB_[2]*c_*c_ - sB_[3]*c_*c_*c_)/r_;
#endif
	default:
		DCHECK(false);
		break;
	}
}

void IIRLowPassFilter::outputDataMoveOneStep()
{
	for (int i = 0; i < N_ - 1; ++i) {
		output_[i] = output_[i+1];
	}
}

double IIRLowPassFilter::filter(double &value)
{
	saveInputData(value);
	outputDataMoveOneStep();

	output_[N_-1] = 0.0;

	for (int i = 0; i < N_; ++i) {
		output_[N_-1] += zA_[i]*input_[N_-1-i];
	}

	for (int i = 1; i < N_; ++i) {
		output_[N_-1] -= zB_[i]*output_[N_-1-i];
	}
	return output_[N_-1];
}

void IIRLowPassFilter::print()
{
	for (int i = 0; i < N_; ++i) {
		fprintf(stderr, " zA_[%d] = %f, zB_[%d] = %f\n", i, zA_[i], i, zB_[i]);
	}
}
