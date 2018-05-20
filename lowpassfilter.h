#ifndef LOW_PASS_FILTER_H
#define LOW_PASS_FILTER_H

#include <vector>
#include <math.h>

#define  PI  3.1415926

class LowPassFilter {
public:
	LowPassFilter() {}
	LowPassFilter(int N, double t, double wc);

	~LowPassFilter() {}

	virtual void   setparam() {}
	virtual double filter(double &value) { return value; }

protected:
	void saveInputData(double &value);

protected:
	int    N_{0};
	double T_{0.0};
	double wc_{0.0};

	std::vector<double> input_;
	std::vector<double> output_;

private:
	//
};

class FirLowPassFilter : public LowPassFilter {
public:
	enum WinType {
		FIR_ORTHOGON,
		FIR_TRIANGLE,
		FIR_HANNING,
		FIR_HAMMING,
		FIR_BLOCKMAN,
	};

public:
	FirLowPassFilter() {}
	~FirLowPassFilter() {}

	FirLowPassFilter(int N, double wc, double t, enum WinType type = FIR_ORTHOGON);

	virtual void   setparam() override;
	virtual double filter(double &value) override;

private:
	double sinc(double x);
	double window(int index);
	void  print();

private:
	enum WinType type_{FIR_ORTHOGON};

	std::vector<double> coefficient_;
};

class IIRLowPassFilter : public LowPassFilter {
public:
	IIRLowPassFilter() {}
	~IIRLowPassFilter() {}

	IIRLowPassFilter(double t, double wc, int N = 2);

	virtual void   setparam() override;
	virtual double filter(double &value) override;

private:
	void doExtractHSCoef();
	void saveOutputData(double &value);
	void outputDataMoveOneStep();
	void print();

private:
	double c_;
	double r_;

	std::vector<double> sA_;
	std::vector<double> sB_;

	std::vector<double> zA_;
	std::vector<double> zB_;

private:
	enum Order {
		FIRST_ORDER = 1,
		SECOND_ORDER,
		THIRD_ORDER
	};

	enum Order order_{SECOND_ORDER};

private:
	static constexpr int    A = 30;
	static constexpr double Q = 0.62;
};

class LowPassFilterFactory {
public:
	enum LowPassType {
		LOW_PASS_FILTER_UNKNOWN,
		LOW_PASS_FILTER_FIR,
		LOW_PASS_FILTER_IIR
	};

	static LowPassFilter* getLowPassFilter(enum LowPassType type, int N, double t, double wc);
};

inline bool isFirLowpassFilterTapsValid(int order, double wc, double fs)
{
	return (order == 2) ? true:false;
}

#endif

