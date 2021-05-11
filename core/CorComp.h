#ifndef CORCOMP_h
#define CORCOMP_h

#include <vector>

class CorComp{
	public:
	CorComp(std::vector<double> v1, std::vector<double> v2);
	double computePearsonCorrelation();
	double computeSpearmanCorrelation();
	const unsigned int getObservationCount();
	std::pair<double, double> getFisherZ(double cor);
	double cdf(double x);
	double getPvalue(double cor);
	double computeMean(const std::vector<double> v);
	double computeVariance(const std::vector<double> v);
	double computeVariance2(const std::vector<double> v);
	void convertToRanks();
	std::vector<double> getRanksV1();
	std::vector<double> getRanksV2();
	private:
	unsigned int observationCount_;
	const std::vector<double> v1_;
	const std::vector<double> v2_;
	std::vector<double> r1_,r2_;
	const double PI = 3.141592653589793;
};
#endif
