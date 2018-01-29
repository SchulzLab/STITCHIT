#ifndef CORCOMP_h
#define CORCOMP_h

#include <vector>

class CorComp{
	public:
	CorComp(std::vector<double> v1, std::vector<double> v2);
	double computePearsonCorrelation();
	//double computeSpearmanCorrelation();

	private:
	double computeMean(const std::vector<double> v);
	double computeVariance(const std::vector<double> v);
	//void convertToRanks();
	unsigned int observationCount_;
	const std::vector<double> v1_;
	const std::vector<double> v2_;
	//std::vector<unsigned int> r1_,r2_;
};
#endif
