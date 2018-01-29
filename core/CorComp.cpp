#include "CorComp.h"
#include <stdexcept>
#include <math.h>
	
CorComp::CorComp(std::vector<double> v1, std::vector<double> v2)
	:v1_(v1),v2_(v2)
{
	if (v1_.size() == v2_.size())
		observationCount_= v1_.size();
	else
		throw std::runtime_error("Size of vectors is not identical");
}


double CorComp::computePearsonCorrelation(){
	double nominator = 0;
	double mv1,mv2;
	mv1 = computeMean(v1_);
	mv2 = computeMean(v2_);
	for (unsigned int i = 0; i < observationCount_; i++){
		nominator+=((v1_[i]-mv1)*(v2_[i]-mv2));
	}
	double denominator = sqrt(computeVariance(v1_))*sqrt(computeVariance(v2_));
	if (denominator == 0)
		return 0.0;
	else
		return nominator/denominator;
}


double CorComp::computeMean(const std::vector<double> v){
	double temp = 0.0;
	for (double ref : v)
		temp+=ref;
	return temp/observationCount_;
}


double CorComp::computeVariance(const std::vector<double> v){
	double mean = computeMean(v);
	double temp = 0.0;
	for (double ref : v){
		temp+=((ref-mean)*(ref-mean));
	}
	return temp;
}
