#include "CorComp.h"
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <map>
#include <limits>	
#include <iostream>

CorComp::CorComp(std::vector<double> v1, std::vector<double> v2)
	:v1_(v1),v2_(v2)
{
	if ((v1_.size() == v2_.size()) and (v1.size() > 0))
		observationCount_= v1_.size();
	else
		if (v1.size() > 0)
			throw std::runtime_error("Size of vectors is not identical");
		else
			throw std::runtime_error("Can not compute correlation on empty vectors");
}

double CorComp::computePearsonCorrelation(){
	double nominator = 0;
	double mv1,mv2;
	mv1 = computeMean(v1_);
	mv2 = computeMean(v2_);
	for (unsigned int i = 0; i < observationCount_; i++){
		nominator+=((v1_[i]-mv1)*(v2_[i]-mv2));
	}
	double denominator = sqrt(computeVariance2(v1_))*sqrt(computeVariance2(v2_));
	if (denominator == 0)
		return 0.0;
	else
		return nominator/denominator;
}

double CorComp::computeSpearmanCorrelation(){
	convertToRanks();
	double nominator = 0;
	double mv1,mv2;
	mv1 = computeMean(r1_);
	mv2 = computeMean(r2_);
	for (unsigned int i = 0; i < observationCount_; i++){
		nominator+=((r1_[i]-mv1)*(r2_[i]-mv2));
	}
	double denominator = sqrt(computeVariance2(r1_))*sqrt(computeVariance2(r2_));
	if (denominator == 0)
		return 0.0;
	else
		return nominator/denominator;
}

const unsigned int CorComp::getObservationCount(){
	return observationCount_;
}

double CorComp::computeMean(const std::vector<double> v){
	double temp = 0.0;
	for (double ref : v)
		temp+=ref;
	return temp/observationCount_;
}

double CorComp::computeVariance2(const std::vector<double> v){
	double mean = computeMean(v);
	double temp = 0.0;
	for (double ref : v){
		temp+=((ref-mean)*(ref-mean));
	}
	return temp;
}

double CorComp::computeVariance(const std::vector<double> v){
	double mean = computeMean(v);
	double temp = 0.0;
	for (double ref : v){
		temp+=((ref-mean)*(ref-mean));
	}
	return temp/(observationCount_-1);
}

void CorComp::convertToRanks(){
	std::multimap<double,double> mv1;
	std::multimap<double,double> mv2;
	for (double element : v1_){
		mv1.insert(std::pair<double,double>(element,0));
	}
	for (double element : v2_){
		mv2.insert(std::pair<double,double>(element,0));
	}

	
	std::multimap<double,double>::iterator rit;
	double rankV=1;
	for (rit=mv1.begin(); rit!=mv1.end(); ++rit){
		rit->second=rankV;
		rankV=rankV+1;
	}

	rankV=1;
	for (rit=mv2.begin(); rit!=mv2.end(); ++rit){
		rit->second=rankV;
		rankV=rankV+1;
	}

	double counter,sum;
	for(double element :v1_){
		std::pair <std::multimap<double,double>::iterator, std::multimap<double,double>::iterator> ret;
		ret = mv1.equal_range(element);
		counter=0;
		sum=0;
		for (std::multimap<double,double>::iterator it=ret.first; it!=ret.second; ++it){
			counter=counter+1;
			sum=sum+it->second;
		}
		r1_.push_back(sum/counter);
	 }

	for(double element :v2_){
		std::pair <std::multimap<double,double>::iterator, std::multimap<double,double>::iterator> ret;
		ret = mv2.equal_range(element);
		counter=0;
		sum=0;
		for (std::multimap<double,double>::iterator it=ret.first; it!=ret.second; ++it){
			counter=counter+1;
			sum=sum+it->second;
		}
		r2_.push_back(sum/counter);
	 }
}

std::vector<double> CorComp::getRanksV1(){
	return r1_;
}

std::vector<double> CorComp::getRanksV2(){
	return r2_;
}

std::pair<double, double> CorComp::getFisherZ(double cor){
	double error;
	double cor2 = std::abs(cor);
	if (observationCount_ < 5)
		error=1.0;
	else
		error=1.0/sqrt(observationCount_-3.0);
	if (std::abs(cor)==1.0)
		return std::make_pair(std::numeric_limits<double>::max(),error);
	else
		return std::make_pair((0.5)*log((1.0+cor2)/(1.0-cor2)),error);
}

double CorComp::cdf(double x){
	if (x > 4.0)
		return 1.0;
	else{
	double sum=x;
	long double value=x;
	unsigned long long int df = 1;
	for (unsigned int i = 1; i< 16; i++){
		value=value*x*x;
		df=df*(2*i+1);
		sum=sum+(value/df);
	}
	return 0.5+(1.0/sqrt(2.0*PI))*exp((-1*x*x)/2.0)*sum;
	}
}

double CorComp::getPvalue(double cor){
	std::pair<double,double> fisher = getFisherZ(cor);
	return 1.0-cdf(fisher.first/fisher.second);
}
