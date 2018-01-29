#ifndef BINSELECTION_H
#define BINSELECTION_H

#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <map>
#include <tuple>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>


//IMPORTANT: The bigwig lib uses half open intervalls, e.g. [A,B)

class BinSelection{
	public:
	BinSelection(const std::string& path);
	void computeMeanSignal(std::string& chrom, const std::vector<std::pair<unsigned int, unsigned int> >& segments);
	std::vector<std::vector<double> >& getMeanSignal();
	std::vector<std::string>& getSampleNames();
	std::vector<double> computeCorrelation(std::map<std::string,double>& expressionMap);
	std::vector<double> getExpressionVectorByNames(std::map<std::string,double>& expressionMap);
	std::vector<double> getSignalVectorBySegment(unsigned int segID);

	private:
	const std::string path_;
	const std::vector<std::pair<unsigned int, unsigned int> > segmentation;
	std::vector<std::string> sampleNames_;
	std::vector<std::vector<double> > meanSignal_;	
};

#endif
