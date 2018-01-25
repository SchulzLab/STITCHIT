#ifndef SPANINPUTGENERATOR_H
#define SPANINPUTGENERATOR_H

#include "../BigWig/bigWig.h"
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


class SPANInputGenerator{
	public:
	SPANInputGenerator(const std::string& path, std::map<std::string, unsigned int> expressionMap);
	void generateSPANInput(const std::tuple<std::string, unsigned int, unsigned int,std::string>& genomicCoordinates);
	std::vector<std::vector<double>>& getInputMatrix();	
	std::vector<std::string>& getSampleNames();
          friend std::ostream& operator<<(std::ostream& os,const SPANInputGenerator& r);

	private:
	std::vector<double> parseIntervals(bwOverlappingIntervals_t *ints, uint32_t start, unsigned int exp, const unsigned int size);
	void printIntervals(bwOverlappingIntervals_t *ints, uint32_t start, unsigned int exp);

	const std::string path_;
	std::vector<std::string> sampleNames_;
	std::map<std::string, unsigned int> expressionMap_;
	std::vector<std::vector<double>> inputMatrix_;

};

#endif
