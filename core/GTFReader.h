#ifndef GTFREADER_H
#define GTFREADER_H

#include <string>
#include <map>
#include <algorithm>
#include <tuple>
#include <stdlib.h> 
#include <fstream>
#include <sstream>

class GTFReader{
	public:
	GTFReader(const std::string& gtfFileName, std::map<std::string, unsigned int>& genomeSize, const unsigned int& window);
	void findGenomicLocation(const std::string& targetGeneID);
	const std::string& getGTFfileName();
	const std::tuple<std::string,unsigned int, unsigned int>& getGenomicLocation();
	const unsigned int& getWindow();

	friend std::ostream& operator<<(std::ostream& os,const GTFReader& r);

	private:
	const std::string gtfFileName_;
	std::tuple<std::string,unsigned int, unsigned int> genomicLocation_;
	std::map<std::string,unsigned int> genomeSize_;
	const unsigned int window_;
};

#endif
