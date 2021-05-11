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
	GTFReader(const std::string& gtfFileName, std::map<std::string, int>& genomeSize, const int& window);
	void findGenomicLocation(const std::string& targetGeneID);
	const std::string& getGTFfileName();
	const std::tuple<std::string,unsigned int, unsigned int,std::string>& getGenomicLocation();
	const int& getWindow();
	void setWindow(int window);
	friend std::ostream& operator<<(std::ostream& os,const GTFReader& r);

	private:
	const std::string gtfFileName_;
	std::tuple<std::string,unsigned int, unsigned int,std::string> genomicLocation_;
	std::map<std::string,int> genomeSize_;
	int window_;
};

#endif
