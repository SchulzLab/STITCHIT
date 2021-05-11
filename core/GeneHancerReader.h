#ifndef GENEHANCERREADER_H
#define GENEHANCERREADER_H

#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <tuple>
#include <vector>

class GeneHancerReader{
	public:
	GeneHancerReader(const std::string& filename);
	void findEnhancers(const std::string& geneID);
	std::vector<std::tuple<std::string, unsigned int, unsigned int> >& getEnhancers();
	const std::string& getFilename();	

	private:	
	const std::string filename_;
	std::vector<std::tuple<std::string, unsigned int, unsigned int> > overlappingEnhancers_;
};

#endif
