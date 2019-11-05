#ifndef REMREADER_H
#define REMREADER_H

#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <tuple>
#include <vector>

class RemReader{
	public:
	RemReader(const std::string& filename);
	std::vector<std::tuple <std::string, unsigned int, unsigned int, std::string, float,float,float,float> >& getREMs();
	const std::string& getFilename();	

	private:	
	const std::string filename_;
	std::vector<std::tuple <std::string, unsigned int, unsigned int, std::string, float,float,float,float> > REMs_;
	void readAllREMs();
};

#endif
