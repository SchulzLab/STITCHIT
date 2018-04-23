#ifndef PEAKREADER_h
#define PEAKREADER_h

#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <tuple>

class PeakReader{
	public:
	PeakReader(const std::string& filename);
	void findOverlappingPeaks(std::tuple<std::string, unsigned int, unsigned int> genomicCoordinates);
	std::vector<std::pair<unsigned int, unsigned int> >& getOverlappingPeaks();
	const std::string& getFilename();	

	private:	
	const std::string filename_;
	std::vector<std::pair<unsigned int, unsigned int> > overlappingPeaks_;
};

#endif
