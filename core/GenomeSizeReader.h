#ifndef GENOMESIZEREADER_h
#define GENOMESIZEREADER_h

#include <map>
#include <string>
#include <sstream>
#include <fstream>

class GenomeSizeReader{
	public:
	GenomeSizeReader(const std::string& filename);
	void loadGenomeSizeFile();
	std::map<std::string,unsigned int>& getGenomeSize();
	const std::string& getFilename();	

          friend std::ostream& operator<<(std::ostream& os,const GenomeSizeReader& r);

	private:	
	const std::string filename_;
	std::map<std::string,unsigned int> genomeSize_;
};

#endif
