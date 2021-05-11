#include "GenomeSizeReader.h"
#include <iostream>

GenomeSizeReader::GenomeSizeReader(const std::string& filename)
	:filename_(filename)
{
	genomeSize_.clear();
}

/*! \brief Brief description.
 *         Brief description continued.
 *  Detailed description starts here.
 * @param
 * @return
 * @throw
 */
void GenomeSizeReader::loadGenomeSizeFile(){
          std::ifstream genomeSizeFile;
          unsigned int size;
          std::string temp;
          std::string chromosome;
	genomeSize_.clear();
          genomeSizeFile.open(filename_);
          if (!genomeSizeFile) throw std::invalid_argument("Genome size file "+filename_+" does not exist.");

          while (!genomeSizeFile.eof()){
                    std::getline(genomeSizeFile,temp);
                    std::stringstream cs(temp);
                    if (temp != ""){
                              if (cs >> chromosome >> size){
                                        genomeSize_[chromosome]=size;
                              }else{
                                        throw std::invalid_argument("Genome size file "+filename_+" is not properly formatted.");
                              }
                    }
          }
          genomeSizeFile.close();
}

std::map<std::string, int>& GenomeSizeReader::getGenomeSize(){
	return genomeSize_;
}


const std::string& GenomeSizeReader::getFilename(){
	return filename_;
} 

std::ostream& operator<<(std::ostream& os, const GenomeSizeReader& r){
	for (const auto& element : r.genomeSize_){
          	os << element.first << " " << element.second << std::endl;
	}
	return os;
} 
