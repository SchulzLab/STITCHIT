#include "GeneHancerReader.h"
#include <iostream>

GeneHancerReader::GeneHancerReader(const std::string& filename)
	:filename_(filename)
{
	overlappingEnhancers_.clear();
}

/*! \brief Brief description.
 *         Brief description continued.
 *  Detailed description starts here.
 * @param
 * @return
 * @throw
 */
void GeneHancerReader::findEnhancers(const std::string& geneID){
          std::ifstream peakFile;
          unsigned int start;
	unsigned int end;
          std::string chromosome;
	std::string currentGeneID;
	std::string temp;
	overlappingEnhancers_.clear();
          peakFile.open(filename_);
          if (!peakFile) throw std::invalid_argument("Peak file "+filename_+" does not exist.");

	auto found = geneID.find("ENSG");
	bool flag = false;
	if (found ==std::string::npos) throw std::invalid_argument("GeneID not properly formatted."); 
          while (!peakFile.eof()){
                    std::getline(peakFile,temp);
                    std::stringstream cs(temp);
                    if (temp != ""){
                              if (cs >> chromosome >> start >> end >> currentGeneID){
				if (currentGeneID==geneID){
					flag=true;
					overlappingEnhancers_.push_back(std::make_tuple(chromosome, start, end));	
				}
				else if (flag==true){
					break;
				}				
                              }else{
                                        throw std::invalid_argument("Enhancer file "+filename_+" is not properly formatted.");
                              }
                    }
          }
          peakFile.close();

	if (overlappingEnhancers_.size()==0) throw std::invalid_argument("No enhancer regions found");
}

std::vector<std::tuple<std::string, unsigned int, unsigned int> >& GeneHancerReader::getEnhancers(){
	return overlappingEnhancers_;
}

const std::string& GeneHancerReader::getFilename(){
	return filename_;
} 
