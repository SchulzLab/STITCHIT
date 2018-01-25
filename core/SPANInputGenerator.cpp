#include "SPANInputGenerator.h"

SPANInputGenerator::SPANInputGenerator(const std::string& path, std::map<std::string, unsigned int> expressionMap)
	:path_(path),
	expressionMap_(expressionMap)
{
	sampleNames_.clear();
	inputMatrix_.clear();
}

/*! \brief Brief description.
 *
 *  Detailed description starts here.
 * @param
 * @param
 */
void SPANInputGenerator::printIntervals(bwOverlappingIntervals_t *ints, uint32_t start, unsigned int exp) {
          uint32_t i;
          if(!ints) return;
          for(i=0; i<ints->l; i++) {
                    if(ints->start && ints->end) {
                              printf("%f\t",ints->value[i]);
                    } else if(ints->start) {
                              printf("%f\t",ints->value[i]);
                    } else {
                              printf("%f\t",ints->value[i]);
                    }
          }
          printf("%u\t",exp);
          printf("\n");
}

/*! \brief Brief description.
 *
 *  Detailed description starts here.
 * @param
 * @param
 * @return
 */
std::vector<double> SPANInputGenerator::parseIntervals(bwOverlappingIntervals_t *ints, uint32_t start, unsigned int exp, const unsigned int size) {
          std::vector<double> temp;
	double tempV;
          temp.reserve(size+1);
          uint32_t i;
          if(!ints) throw std::invalid_argument("Error in retrieving counts from bw files, please ensure the target chromosome is  present in all bw files");
          for(i=0; i<ints->l; i++) {
		if (isnan(ints->value[i])){
			tempV=0;
		} else tempV=ints->value[i];

                    if(ints->start && ints->end) {
                              temp.push_back(tempV);
                    } else if(ints->start) {
                              temp.push_back(tempV);
                    } else {
                              temp.push_back(tempV);
                    }
          }
          return temp;
}



/*! \brief Brief description.
 *         Brief description continued.
 *  Detailed description starts here.
 * @param
 * @return
 * @throw
 */
void SPANInputGenerator::generateSPANInput(const std::tuple<std::string, unsigned int, unsigned int,std::string>& genomicCoordinates){
          //Generating per base input for SPAN
          if(bwInit(1<<17) != 0) throw std::invalid_argument("Error occured in bwInit.");
          //Iterating through all files in the specified directory
          if (boost::filesystem::is_directory(path_)){
                    for(auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(path_), {})){
                              bigWigFile_t *fp = NULL;
                              bwOverlappingIntervals_t *intervals = NULL;
                              std::string filenameS = entry.path().string();
                              char * filename= new char[filenameS.size() +1];
                              std::copy(filenameS.begin(),filenameS.end(), filename);
                              filename[filenameS.size()]='\0';
			if (entry.path().extension().string() != ".bw"){
				throw std::invalid_argument(filenameS+" is not a biwWig file. Expected file type is .bw");
			}
                              //Check for presence of expression data for the current sample, otherwise ignore the sample!
                              if (expressionMap_.find(entry.path().stem().string()) != expressionMap_.end()){
				sampleNames_.push_back(entry.path().stem().string());
                                        std::vector<double> currentVec;
                                        //Openning bw file and printing data
                                        fp = bwOpen(filename, NULL, "r");
                                        if(!fp) throw std::invalid_argument("The bigwig file "+filenameS+" could not be opened.");
                                        
                                        char * chromosome_c_str = new char[std::get<0>(genomicCoordinates).size() + 1];
                                        std::copy(std::get<0>(genomicCoordinates).begin(), std::get<0>(genomicCoordinates).end(), chromosome_c_str);
                                        chromosome_c_str[std::get<0>(genomicCoordinates).size()] = '\0'; // don't forget the terminating 0
                                        intervals = bwGetValues(fp,chromosome_c_str,std::get<1>(genomicCoordinates),std::get<2>(genomicCoordinates),1);
                                        currentVec=parseIntervals(intervals,0,expressionMap_[entry.path().stem().string()],(std::get<2>(genomicCoordinates))-(std::get<1>(genomicCoordinates)));
                                        bwDestroyOverlappingIntervals(intervals);
                                        delete[] chromosome_c_str;
                                        bwClose(fp);
                                        bwCleanup();
				currentVec.push_back(expressionMap_[entry.path().stem().string()]);
                                        inputMatrix_.push_back(currentVec);
                              }else{
                                        throw std::invalid_argument("No expression information available for "+filenameS);
                              }
                    }
          }else{
                    throw std::invalid_argument("You did not specifiy a directory containing bigwig files.");
          }
}

std::vector<std::string>& SPANInputGenerator::getSampleNames(){
	return sampleNames_;
}

std::vector<std::vector<double>>& SPANInputGenerator::getInputMatrix(){
	return inputMatrix_;
}

std::ostream& operator<<(std::ostream& os, const SPANInputGenerator& r){
	unsigned int counter=0;
	for (const auto& firstV : r.inputMatrix_){
		os << r.sampleNames_[counter++]<<" ";
		for (const auto& element : firstV){
			os<<element<<" ";
		}
		os<<std::endl;
	}
	return os;
}
