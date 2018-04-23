#include "PeakReader.h"
#include <iostream>

PeakReader::PeakReader(const std::string& filename)
	:filename_(filename)
{
	overlappingPeaks_.clear();
}

/*! \brief Brief description.
 *         Brief description continued.
 *  Detailed description starts here.
 * @param
 * @return
 * @throw
 */
void PeakReader::findOverlappingPeaks(std::tuple<std::string, unsigned int, unsigned int> genomicCoordinates){
          std::ifstream peakFile;
          unsigned int start;
	unsigned int end;
          std::string temp;
          std::string chromosome;
	overlappingPeaks_.clear();
          peakFile.open(filename_);
          if (!peakFile) throw std::invalid_argument("Peak file "+filename_+" does not exist.");

	std::string queryChromosome=std::get<0>(genomicCoordinates);
	unsigned int queryStart=std::get<1>(genomicCoordinates);
	unsigned int queryEnd=std::get<2>(genomicCoordinates);

          while (!peakFile.eof()){
                    std::getline(peakFile,temp);
                    std::stringstream cs(temp);
                    if (temp != ""){
                              if (cs >> chromosome >> start >> end){
				if (chromosome==queryChromosome){
					if ((start < queryStart) and (end > queryStart)){
						overlappingPeaks_.push_back(std::make_pair (start,end));
					}
					else if((start > queryStart) and (start < queryEnd)){
						overlappingPeaks_.push_back(std::make_pair (start,end));
					}
					if (start > queryEnd){
						break;
					}		
				}	
                              }else{
                                        throw std::invalid_argument("Peak file "+filename_+" is not properly formatted.");
                              }
                    }
          }
          peakFile.close();

	if (overlappingPeaks_.size()==0) throw std::invalid_argument("No overlapping peaks found");
}

std::vector<std::pair<unsigned int, unsigned int> >& PeakReader::getOverlappingPeaks(){
	return overlappingPeaks_;
}


const std::string& PeakReader::getFilename(){
	return filename_;
} 

