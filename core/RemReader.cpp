#include "RemReader.h"
#include <iostream>

RemReader::RemReader()
{

}


RemReader::RemReader(const std::string& filename)
	:filename_(filename)
{
	readAllREMs();
}

void RemReader::readBEDFile(const std::string& filename){
	std::ifstream remFile;
	REMs_.clear();
	std::string chromosome;
	unsigned int start;
	unsigned int end;
	std::string temp;
	remFile.open(filename);
	if (!remFile) throw std::invalid_argument("REM file "+filename+" does not exist.");
	while (!remFile.eof()){
		std::getline(remFile,temp);
		std::stringstream cs(temp);
		if (temp != ""){
			if (cs >> chromosome >> start >> end){
				BED_.push_back(std::make_tuple(chromosome, start, end));
			}else{
				throw std::invalid_argument("REM file "+filename+" is not properly formatted.");
			}
		}
	}
	remFile.close();
	if (BED_.size()==0) throw std::invalid_argument("No entries found in bed file");
}

	
std::vector<std::tuple <std::string, unsigned int, unsigned int> >& RemReader::getSelectedREMs(){
	return BED_;
}



/*! \brief Brief description.
 *         Brief description continued.
 *  Detailed description starts here.
 * @param
 * @return
 * @throw
 */
void RemReader::readAllREMs(){
	std::ifstream remFile;
	REMs_.clear();
	std::string chromosome;
	unsigned int start;
	unsigned int end;
	std::string geneID;
	float score;
	float pValue;
	std::string temp;
	remFile.open(filename_);
	if (!remFile) throw std::invalid_argument("REM file "+filename_+" does not exist.");
	while (!remFile.eof()){
		std::getline(remFile,temp);
		std::stringstream cs(temp);
		if (temp != ""){
			if (cs >> chromosome >> start >> end >> geneID >> score >> pValue){
				REMs_.push_back(std::make_tuple(chromosome, start, end, geneID, score, pValue,0.0,0.0));
			}else{
				throw std::invalid_argument("REM file "+filename_+" is not properly formatted.");
			}
		}
	}
	remFile.close();
	if (REMs_.size()==0) throw std::invalid_argument("No REMs found");
}

std::vector<std::tuple<std::string, unsigned int, unsigned int, std::string, float,float,float,float> >& RemReader::getREMs(){
	return REMs_;
}


const std::string& RemReader::getFilename(){
	return filename_;
} 

