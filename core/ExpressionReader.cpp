#include "ExpressionReader.h"
#include <iostream>

ExpressionReader::ExpressionReader(const std::string& expFileName)
	:expFileName_(expFileName)
{
	expressionMap_.clear();
}
          
/*! \brief Brief description.
 *         Brief description continued.
 *  Detailed description starts here.
 * @param
 * @return
 * @throw
 */
void ExpressionReader::loadExpressionData(const std::string& targetGeneID){
          std::ifstream expressionFile;
          expressionFile.open(expFileName_);
          if (!expressionFile) throw std::invalid_argument("Expression file "+expFileName_+" could not be opened");

	expressionMap_.clear();
          //Reading sample Names
          std::string temp, buf;
          std::getline(expressionFile,temp);
          std::vector<std::string> sampleNames;
          std::stringstream sN(temp);
          while(sN >> buf){
                    sampleNames.push_back(buf);
          }         
	unsigned int numberSamples = sampleNames.size();
          //Generating expression map
          unsigned int value;
          unsigned int counter=0;
	bool flag=false;
          while (!expressionFile.eof()){
                    std::getline(expressionFile,temp);
                    std::stringstream sE(temp);
                    if (temp != ""){
                              if (sE >> buf){
			          if (buf == targetGeneID){
                                                  while (sE >> value){
                                                            expressionMap_[sampleNames[counter]]=value;
                                                            counter+=1;
                                                  }                                                 
					flag=true;
                                                  expressionFile.close();
					if (counter != numberSamples){
						throw std::invalid_argument("The number of gene expression values for gene "+targetGeneID+" does not match the number of expected entries.");
                                        	}
                              	}
			}
                              else{
                                        throw std::invalid_argument("Expression file "+expFileName_+" is not properly formatted");
                              	}
                    	}
          	}
	if (!flag){
	          expressionFile.close();
	          throw std::invalid_argument("Expression data for the gene "+targetGeneID+" could not be found");
	}
}

std::map<std::string,unsigned int>& ExpressionReader::getExpressionMap(){
	return expressionMap_;
}

const std::string& ExpressionReader::getFilename(){
	return expFileName_;

}

std::ostream& operator<<(std::ostream& os, const ExpressionReader& r){
          for (const auto& element : r.expressionMap_){
                    os << element.first << " " << element.second << std::endl;
          }
	return os;
}   
