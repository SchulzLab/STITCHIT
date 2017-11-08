#include "ExpressionReader.h"

ExpressionReader::ExpressionReader(const std::string& expFileName)
	:expFileName_(expFileName)
{
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

const std::string& ExpressionReader::getExpFileName(){
	return expFileName_;

}

std::ostream& operator<<(std::ostream& os, const ExpressionReader& r){
          for (const auto& element : r.expressionMap_){
                    os << element.first << " " << element.second << std::endl;
          }
}   
