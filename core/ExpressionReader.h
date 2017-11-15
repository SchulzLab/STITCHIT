#ifndef EXPRESSIONREADER_H
#define EXPRESSIONREADER_H

#include <string>
#include <map>
#include <algorithm>
#include <stdlib.h>
#include <fstream>
#include <sstream>

class ExpressionReader{
	public:
	ExpressionReader(const std::string& expFileName);
	void loadExpressionData(const std::string& targetGeneID);
	std::map<std::string,unsigned int>& getExpressionMap();
	const std::string& getFilename();

          friend std::ostream& operator<<(std::ostream& os,const ExpressionReader& r);


	private:
	const std::string expFileName_;
	std::map<std::string,unsigned int> expressionMap_; 
};

#endif
