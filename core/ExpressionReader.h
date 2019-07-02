#ifndef EXPRESSIONREADER_H
#define EXPRESSIONREADER_H

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <stdlib.h>
#include <fstream>
#include <sstream>

class ExpressionReader{
	public:
	ExpressionReader(const std::string& expFileName);
	void loadExpressionData(const std::string& targetGeneID,bool log2Transform);
	std::map<std::string, double>& getExpressionMap();
	const std::string& getFilename();
	void checkDiversity();
          friend std::ostream& operator<<(std::ostream& os,const ExpressionReader& r);


	private:
	const std::string expFileName_;
	std::map<std::string, double> expressionMap_; 
};

#endif
