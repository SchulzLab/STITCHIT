//#include "SPANInputGenerator.h"
//#include "GenomeSizeReader.h"
//#include "GTFReader.h"
//#include "ExpressionReader.h"
#include "BinSelection.h"
#include "RemReader.h"
#include "boost/program_options/options_description.hpp"
#include "boost/program_options/variables_map.hpp"
#include "boost/program_options/parsers.hpp"
#include <map>
/*
*/

int main(int argc, char *argv[]){
	//Check the provided arguments
	//storing command args in self explanatory variables for better readability
	std::string bigWigFile;
	std::string expressionOriginal;
	std::string corM;
	std::string outputPrefix;
	std::string regionFiles;
	bool verbose;
	std::vector<std::string> fileNames;
	std::map<std::string,std::vector<std::tuple<std::string, unsigned int, unsigned int, std::string, float, float, float, float> > > regulatoryElements;
	boost::program_options::options_description desc("Options to compute the coverage");

	desc.add_options()
		("help","Show help message")
		("regionFiles,k",boost::program_options::value<std::string>(&regionFiles),"Path to gene-specific regulatory region file.")
		("bigWigFile,b", boost::program_options::value<std::string>(&bigWigFile), "Path to big wig file")
//		("expressionOriginal,o",boost::program_options::value<std::string>(&expressionOriginal), "File containng the original expression information")
//		("pvalue,p",boost::program_options::value<float>(&pvalue)->default_value(0.05),"Signifcance level for the correlation of a Peak to be considered in the output, default is 0.05")
//		("correlationMeasure,m",boost::program_options::value<std::string>(&corM)->default_value("Both"),"Method used to compute correlation between expression and epigenetic signal. Can be Both (default), Pearson, or Spearman")
		("prefix,f",boost::program_options::value<std::string>(&outputPrefix)->default_value(""),"Path were resulting files should be stored, defaults to the source directory")
		("verbose,v", boost::program_options::value<bool>(&verbose)->default_value(false), "True if additional status reports should be generated, false otherwise (default)")
	;
	
	boost::program_options::variables_map vm;
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
	boost::program_options::notify(vm);

	if (argc <= 8){
		std::cout<<desc<<std::endl;
		return 1;
	}

	if (vm.count("help")){
		std::cout<<desc<<std::endl;
		return 1;
	}

	if (not(vm.count("regionFiles"))){
		std::cout<<"--regionFiles is not specified. Please provide the path to a folder containing one bed file per sample"<<std::endl;;		
		return 1;
	}

	if (vm.count("prefix")){
		if ((outputPrefix[outputPrefix.size()]!='/') and (outputPrefix.size()>1)){
			std::cout<<"Results will be stored in "<<outputPrefix<<std::endl;
			outputPrefix+="/";
		}
	}

	if (not(vm.count("bigWigFile"))){
		std::cout<<"--bigWigFile is not specified. Please provide a file containing the read counts"<<std::endl;
		return 1;
	}

	if (boost::filesystem::is_directory(regionFiles)){
		for(auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(regionFiles), {})){
			std::string filenameS = entry.path().string();
			char * filename= new char[filenameS.size() +1];
			std::copy(filenameS.begin(),filenameS.end(), filename);
			filename[filenameS.size()]='\0';
			if (entry.path().extension().string() != ".bed") throw std::invalid_argument(filenameS+" is not a biwWig file. Expected file type is .ed");
			if (verbose) std::cout<<"Reading from "<<entry.path().stem()<<std::endl;
			RemReader pr(filename);
			regulatoryElements[filename]=pr.getREMs();		
			fileNames.push_back(filename);		
			std::cout<<"Identified "<<regulatoryElements[filename].size()<<" regulatory elements."<<std::endl;
			}
	}else{
		if (boost::filesystem::exists(regionFiles)){
			RemReader pr(regionFiles);
			fileNames.push_back(regionFiles);
			regulatoryElements[regionFiles]=pr.getREMs();
			std::cout<<"Identified "<<regulatoryElements[regionFiles].size()<<" regulatory elements within "<<regionFiles<<std::endl;
		}
	}

	if (verbose){
		std::cout<<"Chromosome	Start	End	GeneID	Reg.Coef	pValue"<<std::endl;
		for (unsigned int i=0; i<fileNames.size();i++){
			unsigned int nRem=regulatoryElements[fileNames[i]].size();
			for (unsigned int j=0; j<nRem; j++){
				std::cout<<std::get<0>(regulatoryElements[fileNames[i]][j])<<"	"<<std::get<1>(regulatoryElements[fileNames[i]][j])<<"	"<<std::get<2>(regulatoryElements[fileNames[i]][j])<<"	"<<std::get<3>(regulatoryElements[fileNames[i]][j])<<"	"<<std::get<4>(regulatoryElements[fileNames[i]][j])<<"	"<<std::get<5>(regulatoryElements[fileNames[i]][j])<<"	"<<std::endl;
			}
		}
	}

	std::map<std::string, double> expressionMapO;// = expO.getExpressionMap();


	//Compute signal in selected peak regions
	std::cout<<"Computing mean signal per unified peak and sample from "<< bigWigFile<<std::endl;
	for (unsigned int i=0; i<fileNames.size();i++){
		BinSelection bs = BinSelection(bigWigFile,expressionMapO);
		bs.computeMeanSignal(regulatoryElements[fileNames[i]]);
	}

	std::cout<<"Chromosome	Start	End	GeneID	Reg.Coef	pValue	Signal	Combined"<<std::endl;
	for (unsigned int i=0; i<fileNames.size();i++){
		unsigned int nRem=regulatoryElements[fileNames[i]].size();
		for (unsigned int j=0; j<nRem; j++){
			std::cout<<std::get<0>(regulatoryElements[fileNames[i]][j])<<"	"<<std::get<1>(regulatoryElements[fileNames[i]][j])<<"	"<<std::get<2>(regulatoryElements[fileNames[i]][j])<<"	"<<std::get<3>(regulatoryElements[fileNames[i]][j])<<"	"<<std::get<4>(regulatoryElements[fileNames[i]][j])<<"	"<<std::get<5>(regulatoryElements[fileNames[i]][j])<<"	"<<std::get<6>(regulatoryElements[fileNames[i]][j])<<"	"<<std::get<7>(regulatoryElements[fileNames[i]][j])<<std::endl;
		}
	}
	return 0;
} 
