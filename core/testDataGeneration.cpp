#include "ExpressionReader.h"
#include "BinSelection.h"
#include "RemReader.h"
#include "boost/program_options/options_description.hpp"
#include "boost/program_options/variables_map.hpp"
#include "boost/program_options/parsers.hpp"
#include <random>
#include <algorithm>
#include <cctype>

/*
argv1 path to bigwig (bw) files
argv2 gene annotation in gencode format
argv3 geneID
argv4 gene expression matrix
argv5 genome size file 
*/

int main(int argc, char *argv[]){
	//Check the provided arguments
	//storing command args in self explanatory variables for better readability
	std::string regionFile;
	std::string bigWigPath;
	std::string geneID;
	std::string expressionOriginal;
	std::string outputPrefix;
	std::string sampleVector;
	bool verbose;
	std::vector<std::tuple<std::string, unsigned int, unsigned int> > regulatoryElements;

	boost::program_options::options_description desc("STITCHIT options");
	desc.add_options()
		("help","Show help message")
		("regionFile,k",boost::program_options::value<std::string>(&regionFile),"Path to a gene-specific regulatory region file.")
		("bigWigPath,b", boost::program_options::value<std::string>(&bigWigPath), "Path to big wig files")
		("geneID,g",boost::program_options::value<std::string>(&geneID), "ID of the gene that should be segmented")
		("expressionOriginal,o",boost::program_options::value<std::string>(&expressionOriginal), "File containng the original expression information")
		("prefix,f",boost::program_options::value<std::string>(&outputPrefix)->default_value(""),"Path were resulting files should be stored, defaults to STITCH source directory")
		("verbose,v", boost::program_options::value<bool>(&verbose)->default_value(false), "True if additional status reports should be generated, false otherwise (default)")
		("sampleVector,y",boost::program_options::value<std::string>(&sampleVector),"Binary vector indicating which samples are to be used for training. All samples are used if this option is not specified (Default NULL).")
	;
	
	boost::program_options::variables_map vm;
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
	boost::program_options::notify(vm);

	if (argc == 1){
		std::cout<<desc<<std::endl;
		return 1;
	}

	if (vm.count("help")){
		std::cout<<desc<<std::endl;
		return 1;
	}
	
	if (not(vm.count("regionFile"))){
		std::cout<<"--regionFile is not specified. Please provide the path to a folder containing one bed file per sample"<<std::endl;;		
		return 1;
	}

	if (not(vm.count("bigWigPath"))){
		std::cout<<"--bigWigPath is not specified. Please provide the path to a folder containing one big wig file per sample"<<std::endl;;		
		return 1;
	}

	if (not(vm.count("geneID"))){
		std::cout<<"--geneID is not specified. Please provide the geneID of the gene that should be segmented."<<std::endl;;		
		return 1;
	}

	if (not(vm.count("expressionOriginal"))){
		std::cout<<"--expressionOriginal is not specified. Please provide the path to a folder containing one big wig file per sample"<<std::endl;;		
		return 1;
	}

	if (vm.count("prefix")){
		if ((outputPrefix[outputPrefix.size()]!='/') and (outputPrefix.size()>1)){
			if (verbose){
				std::cout<<"Results will be stored in "<<outputPrefix<<std::endl;
			}
			outputPrefix+="/";
		}
	}

	if (not(vm.count("sampleVector"))){
		std::cout<<"sampleVector is not specified. Please provide a binary vector indicating which samples in the bigWigPath are to be used as test cases for the linear model"<<std::endl;;		
		return 1;

	}

	//Loading original gene expression data
	if(verbose){
		std::cout<<"Extracting original gene expression information for "<<geneID<<std::endl;
	}
	ExpressionReader expO(expressionOriginal);
	expO.loadExpressionData(geneID,false);
	std::map<std::string, double> expressionMapO = expO.getExpressionMap();

	//Determining sample numbers of test samples
	std::vector<unsigned int> sVector;
	sVector.resize(0);
	for (unsigned int i=0; i<sampleVector.size(); i++){
		if (sampleVector[i]=='1'){
			sVector.push_back(i);
		}
	}

	//Read regulatory elements from REM file
	if (boost::filesystem::is_directory(regionFile)) throw std::invalid_argument(regionFile+" is not a bed file but a directory");
	else{
		if (boost::filesystem::exists(regionFile)){
			RemReader pr;
			pr.readBEDFile(regionFile);
			// New function RemReader getSelectedREMs
			regulatoryElements=pr.getSelectedREMs();
			if(verbose)std::cout<<"Identified "<<regulatoryElements.size()<<" regulatory elements within "<<regionFile<<std::endl;
		}
	}


	//Compute signal in segmented regions
	if (verbose){
		std::cout<<"Computing mean signal per bin and sample"<<std::endl;
	}
	BinSelection bs = BinSelection(bigWigPath,expressionMapO);
	bs.computeMeanSignal(regulatoryElements, sVector);
	if (verbose){
		std::cout<<"Computing correlation between signal and gene expression"<<std::endl;
	}
	//New function BinSelection storeSignal
	bs.storeSignal(outputPrefix+"Segmentation_"+geneID+"_TestData.txt", regulatoryElements);

	return 0;
} 


