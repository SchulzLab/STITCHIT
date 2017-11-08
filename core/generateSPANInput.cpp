#include "../BigWig/bigWig.h"
#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <map>
#include <tuple>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
/*
argv1 path to bigwig (bw) files
argv2 gene annotation in gencode format
argv3 geneID
argv4 gene expression matrix
argv5 genome size file 
*/

/*! \brief Brief description.
 *
 *  Detailed description starts here.
 * @param
 * @param
 */
void printIntervals(bwOverlappingIntervals_t *ints, uint32_t start, unsigned int exp) {
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
	printf("%f\t",exp);
	printf("\n");
}

/*! \brief Brief description.
 *
 *  Detailed description starts here.
 * @param
 * @param
 * @return
 */
std::vector<double> parseIntervals(bwOverlappingIntervals_t *ints, uint32_t start, unsigned int exp, const unsigned int size) {
	std::vector<double> temp;
	temp.reserve(size+1);
	uint32_t i;
	if(!ints) throw std::invalid_argument("Error in the provided genomic intervals");
	for(i=0; i<ints->l; i++) {
		if(ints->start && ints->end) {
			temp.push_back(ints->value[i]);
		} else if(ints->start) {
			temp.push_back(ints->value[i]);
		} else {
			temp.push_back(ints->value[i]);
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
const std::map<std::string,unsigned int> loadGenomeSizeFile(const std::string& filename){
	std::map<std::string,unsigned int> tempGenomeSize;
	std::ifstream genomeSizeFile;
	unsigned int size;
	std::string temp;
	std::string chromosome;
	genomeSizeFile.open(filename);
	if (!genomeSizeFile) throw std::invalid_argument("Genome size file "+filename+" is not properly formatted.");
	
	while (!genomeSizeFile.eof()){
		std::getline(genomeSizeFile,temp);
		std::stringstream cs(temp);
		if (temp != ""){
			if (cs >> chromosome >> size){
				tempGenomeSize[chromosome]=size;
			}else{
				throw std::invalid_argument("Genome size file "+filename+" is not properly formatted.");
			}
		}
	}
	genomeSizeFile.close();
	return tempGenomeSize;
}

/*! \brief Brief description.
 *         Brief description continued.
 *  Detailed description starts here.
 * @param
 * @return
 * @throw
 */
std::tuple<std::string,unsigned int, unsigned int> getGenomicLocation(const std::string& targetGeneID, const std::string& gtfFileName, std::map<std::string,unsigned int>& genomeSize, const unsigned int& window ){
	unsigned int pos1, pos2;
	std::string chromosome, buf, temp, convertedGeneID, geneID;
	std::ifstream annotationFile;
	const unsigned int zero = 0;
	annotationFile.open(gtfFileName);
	if (!annotationFile) throw std::invalid_argument("Gene annotation file "+gtfFileName+" could not be loaded");

	while(!annotationFile.eof()){
		std::getline(annotationFile,temp);
		std::stringstream ss(temp); // Insert the string into a stream
		if (temp != ""){
			if (ss >> chromosome >> buf >> buf){
				if (buf == "gene"){
					if (ss >> pos1 >> pos2 >> buf >> buf >> buf >> buf >> geneID){
						convertedGeneID=geneID.substr(1,15);	
						if (convertedGeneID==targetGeneID){
							pos1=std::max(zero,pos1-window);
							pos2=std::min(genomeSize[chromosome],pos2+window);
							return std::make_tuple(chromosome, pos1, pos2);
						}
					}
					else{
						throw std::invalid_argument("Annotation file "+gtfFileName+" does not match specifications");
					}
				}
			}
			else {
				throw std::invalid_argument("Annotation file "+gtfFileName+" does not match specifications");
			}
		}
	}		
	annotationFile.close();
	//Exit if the gene could not be found
	throw std::invalid_argument("The provided geneID "+targetGeneID+" could not be found in the annotation file");
}



/*! \brief Brief description.
 *         Brief description continued.
 *  Detailed description starts here.
 * @param
 * @return
 * @throw
 */
std::map<std::string,unsigned int> getExpressionMap(const std::string& targetGeneID, const std::string& expFileName){
	std::ifstream expressionFile;
	expressionFile.open(expFileName);
	if (!expressionFile) throw std::invalid_argument("Expression file "+expFileName+" could not be opened");


	//Reading sample Names
	std::string temp, buf;
	std::getline(expressionFile,temp);
	std::vector<std::string> sampleNames;
	std::stringstream sN(temp);
	while(sN >> buf){
		sampleNames.push_back(buf);
	}	
	//Generating expression map
	std::map<std::string,unsigned int> expressionMapTemp;
	unsigned int value;
	unsigned int counter=0;
	while (!expressionFile.eof()){
		std::getline(expressionFile,temp);
		std::stringstream sE(temp);
		if (temp != ""){
			if (sE >> buf){
				if (buf == targetGeneID){
					while (sE >> value){
						expressionMapTemp[sampleNames[counter]]=value;
						counter+=1;
					}					
					expressionFile.close();
					return expressionMapTemp;
				}
			}
			else{
				throw std::invalid_argument("Expression file "+expFileName+" is not properly formatted");
			}
		}
	}
	expressionFile.close();
	throw std::invalid_argument("Expression data for the gene "+targetGeneID+" could not be found");
	}

/*! \brief Brief description.
 *         Brief description continued.
 *  Detailed description starts here.
 * @param
 * @return
 * @throw
 */
std::vector<std::vector<double>> generateSPANInput(const std::string& path, const std::tuple<std::string, unsigned int, unsigned int>& genomicCoordinates, std::map<std::string, unsigned int> expressionMap){
	//Generating per base input for SPAN
	if(bwInit(1<<17) != 0) throw std::invalid_argument("Error occured in bwInit.");
	std::vector<std::vector<double> > matrix;
	//Iterating through all files in the specified directory
	if (boost::filesystem::is_directory(path)){
		for(auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(path), {})){
			bigWigFile_t *fp = NULL;
			bwOverlappingIntervals_t *intervals = NULL;
			std::string filenameS = entry.path().string();
			char * filename= new char[filenameS.size() +1];
			std::copy(filenameS.begin(),filenameS.end(), filename);
			filename[filenameS.size()]='\0';

			//Check for presence of expression data for the current sample, otherwise ignore the sample!
			if (expressionMap.find(entry.path().stem().string()) != expressionMap.end()){
				std::vector<double> currentVec;
				//Openning bw file and printing data
				fp = bwOpen(filename, NULL, "r");
				if(!fp) throw std::invalid_argument("The bigwig file "+filenameS+" could not be opened.");
				
				char * chromosome_c_str = new char[std::get<0>(genomicCoordinates).size() + 1];
				std::copy(std::get<0>(genomicCoordinates).begin(), std::get<0>(genomicCoordinates).end(), chromosome_c_str);
				chromosome_c_str[std::get<0>(genomicCoordinates).size()] = '\0'; // don't forget the terminating 0
				intervals = bwGetValues(fp,chromosome_c_str,std::get<1>(genomicCoordinates),std::get<2>(genomicCoordinates),1);
				currentVec=parseIntervals(intervals,0,expressionMap[entry.path().stem().string()],std::get<2>(genomicCoordinates)-std::get<1>(genomicCoordinates));
				printIntervals(intervals,0,expressionMap[entry.path().stem().string()]);
				bwDestroyOverlappingIntervals(intervals);
				delete[] chromosome_c_str;
				bwClose(fp);
				bwCleanup();
				matrix.push_back(currentVec);
			}else{
				std::cerr<<"No expression information available for "<<filenameS<<std::endl;
			}
		}
	}else{
		throw std::invalid_argument("You did not specifiy a directory containing bigwig files.");
	}
	return matrix;
}

int main(int argc, char *argv[]){
	//Check the provided arguments
	if(argc != 6) {
		std::cerr<<"Incorrect number of parameters. Usage:"<<std::endl;
		std::cerr<<"<Path to bigwig files> <Gene annotation file in gencode format> <geneID> <Discretised gene expression matrix> <Genome size file>"<<std::endl;
		return 1;
	}

	//storing command args in self explanatory variables for better readability
	const std::string bigWigPath = argv[1];
	const std::string annotationFile = argv[2];
	const std::string geneID = argv[3];
	const std::string expressionDiscretised = argv[4];
	const std::string genomeSizeFile = argv[5];
	const unsigned int window = 5000;

	//Loading genome size file
	std::map<std::string,unsigned int> genomeSize;
	genomeSize = loadGenomeSizeFile(genomeSizeFile);

	//Reading the annotation file to retrieve gene coordinates
	std::tuple<std::string, unsigned int, unsigned int> genomicCoordinates;
	genomicCoordinates = getGenomicLocation(geneID,annotationFile,genomeSize,window);

	//Generating expression map
	std::map<std::string,unsigned int> expressionMap;
	expressionMap = getExpressionMap(geneID,expressionDiscretised);

	//generateSPANInput
	std::vector<std::vector<double> > perBaseInputData;
	perBaseInputData = generateSPANInput(bigWigPath,genomicCoordinates,expressionMap);

	return 0;
}
