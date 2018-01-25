#include "SPANInputGenerator.h"
#include "GenomeSizeReader.h"
#include "GTFReader.h"
#include "ExpressionReader.h"
#include "../SPAN/Data.h"
#include "../SPAN/Span.h"
/*
argv1 path to bigwig (bw) files
argv2 gene annotation in gencode format
argv3 geneID
argv4 gene expression matrix
argv5 genome size file 
*/

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
          const unsigned int window = 100000;
	const unsigned int stepSize = 1;
	const unsigned int maxCores = 3;
	bool verbose = false;
          //Loading genome size file
	GenomeSizeReader gsr(genomeSizeFile);
	gsr.loadGenomeSizeFile();
	std::map<std::string,int> genomeSize;
	genomeSize = gsr.getGenomeSize();

          //Reading the annotation file to retrieve gene coordinates
	std::cout<<"Looking up genomic coordinates for "<<geneID<<std::endl;
	GTFReader gtf(annotationFile,genomeSize,window);
	gtf.findGenomicLocation(geneID);
	std::tuple<std::string, unsigned int, unsigned int,std::string> genomicCoordinates;
	genomicCoordinates = gtf.getGenomicLocation();
	std::cout<<"Coordinates found: "<<std::get<0>(genomicCoordinates)<<" "<<std::get<1>(genomicCoordinates)+window<<" "<<std::get<2>(genomicCoordinates)-window<<",gene length: "<<std::get<2>(genomicCoordinates)-std::get<1>(genomicCoordinates)<<std::endl;
	
          //Generating expression map
	std::cout<<"Extracting gene expression information for "<<geneID<<std::endl;
	ExpressionReader expR(expressionDiscretised);
	expR.loadExpressionData(geneID);
	std::map<std::string,unsigned int> expressionMap;
	expressionMap = expR.getExpressionMap();

          //GenerateSPANInput
	std::cout<<"Generating input matrix using a search window extended by "<<window<<"bp up- and downstream of the gene"<<std::endl;
	SPANInputGenerator SPIG(bigWigPath,expressionMap);
	SPIG.generateSPANInput(genomicCoordinates);
	std::vector<std::vector<double> > perBaseInputData;
	perBaseInputData = SPIG.getInputMatrix();
	Data input = Data();
	input.setData(perBaseInputData,true,false,'g',stepSize,false);

	//Feed input to SPAN and execute
	std::cout<<"Segmentation in progress..."<<std::endl;
	SPAN sp= SPAN();
	std::vector<std::pair<unsigned int, unsigned int> > segments = sp.runSpan(input,stepSize,maxCores,verbose);
	std::vector<std::pair<unsigned int, unsigned int> > genomeConv= sp.convertSegmentationToGenomicCoordinates(segments,genomicCoordinates);
	if (verbose){
		std::cout<<SPIG<<std::endl;
		std::cout<<gsr<<std::endl;
		std::cout<<expR<<std::endl;
		std::cout<<gtf<<std::endl;
	}
	std::cout<<"Segmentation into "<<segments.size()<<" bins completed."<<std::endl;
	//Print binning
//	for (auto& element : segments){
//		std::cout<<"["<<element.first<<","<<element.second<<")"<<std::endl;
//	}

	//Print binning
	for (auto& element : genomeConv){
		std::cout<<"["<<element.first<<","<<element.second<<")"<<std::endl;
	}
          return 0;
} 
