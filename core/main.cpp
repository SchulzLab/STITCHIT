#include "SPANInputGenerator.h"
#include "GenomeSizeReader.h"
#include "GTFReader.h"
#include "ExpressionReader.h"

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
          const unsigned int window = 5000;

          //Loading genome size file
	GenomeSizeReader gsr(genomeSizeFile);
	gsr.loadGenomeSizeFile();
	std::map<std::string,unsigned int> genomeSize;
	genomeSize = gsr.getGenomeSize();

          //Reading the annotation file to retrieve gene coordinates
	GTFReader gtf(annotationFile,genomeSize,window);
	gtf.findGenomicLocation(geneID);
	std::tuple<std::string, unsigned int, unsigned int> genomicCoordinates;
	genomicCoordinates = gtf.getGenomicLocation();

          //Generating expression map
	ExpressionReader expR(expressionDiscretised);
	expR.loadExpressionData(geneID);
	std::map<std::string,unsigned int> expressionMap;
	expressionMap = expR.getExpressionMap();

          //generateSPANInput
	SPANInputGenerator SPIG(bigWigPath,expressionMap);
	SPIG.generateSPANInput(genomicCoordinates);
	std::vector<std::vector<double> > perBaseInputData;
	perBaseInputData = SPIG.getInputMatrix();

	std::cout<<SPIG;
	std::cout<<gsr;
	std::cout<<expR;
	std::cout<<gtf<<std::endl;

          return 0;
} 
