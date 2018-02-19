#include "SPANInputGenerator.h"
#include "GenomeSizeReader.h"
#include "GTFReader.h"
#include "ExpressionReader.h"
#include "BinSelection.h"
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
          if(argc != 7) {
                    std::cerr<<"Incorrect number of parameters. Usage:"<<std::endl;
                    std::cerr<<"<Path to bigwig files> <Gene annotation file in gencode format> <geneID> <Discretised gene expression matrix> <Original gene expression matrix> <Genome size file>"<<std::endl;
                    return 1;
          }

          //storing command args in self explanatory variables for better readability
          const std::string bigWigPath = argv[1];
          const std::string annotationFile = argv[2];
          const std::string geneID = argv[3];
          const std::string expressionDiscretised = argv[4];
          const std::string expressionOriginal = argv[5];
	const std::string genomeSizeFile = argv[6];
          const unsigned int window = 1000;
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
	std::cout<<"Extracting discretised gene expression information for "<<geneID<<std::endl;
	ExpressionReader expR(expressionDiscretised);
	expR.loadExpressionData(geneID,false);
	std::map<std::string, double> expressionMap;
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
	//Convert to genomic coordinates
	std::vector<std::pair<unsigned int, unsigned int> > genomeConv= sp.convertSegmentationToGenomicCoordinates(segments,genomicCoordinates);
	
	std::cout<<"Segmentation into "<<segments.size()<<" bins completed."<<std::endl;
	//Print binning
	//	for (auto& element : segments){
	//		std::cout<<"["<<element.first<<","<<element.second<<")"<<std::endl;
	//	}

	//Print binning
//	for (auto& element : genomeConv){
//		std::cout<<"["<<element.first<<","<<element.second<<"]"<<std::endl;
//	}

	//Loading original gene expression data
	std::cout<<"Extracting original gene expression information for "<<geneID<<std::endl;
	ExpressionReader expO(expressionOriginal);
	expO.loadExpressionData(geneID,false);
	std::map<std::string, double> expressionMapO = expO.getExpressionMap();

	//Compute signal in segmented regions
	std::cout<<"Computing mean signal per bin and sample"<<std::endl;
	BinSelection bs = BinSelection(bigWigPath,expressionMapO);
	bs.computeMeanSignal(std::get<0>(genomicCoordinates), genomeConv);
	std::cout<<"Computing correlation between signal and gene expression"<<std::endl;
	//Assess correlation of signal in bins to gene expression
	std::vector<std::pair<double,double> > corP = bs.computePearsonCorrelation();
//	std::cout<<"Pearson based correlation"<<std::endl;
//	for (auto& element : corP){
//		std::cout<<element.first<<" "<<element.second<<std::endl;
//	}
//	std::cout<<"Spearman based correlation"<<std::endl;
	std::vector<std::pair<double,double> > corS = bs.computeSpearmanCorrelation();
//	for (auto& element : corS){
//		std::cout<<element.first<<" "<<element.second<<std::endl;
//	}
	std::cout<<"Start	End	Pearson	pValue	Spearman	pValue"<<std::endl;
	for (unsigned int i=0; i<genomeConv.size();i++){
		std::cout<<genomeConv[i].first<<"	"<<genomeConv[i].second<<"	"<<corP[i].first<<"	"<<corP[i].second<<"	"<<corS[i].first<<"	"<<corS[i].second<<std::endl;
	}
	//Generate a txt file with DNase signal and gene expression across sample for the gene of interest including sample IDs and genomic location
          if (verbose){
		std::cout<<gtf<<std::endl;
		std::cout<<gsr<<std::endl;
		std::cout<<expR<<std::endl;
		std::cout<<expO<<std::endl;
		std::cout<<SPIG<<std::endl;	
		std::cout<<bs<<std::endl;	
	}
	
	return 0;
} 
