#include "SPANInputGenerator.h"
#include "GenomeSizeReader.h"
#include "GTFReader.h"
#include "ExpressionReader.h"
#include "BinSelection.h"
#include "../SPAN/Data.h"
#include "../SPAN/Span.h"
#include "boost/program_options/options_description.hpp"
#include "boost/program_options/variables_map.hpp"
#include "boost/program_options/parsers.hpp"

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
          std::string bigWigPath;
          std::string annotationFile;
          std::string geneID;
          std::string expressionDiscretised;
          std::string expressionOriginal;
	std::string genomeSizeFile;
	std::string corM;
	unsigned int window;
	unsigned int stepSize;
	unsigned int maxCores;
	float pvalue;
	bool verbose;
	
	boost::program_options::options_description desc("STITCH options");
	desc.add_options()
		("help","Show help message")
		("bigWigPath,b", boost::program_options::value<std::string>(&bigWigPath), "Path to big wig files")
		("annotationFile,a",boost::program_options::value<std::string>(&annotationFile), "Path to the annotation file that should be used")
		("geneID,g",boost::program_options::value<std::string>(&geneID), "ID of the gene that should be segmented")
		("expressionDiscretised,d",boost::program_options::value<std::string>(&expressionDiscretised), "File containing the discretised expression information")
		("expressionOriginal,o",boost::program_options::value<std::string>(&expressionOriginal), "File containng the original expression information")
		("genomeSizeFile,s",boost::program_options::value<std::string>(&genomeSizeFile), "File containig the size of the chromosomes of the reference genome")
		("window,w",boost::program_options::value<unsigned int>(&window)->default_value(5000), "Size of the window considered upstream and downstream of the gene body")
		("stepSize,z",boost::program_options::value<unsigned int>(&stepSize)->default_value(1), "Resolution parameter used by SPAN to merge the initial binning")
		("cores,c", boost::program_options::value<unsigned int>(&maxCores)->default_value(2), "Number of cores to be used for the computation")
		("pvalue,p",boost::program_options::value<float>(&pvalue)->default_value(0.05),"Signifcance level for the correlation of a segment to be considered in the output")
		("correlationMeasure,m",boost::program_options::value<std::string>(&corM)->default_value("Pearson"),"Method used to compute correlation between expression and epigenetic signal. Can be Pearson(default) or Spearman")
		("verbose,v", boost::program_options::value<bool>(&verbose)->default_value(false), "True if additional status reports should be generated, false otherwise")
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

	if (not(vm.count("bigWigPath"))){
		std::cout<<"--bigWigPath is not specified. Please provide the path to a folder containing one big wig file per sample"<<std::endl;;		
		return 1;
	}

	if (not(vm.count("annotationFile"))){
		std::cout<<"--annotationFile is not specified. Please provide the path to a gene annotation file in gtf format"<<std::endl;;		
		return 1;
	}

	if (not(vm.count("geneID"))){
		std::cout<<"--geneID is not specified. Please provide the geneID of the gene that should be segmented."<<std::endl;;		
		return 1;
	}

	if (not(vm.count("expressionDiscretised"))){
		std::cout<<"--expressionDiscretised is not specified. Please provide the path to a file containing discretised expression data with genes in the rows and samples in the columns. Note that multiclass classification is supported."<<std::endl;;		
		return 1;
	}

	if (not(vm.count("expressionOriginal"))){
		std::cout<<"--expressionOriginal is not specified. Please provide the path to a folder containing one big wig file per sample"<<std::endl;;		
		return 1;
	}

	if (not(vm.count("genomeSizeFile"))){
		std::cout<<"--genomeSizeFile is not specified. Please provide a file holding the total length per chromosome."<<std::endl;;		
		return 1;
	}

	if (vm.count("correlationMeasure")){
		if ((corM != "Pearson") and (corM != "Spearman")){
			std::cout<<"Please specify a valid correlation measure: Pearson or Spearman"<<std::endl;;		
			return 1;
		}
	}

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
	std::cout<<"Coordinates found: "<<std::get<0>(genomicCoordinates)<<" "<<std::get<1>(genomicCoordinates)+window<<" "<<std::get<2>(genomicCoordinates)-window<<", gene length: "<<std::get<2>(genomicCoordinates)-std::get<1>(genomicCoordinates)<<std::endl;
	
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
	input.setData(perBaseInputData,true,'g',stepSize,false);

	//Feed input to SPAN and execute
	std::cout<<"Segmentation in progress..."<<std::endl;
	SPAN sp= SPAN();
	std::vector<std::pair<unsigned int, unsigned int> > segments = sp.runSpan(input,stepSize,maxCores,verbose);
	//Convert to genomic coordinates
	std::vector<std::pair<unsigned int, unsigned int> > genomeConv= sp.convertSegmentationToGenomicCoordinates(segments,genomicCoordinates);
	
	std::cout<<"Segmentation into "<<segments.size()<<" bins completed."<<std::endl;

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
	std::vector<std::pair<double,double> > corS = bs.computeSpearmanCorrelation();
	std::cout<<"Start	End	Pearson	pValue	Spearman	pValue"<<std::endl;
	for (unsigned int i=0; i<genomeConv.size();i++){
		std::cout<<genomeConv[i].first<<"	"<<genomeConv[i].second<<"	"<<corP[i].first<<"	"<<corP[i].second<<"	"<<corS[i].first<<"	"<<corS[i].second<<std::endl;
	}

	//Generate a txt file with DNase-seq signal and gene expression across samples for the gene of interest in the significant segments including sample IDs and genomic location
	if (corM=="Pearson"){
		bs.storeSignificantSignal("Segmentation_"+geneID+".txt",pvalue,corP,genomeConv);
	}else{
		bs.storeSignificantSignal("Segmentation_"+geneID+".txt",pvalue,corS,genomeConv);
	}
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
