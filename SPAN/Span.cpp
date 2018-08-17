#include "Span.h"

//Adapted version of the SPAN code from Alex which generates the binning and returns it as a vector of tupels
std::vector<std::pair<unsigned int, unsigned int> > SPAN::runSpan(Data& d, int s = 1, unsigned int maxCores = 1,bool verbose = false,unsigned int stitchSegmentLength = 5000){
	int k = 0;
	int parts = 1;
	int stdwidth = stitchSegmentLength;
	Binning bins(d, verbose);
	// set up fractions
	std::vector<Fraction> fractions;
	parts=ceil(d.n/float(stdwidth));//round(log2(d.n)); //Added by fschmidt
	fractions.resize(parts); //Added by fschmidt
    // track the costs for fractions

    int numberOfInitialBins = 0;
    double compressedScore= 0.0;
    double initialScore = 0.0;
	if (verbose){
		std::cout<<"Spitting into "<<parts<<" subproblems"<<std::endl;
	}	
	if(parts > 1){
		omp_set_num_threads(fmin(parts,maxCores)); //Added by fschmidt
		#pragma omp parallel for
		for(int i = 0; i < parts; i++){
			int begin = i * stdwidth;
			int end = i == (parts-1) ? d.n : (i+1) * stdwidth;
			Fraction f = Fraction(begin, end, s, d.getCategories());
			bins.runSPAN(k, f, false);
	  		fractions[i]=f; //Added by fschmidt
		          numberOfInitialBins += f.initialBins;
		          initialScore += f.initial;//Data;
		 	compressedScore += f.compressed;
		}
	}else{
		Fraction f = Fraction(0, 0 == (parts-1) ? d.n : stdwidth, s, d.getCategories());
          bins.runSPAN(k, f, false);
          fractions[0]=f; //Added by fschmidt
	numberOfInitialBins = f.initialBins;
          initialScore = f.initial;//Data;
	compressedScore = f.compressed;
	}

	Fraction fAll = fractions[0];	
	for(unsigned int i = 1; i < fractions.size(); i++){
		fAll.mergeIn(fractions[i]);
		}
	double compressionRatio = compressedScore/initialScore; //finalCosts / initialCosts;

	if (verbose){
		std::cout<<"Compression ratio: "<<compressionRatio<<std::endl;// "(from  "<<initialCosts << " to "<< finalCosts<<")"<<std::endl;
	}
	std::vector<std::pair<unsigned int, unsigned int> > resultVector;
	for(unsigned int i = 0; i < fAll.seg.size(); i++){
		resultVector.push_back(std::make_pair((fAll.seg)[i].start + 1,(fAll.seg)[i].end));
	}
	return resultVector;
}

std::vector<std::pair<unsigned int, unsigned int> > SPAN::convertSegmentationToGenomicCoordinates(std::vector<std::pair<unsigned int, unsigned int> >& segmentation, std::tuple<std::string, unsigned int, unsigned int,std::string> genomicCoordinates){
	std::vector<std::pair<unsigned int, unsigned int> > transformedSegments;
	unsigned int correctionPosition = std::get<1>(genomicCoordinates)-1;
	for (auto& posPair : segmentation){
		transformedSegments.push_back(std::make_pair(posPair.first+correctionPosition,posPair.second+correctionPosition));
	}
	return transformedSegments;
}
