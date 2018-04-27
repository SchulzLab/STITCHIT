#include "../BigWig/bigWig.h"
#include "BinSelection.h"
#include "CorComp.h"

BinSelection::BinSelection(const std::string& path,std::map<std::string,double>& expressionMap)
	:path_(path), expressionMap_(expressionMap)
{
	meanSignal_.resize(0);

}



//IMPORTANT NOTE: This method transforms the first coordinate to the 0-based index used in bigWig.h. As bigWig.h works with half open intervalls, but SPAN provides closed intervalls, the second coordinate is not changed.
void BinSelection::computeMeanSignal(std::string& chrom, const std::vector<std::pair<unsigned int, unsigned int> >& segments){
	 //Generating per base input for SPAN
	if(bwInit(1<<17) != 0) throw std::runtime_error("Error occured in bwInit.");
	//Iterating through all files in the specified directory
	if (boost::filesystem::is_directory(path_)){
		for(auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(path_), {})){
			double *stats = NULL;
			bigWigFile_t *fp = NULL;
			std::string filenameS = entry.path().string();
			char * filename= new char[filenameS.size() +1];
			char * chromC= new char[chrom.size()+1];
			std::copy(filenameS.begin(),filenameS.end(), filename);
			std::copy(chrom.begin(),chrom.end(), chromC);
			filename[filenameS.size()]='\0';
			chromC[chrom.size()]='\0';
			if (entry.path().extension().string() != ".bw") throw std::invalid_argument(filenameS+" is not a biwWig file. Expected file type is .bw");
			//Check for presence of expression data for the current sample, otherwise ignore the sample!
			std::vector<double> currentSample;
			fp = bwOpen(filename, NULL, "r");
			if(!fp) {
				throw std::invalid_argument("An error occured while opening the bw files");
			}
			//Get an example statistic - standard deviation
			//We want ~1 bins in the range
			sampleNames_.push_back(entry.path().stem().string());
			for (auto& seg : segments){
			stats = bwStats(fp, chromC, seg.first-1, seg.second, 1, mean);
			if(stats) {
					currentSample.push_back(stats[0]);
				}
			}
			free(stats);
			bwClose(fp);
			meanSignal_.push_back(currentSample);
		}
	}
	bwCleanup();
}

//IMPORTANT NOTE: This method transforms the first coordinate to the 0-based index used in bigWig.h. As bigWig.h works with half open intervalls, but SPAN provides closed intervalls, the second coordinate is not changed.
void BinSelection::computeMeanSignal(const std::vector<std::tuple<std::string, unsigned int, unsigned int> >& segments){
	 //Generating per base input for SPAN
	if(bwInit(1<<17) != 0) throw std::runtime_error("Error occured in bwInit.");
	//Iterating through all files in the specified directory
	if (boost::filesystem::is_directory(path_)){
		for(auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(path_), {})){
			double *stats = NULL;
			bigWigFile_t *fp = NULL;
			std::string filenameS = entry.path().string();
			char * filename= new char[filenameS.size() +1];
			std::copy(filenameS.begin(),filenameS.end(), filename);
			filename[filenameS.size()]='\0';
			if (entry.path().extension().string() != ".bw") throw std::invalid_argument(filenameS+" is not a biwWig file. Expected file type is .bw");
			//Check for presence of expression data for the current sample, otherwise ignore the sample!
			std::vector<double> currentSample;
			fp = bwOpen(filename, NULL, "r");
			if(!fp) {
				throw std::invalid_argument("An error occured while opening the bw files");
			}
			//We want ~1 bins in the range
			sampleNames_.push_back(entry.path().stem().string());
			for (auto& seg : segments){
				std::string chrom = std::get<0>(seg);
				char * chromC= new char[chrom.size()+1];
				std::copy(chrom.begin(),chrom.end(), chromC);
				chromC[chrom.size()]='\0';
				stats = bwStats(fp, chromC, std::get<1>(seg)-1, std::get<2>(seg), 1, mean);
				if(stats) {
					currentSample.push_back(stats[0]);
				}
			}
			free(stats);
			bwClose(fp);
			meanSignal_.push_back(currentSample);
		}
	}
	bwCleanup();
}


std::vector<double> BinSelection::getExpressionVectorByNames(){
	std::vector<double> expVecTemp;
	for (auto& sample : sampleNames_){
		expVecTemp.push_back(expressionMap_[sample]);
	}
	return expVecTemp;
}

std::vector<double> BinSelection::getSignalVectorBySegment(unsigned int segID){
	std::vector<double> segmentSignal;
	for (unsigned int sam = 0; sam < meanSignal_.size(); sam++){
		segmentSignal.push_back(meanSignal_[sam][segID]);
	}			
	return segmentSignal;
}

std::vector<std::pair<double,double> > BinSelection::computePearsonCorrelation(){
	std::vector<std::pair<double, double> > correlation;
	unsigned int numSeg = meanSignal_[0].size();
	std::vector<double> expVec=getExpressionVectorByNames();
	for (unsigned int seg = 0; seg < numSeg; seg++){
		std::vector<double> signalVec = getSignalVectorBySegment(seg);
		CorComp cC(expVec,signalVec);
		double cor = cC.computePearsonCorrelation();
		double pValue = cC.getPvalue(cor);
		correlation.push_back(std::make_pair(cor,pValue));
	}
	return correlation;
}

std::vector<std::pair<double, double> > BinSelection::computeSpearmanCorrelation(){
	std::vector<std::pair<double, double> > correlation;
	unsigned int numSeg = meanSignal_[0].size();
	std::vector<double> expVec=getExpressionVectorByNames();
	for (unsigned int seg = 0; seg < numSeg; seg++){
		std::vector<double> signalVec = getSignalVectorBySegment(seg);
		CorComp cC(expVec,signalVec);
		double cor = cC.computeSpearmanCorrelation();
		double pValue = cC.getPvalue(cor);
		correlation.push_back(std::make_pair(cor,pValue));
	}
	return correlation;
}



std::vector<std::vector<double> >& BinSelection::getMeanSignal(){
	return meanSignal_;
}

std::vector<std::string>& BinSelection::getSampleNames(){
	return sampleNames_;
}

/*! \brief Brief description.
 *
 *  Detailed description starts here.
 * @param
 * @param
 */
std::ostream& operator<<(std::ostream& os, const BinSelection& r){
	unsigned int numSam = r.meanSignal_.size();
	std::vector<double> expVecTemp;
	for (auto& sample : r.sampleNames_){
		expVecTemp.push_back(r.expressionMap_[sample]);
	}
	for (unsigned int sam = 0; sam < numSam; sam++){
		os << r.sampleNames_[sam];
		for (double element : r.meanSignal_[sam]){
			os << "	" << element;
		}	
	os <<"	"<< expVecTemp[sam]<< '\n';
	}
	return os;
}

void BinSelection::storeSignificantSignal(const std::string& filename, float threshold, std::vector<std::pair<double,double> > correlation,  std::vector<std::pair<unsigned int, unsigned int> > intervalPosition,std::tuple<std::string, unsigned int, unsigned int,std::string> genomePos){
	unsigned int numSam = meanSignal_.size();
	std::vector<double> expVecTemp;
	for (auto& sample : sampleNames_){
		expVecTemp.push_back(expressionMap_[sample]);
	}
	std::vector<unsigned int> validIndex;
	for (unsigned int index = 0; index< correlation.size(); index++){
		if (correlation[index].second <= threshold){
			validIndex.push_back(index);
		}
	}

	if (not(validIndex.empty())){
		std::ofstream outfile;
		outfile.open(filename);	
		outfile<<"\t";
		for (auto& item : validIndex){
			outfile <<std::get<0>(genomePos)<<":"<<intervalPosition[item].first << "-" <<intervalPosition[item].second<<"\t";
		}
		outfile<<"Expression\n";
		for (unsigned int sam = 0; sam < numSam; sam++){
			outfile << sampleNames_[sam];
			for (auto& item : validIndex){
				outfile << "\t" <<  meanSignal_[sam][item];
			}	
		outfile <<"\t"<< expVecTemp[sam]<< '\n';
		}
		outfile.close();
	}
}

void BinSelection::storeSignificantSignal(const std::string& filename, float threshold, std::vector<std::pair<double,double> > correlation,  std::vector<std::tuple<std::string, unsigned int, unsigned int> > intervalPosition,std::tuple<std::string, unsigned int, unsigned int,std::string> genomePos){
	unsigned int numSam = meanSignal_.size();
	std::vector<double> expVecTemp;
	for (auto& sample : sampleNames_){
		expVecTemp.push_back(expressionMap_[sample]);
	}
	std::vector<unsigned int> validIndex;
	for (unsigned int index = 0; index< correlation.size(); index++){
		if (correlation[index].second <= threshold){
			validIndex.push_back(index);
		}
	}

	if (not(validIndex.empty())){
		std::ofstream outfile;
		outfile.open(filename);	
		outfile<<"\t";
		for (auto& item : validIndex){
			outfile <<std::get<0>(genomePos)<<":"<<std::get<1>(intervalPosition[item]) << "-" <<std::get<2>(intervalPosition[item])<<"\t";
		}
		outfile<<"Expression\n";
		for (unsigned int sam = 0; sam < numSam; sam++){
			outfile << sampleNames_[sam];
			for (auto& item : validIndex){
				outfile << "\t" <<  meanSignal_[sam][item];
			}	
		outfile <<"\t"<< expVecTemp[sam]<< '\n';
		}
		outfile.close();
	}
}
	
