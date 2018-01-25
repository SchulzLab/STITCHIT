#include "GTFReader.h"

GTFReader::GTFReader(const std::string& gtfFileName, std::map<std::string, int>& genomeSize, const int& window)
	:gtfFileName_(gtfFileName),
          genomeSize_(genomeSize),                                                                                                                                                                                                             
          window_(window)
{
	genomicLocation_= std::make_tuple("chr0",0,0,"NA");
}

/*! \brief Brief description.
 *         Brief description continued.
 *  Detailed description starts here.
 * @return
 * @throw
 */
void GTFReader::findGenomicLocation(const std::string& targetGeneID){
	int pos1, pos2;
	std::string chromosome, buf, temp, strand, convertedGeneID, geneID;
	std::ifstream annotationFile;
	const int zero = 0;
	annotationFile.open(gtfFileName_);
	if (!annotationFile) throw std::invalid_argument("Gene annotation file "+gtfFileName_+" could not be loaded");
	
	bool flag=false;
	while(!annotationFile.eof()){
		std::getline(annotationFile,temp);
		std::stringstream ss(temp); // Insert the string into a stream
		if (temp != ""){
			if (ss >> chromosome >> buf >> buf){
				if (buf == "gene"){
					if (ss >> pos1 >> pos2 >> buf >> strand >> buf >> buf >> geneID){
						convertedGeneID=geneID.substr(1,15);    
						if (convertedGeneID==targetGeneID){
							pos1=std::max(zero,pos1-window_);
							pos2=std::min(genomeSize_[chromosome],pos2+window_);
							genomicLocation_=std::make_tuple(chromosome, pos1, pos2, strand);
							flag=true;
							break;
						}
					}	
					else{
						throw std::invalid_argument("Annotation file "+gtfFileName_+" does not match specifications");
					}
				}
			}
			else {
				throw std::invalid_argument("Annotation file "+gtfFileName_+" does not match specifications");
			}	
		}
	}          
	if (!flag){         
		annotationFile.close();
		//Exit if the gene could not be found
		throw std::invalid_argument("The provided geneID "+targetGeneID+" could not be found in the annotation file");
	}
}


const std::string& GTFReader::getGTFfileName(){
	return gtfFileName_;
}

const std::tuple<std::string,unsigned int, unsigned int,std::string>& GTFReader::getGenomicLocation(){
	return genomicLocation_;
}

const int& GTFReader::getWindow(){
	return window_;
}

void GTFReader::setWindow(int window){
	window_ = window;
}


std::ostream& operator<<(std::ostream& os, const GTFReader& r){
	os << std::get<0>(r.genomicLocation_) << " " << std::get<1>(r.genomicLocation_) << " " << std::get<2>(r.genomicLocation_) <<std::endl;
	return os;
}
