#include "Span.h"

//Adapted version of the SPAN code from Alex which generates the binning and returns it as an Fraction object
std::vector<std::pair<unsigned int, unsigned int> > SPAN::runSpan(Data& d, int s = 1, unsigned int maxCores = 1){
    int k = 0;
    bool useClasses = true;
    bool useCRE = false;
    char method = 'g';
    bool verbose = false;
    bool ignoreLabel = false;
    int parts = 1;
    
    //S::costfunptr ptr;
    //ptr = gaussian;
    Binning bins(d, verbose);
    // set up fractions
    std::vector<Fraction*> fractions;
    parts=round(log2(d.n)); //Added by fschmidt
    fractions.resize(parts); //Added by fschmidt
    std::cout<<"Spitting into "<<parts<<" subproblems"<<std::endl;
    if(parts > 1){
        omp_set_num_threads(fmin(parts,maxCores)); //Added by fschmidt
        int stdwidth = floor((double)d.n / (double)parts);
        #pragma omp parallel for
        for(int i = 0; i < parts; i++){
            int begin = i * stdwidth;
            int end = i == (parts-1) ? d.n : (i+1) * stdwidth;
            Fraction* f = new Fraction(begin, end, s, useClasses);
            bins.runSPAN(k, f);
  	  //fractions.push_back(f);
	  fractions[i]=f; //Added by fschmidt
    }
   }
	
	std::cout<<"Merging into a second layer"<<std::endl;
	//Test if a second division layer to merge smaller cluster improves the required time.
	std::vector<Fraction*> second_layer;
	Fraction* current = NULL;
	for(int i =0; i < parts; i++){
		if(i % 5 == 0){ // oder welche Bedingung auch immer du magst
			if(i != 0){
				second_layer.push_back(current);
			}
			current = fractions[i];
		}else{
			current->mergeIn(fractions[i]);
		}
	}
	second_layer.push_back(current);  // last element
	std::cout<<"Second layer filled with "<<second_layer.size()<<std::endl;
	// algorithmus ausführen
	#pragma omp parallel for
	for(int i = 0; i < second_layer.size(); i++){
		bins.runSPAN(k, second_layer[i]);
	}

	std::cout<<"Merging results of second layer"<<std::endl;
        Fraction* fAll = second_layer.size() > 1 ? second_layer[0] : new Fraction(0, d.n, s, useClasses);
    for(int i = 1; i < second_layer.size(); i++){
        fAll->mergeIn(second_layer[i]);
    }
    //END second division layer

   //Generate vector of tupels holding final binning
   // run last run
    fprintf(stdout, "Final run:\n");
    bins.runSPAN(k, fAll);
	
	std::vector<std::pair<unsigned int, unsigned int> > resultVector;
	for(int i = 0; i < fAll->seg->size(); i++){
		std::cout<<"Adding "<<(*(fAll->seg))[i].start + 1<<", "<<(*(fAll->seg))[i].end<<std::endl;
 		resultVector.push_back(std::make_pair((*(fAll->seg))[i].start + 1,(*(fAll->seg))[i].end));
 	}

    if(parts > 1){
        for(int i = 0; i < parts; i++){
            delete fractions[i];
        }
    }

	return resultVector;
}

