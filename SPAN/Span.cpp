#include <getopt.h> //options
#include <stdlib.h> //set handler
#include <cstdio>
#include <fstream>
#include <cmath>
#include <string>
#include <omp.h> //Added by fschmidt
#include "Wall_Time.h"
#include "Data.h"
#include "Segment.h"

/**
 * Main class; reads params and executes Span.
 **/

void outofmemory(){
    std::cerr << "Oh noes! out of memory\n";
    exit(1);
}

double inline gaussian(double sum, double count){
    if(count <= 0){
        return 0;
    }
    return - sum * sum / (2 * count * log(2.0f));
}

double inline poisson(double sum, double count){
    if(count <= 0){
        return 0;
    }
    double fr = sum / count;
    double res = -sum / log(2);
    if(sum > 0){
        res += sum * log2(fr);
    }
    return -res;
}



//Adapted version of the SPAN code from Alex which generates the binning and returns it as an Fraction object
Fraction runSpan(Data& d, int s = 1, unsigned int maxCores = 1){
    int k = 0;
    bool useClasses = true;
    bool useCRE = false;
    char method = 'g';
    bool verbose = false;
    bool ignoreLabel = false;
    int parts = 1;
    std::set_new_handler(outofmemory);
    
    S::costfunptr ptr;
    ptr = gaussian;
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
	
	std::cout<<"Merging into a sceond layer"<<std::endl;
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

    // run last run
    fprintf(stdout, "Final run:\n");
    bins.runSPAN(k, fAll);
    if(parts > 1){
        for(int i = 0; i < parts; i++){
            delete fractions[i];
        }
    }
   return (*fAll);
     
    // delete all
    // delete fAll;
}


void printFraction(Fraction& F, std::string outname=""){
	if (outname=="")
		throw std::invalid_argument("No file name for results specified");
	if (&F == NULL)
		throw std::invalid_argument("No Fraction object provided");
    // write output
        std::ofstream myfile;
        myfile.open(outname);
        if(myfile.is_open()){
            for(int i = 0; i < F.seg->size(); i++){
                myfile << std::to_string((*(F.seg))[i].start + 1) + "\t" + std::to_string((*(F.seg))[i].end) + "\n";
            }
        }
        myfile.close();
  }
