#include <assert.h>
#include <vector>
#include <cstdio>
#include <iostream>
#include <limits>
#include <algorithm>
#include <utility> // needed for pair
#include <cmath> // needed for abs and log2

#include "Segment.h"

Binning::Binning(Data& d, bool v){
	data = d;
	verbose = v;
	// determine number of columns
	length = d.n;
}

double Binning::modelCost(int bins, int total) const{
	assert(bins > 0);
	double weightCutpoints = 0.0;
	if(total == length)
		weightCutpoints = S::log2nChoosek(total - 1, (bins - 1));
	double w_bins = S::log2N(bins);
	double means = 0.0;
	if(total == length) 
		means = data.categoryCount() * bins * log2(std::abs(data.max - data.min) / (1.0 / (double) data.resolution));
	return w_bins + means + weightCutpoints;
}


void Binning::printSegmentation(S::segmentation seg) const {
    std::cout << "Number of bins: " << seg.size() << " --> ";
    for (unsigned long i = 0; i < seg.size(); i++) {
        std::cout << "[" << seg[i].start << ", " << seg[i].end << ") ";
    }
    std::cout << std::endl;
}

std::vector<int> Binning::Merge(std::vector<int>& a, std::vector<int>& b) const {
    std::vector<int> newone;
    unsigned long n = 0;
    unsigned long r = 0;
    while( n < a.size() || r < b.size()){
        if (r == b.size()) {
            newone.push_back(a[n]);
            n++;
        }else if(n == a.size()) {
            newone.push_back(b[r]);
            r++;
        }else if (a[n] == b[r]) {
            newone.push_back(a[n]);
            n++;
            r++;
        }else if (a[n] < b[r]) {
            newone.push_back(a[n]);
            n++;
        }else{
            newone.push_back(b[r]);
            r++;
        }
    }
    return newone;
}

void Binning::applyDPFlexi(int kbins, Fraction& f,bool exactK) const{
    int beta = (int) f.seg.size();
    int fSize = f.size();
    // precompute weights
    S::weights w = data.precomputeWeightsGaussFast(f);
    

    int depth = kbins <= 1 ? beta : std::min(kbins, beta);
    
    // init result vector
    std::vector< S::binning > results(beta);
    std::vector< double > compressionCosts(beta);
    S::binning bestBinning;

    // fill first row
    for(int i = 0; i < beta; i++){
        S::interval v;
        if(i == 0){
            v = (f.seg)[i];
        }else{
            v = results[i-1].bins[0];
            v.joinWith((f.seg)[i]);
        }
        results[i].bins.push_back(v);
        results[i].totalCosts = w[i][0];
    }
    // determine costs
    bestBinning = results[beta-1];
    double mc1 = modelCost(1, fSize);
    bestBinning.totalCosts += mc1;
    f.initial = bestBinning.totalCosts;
    if(verbose)
        std::cout << "Size: " << (1) << " L(D|M): " << results[beta-1].totalCosts << " L(M): " << mc1 << std::endl;
    // fill rest
    for(int l = 1; l < depth; l++){
        std::vector< S::binning > new_row(beta);
        for(int i = l; i < beta; i++){
            // find max pos (minimizing L(D|M))
            int pos = 0;
            double best = std::numeric_limits<double>::max();
            S::interval best_iv;
            for(int j = l-1; j < i; j++){
                S::interval rest_iv;
                rest_iv.start = (f.seg)[j+1].start;
                rest_iv.end = (f.seg)[i].end;
                double score = results[j].totalCosts + w[i][j+1];
                if(score < best){
                    pos = j;
                    best = score;
                    best_iv = rest_iv;
                }
            }
            // set entry
            new_row[i].totalCosts = best;
            S::segmentation prev_bins = results[pos].bins;
            prev_bins.push_back(best_iv);
            new_row[i].bins = prev_bins;
        }
        // update costs and last row
        results = new_row;
        // calculate costs
        double currMC = modelCost(l+1, fSize);
        double new_costs = results[beta-1].totalCosts + currMC;
        if(verbose)
            std::cout << "Size: " << (l+1) << " L(D|M): " << results[beta-1].totalCosts << " L(M): " << currMC << std::endl;
        if(new_costs < bestBinning.totalCosts){
            bestBinning = results[beta-1];
            bestBinning.totalCosts = new_costs;
        }
    }
    // find best
    if(kbins > 0 && exactK){
        bestBinning = results[beta-1];
    }
    f.compressed = bestBinning.totalCosts;
    double compression = f.compressed / f.initial * 100;
    if (verbose){
    std::cout << "----------------RESULTS----------------\n";
    printf("The relative compression is %.4f%%.\n", compression);
}
    f.updateSegmentation(bestBinning.bins);
}

void Binning::runSPAN(int k, Fraction& fraction, bool exactK) const{
	if(!fraction.merged){
        // apply initial equal width binning if stepSize > 1; otw. every feature is in its own bin
		std::vector<int> merge;
		for(int i = fraction.start; i < fraction.end; i += fraction.stepSize){
			merge.push_back(i);
		}
		int pos = 0;
		S::interval actual;
		for(int i = 0; i < (int) merge.size(); i++){
			if (pos == 0) {
				actual.start = merge[i];
			}else{
				actual.end = merge[i];
			}
			pos++;
			if (pos == 2) {
				pos = 1;
				fraction.seg.push_back(actual);
				actual.start = actual.end;
			}
		}
		if((fraction.seg)[fraction.seg.size()-1].end != fraction.end){
			actual.start = actual.end;
			actual.end = fraction.end;
			fraction.seg.push_back(actual);
		}
	} // if merged segments already exist
    
	bool autorun = k <= 1;
	if (verbose){
		std::cout << "Segment-count before execution: " << fraction.seg.size() << " (initial step size=" << fraction.stepSize << ") \n";
		if (autorun) {
			std::cout << "--> Use MDL to select number of segments. \n";
		} else {
		std::cout << "--> Don't trust MDL and find k=" << k << " segments. \n";
		}
	}
	applyDPFlexi(k,fraction,exactK);
}
