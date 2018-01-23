#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <assert.h>
#include <utility>
#include <limits>

#include "Data.h"

#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif

/**
 * This file includes the method to read the data and several scoring methods.
 **/

unsigned int Data::getRowCount(){
	return data.size();
}

unsigned int Data::getColCount(){
	return data[0].size();
}

double& Data::getElement(unsigned int row,unsigned int col){
	return data[row][col];
}

std::vector<double>& Data::getRow(unsigned int row){
	return data[row];
}

void Data::setData(std::vector<std::vector<double> > input, bool binaryClassifierInLastRow, bool comulat, char meth, int s,  bool ignoreLabel){
	commulative = comulat;
	method = meth;
	eval = NULL;
	stepSize = s;
	log2sqrt2pi = log2(std::sqrt(2 * 3.14159265358979323846));
	logbase = log(2);
	resolution = 10;
	data=input;

	m = data.size();
	n = data[0].size();
	
	min = std::numeric_limits<double>::max();
	max = -min;
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			double v = data[j][i];
			if(v < min){
				min = v;
				}
			if(v > max){
				max = v;
				}
		}
	}
	resolution = 10;
}



//This function assumes that no sample IDs are provided
void Data::read_file(const char* filename, bool binaryClassifierInLastRow, bool comulat, char meth, int s,  bool ignoreLabel){
    // set commulative
    commulative = comulat;
    method = meth;
    eval = NULL;
    stepSize = s;
    log2sqrt2pi = log2(std::sqrt(2 * 3.14159265358979323846));
    logbase = log(2);
    
    std::ifstream f;
    std::string line;
    
    f.open(filename);
    if (f.is_open()) {
        int counter = 0;
        while (std::getline(f, line)) {
            std::vector<double> new_vec;
            std::istringstream is(line);
            double actual;
            while (is >> actual) {
                new_vec.push_back(actual);
            }
            if(binaryClassifierInLastRow || ignoreLabel){
                int classifier = new_vec[new_vec.size() - 1];
                if(classifier == 0 && !ignoreLabel){
                    cl0.push_back(counter);
                }else{
                    cl1.push_back(counter);
                }
                new_vec.pop_back();
            }else{
                cl1.push_back(counter);
            }
            all.push_back(counter);
            data.push_back(new_vec);
            counter++;
        }
    }
    
    // Set n (cols) and m (rows) and check if each row has the same amount of columns.
    m = data.size();
    if ( m > 0 ){
        n = data[0].size();
    } else {
        f.close();
       throw std::invalid_argument("Invalid number of rows");
	}
    for ( int i = 0; i < m; i++) {
        if ( data[i].size() != (unsigned long) n) {
                    f.close();
		throw std::invalid_argument("Invalid number of columns");
        }
    }
    min = std::numeric_limits<double>::max();
    max = -min;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            double v = data[j][i];
            if(v < min){
                min = v;
            }
            if(v > max){
                max = v;
            }
        }
    }
    
    resolution = 10;
    
    f.close();
}


Data::Data(const char* filename, bool binaryClassifierInLastRow, bool comulat, char m, int s, bool ignoreLabel){
    read_file(filename, binaryClassifierInLastRow, comulat, m, s, ignoreLabel);
}

sc_pair Data::getDensities(S::interval iv, int label){
    std::vector<int> rows;
    if(label == 0){
        rows = cl0;
    }else if(label == 1){
        rows = cl1;
    }else{
        rows = all;
    }
    double sum = 0;
    double count = (iv.end - iv.start) * rows.size();
    for(int i = iv.start; i < iv.end; i++){
        for(unsigned long k = 0; k < rows.size(); k++){
            int j = rows[k];
            sum += data[j][i];
        }
    }
    sc_pair ret;
    ret.first = sum;
    ret.second = count;
    return ret;
}

count_entry Data::getEntry(S::interval iv){
    bool binary = (cl0.size() > 0);
    count_entry ce;
    ce.init();
    if(binary){
        sc_pair p0 = getDensities(iv, 0);
        sc_pair p1 = getDensities(iv, 1);
        ce.ps0 = p0.first;
        ce.ps1 = p1.first;
        ce.pc0 = p0.second;
        ce.pc1 = p1.second;
    }else{
        sc_pair p = getDensities(iv, -1);
        ce.pc0 = 0;
        ce.ps0 = 0;
        ce.ps1 = p.first;
        ce.pc1 = p.second;
    }
    return ce;
}

omega* Data::extractDataMap(S::interval iv, int label){
    std::vector<int> rows;
    if(label == 0){
        rows = cl0;
    }else if(label == 1){
        rows = cl1;
    }else{
        rows = all;
    }
    omega* elementCount = new omega();
    for(int i = iv.start; i < iv.end; i++){
        for(unsigned long k = 0; k < rows.size(); k++){
            int j = rows[k];
            elementCount->operator[](data[j][i]) += 1;
        }
    }
    return elementCount;
}

double Data::getComplexity(count_entry e){
    double comp = eval(e.ps1, e.pc1);
    if(cl0.size() > 0){
        comp += eval(e.ps0, e.pc0);
    }
    return comp;
}


double Data::ModelWeight(S::interval iv, data_segment e, int i, int j){
    sc_pair p1;
    p1.first = e.ps1;
    p1.second = e.pc1;
    omega* o1 = e.data1;
    int N1 = cl1.size() * iv.length();
    if(cl0.size() > 0){
        sc_pair p0;
        p0.first = e.ps0;
        p0.second = e.pc0;
        omega* o0 = e.data0;
        int N0 = cl0.size() * iv.length();
        double c1 = 0;
        double c0 = 0;
        if(method == 'g'){
            c0 = Gauss(iv, p0, o0, N0, i, j);
            c1 = Gauss(iv, p1, o1, N1, i, j);
        }else if(method == 'l'){
            c0 = Laplace(iv, p0, o0, N0);
            c1 = Laplace(iv, p1, o1, N1);
        }else{
            c1 = SAE(iv, p1, o1);
            c0 = SAE(iv, p0, o0);
        }
        double score = (c1 / (double) cl1.size()) + (c0 / (double) cl0.size());
        return score;
    }else{
        int N = all.size() * iv.length();
        double c1 = 0;
        if(method == 'g'){
            c1 = Gauss(iv, p1, o1, N, i, j);
        }else if(method == 'l'){
            c1 = Laplace(iv, p1, o1, N);
        }else{
            double c1 = SAE(iv, p1, o1);
        }
        return (c1 / (double)all.size());
    }
}

double Data::GaussApprox(int i, int j, Fraction* f){
    double c1 = GaussApprox(i, j, 1, f) / (double) cl1.size();
    double c0 = 0.0;
    if(cl0.size() > 0){
        c0 = GaussApprox(i, j, 0, f) / (double) cl0.size();
    }
    double score = c1 + c0;
    return score;
}

double Data::GaussApprox(int i, int j, int label, Fraction* f){
    double sigma = 0.0;
    int ii = (j + 1);
    int jj = (i + 1);
    int diff = (*(f->seg))[i].end - (*(f->seg))[j].start;
    std::vector<double> css;
    std::vector<double> cs;
    int N = 0;
    if(label == 0){
        N = cl0.size();
        css = *(f->sse0);
        cs = *(f->se0);
    }else{
        N = cl1.size();
        css = *(f->sse1);
        cs = *(f->se1);
    }
    double myweight = (double) diff * (double) N;
    double sse = (css[jj] - css[ii-1]) - (1/myweight) * ((cs[jj] - cs[ii-1]) * (cs[jj] - cs[ii-1]));
    if(sse <= 0.0)
        return 0.0;
    sigma = std::sqrt(sse / (double) myweight);
    double score = (myweight/2) * log2(sigma*sigma*2.0*M_PI) + ((myweight/2.0)*log(2.0)) + N * log2(resolution);
    score = score < 0.0 ? 0.0 : score;
    return score;
}

double Data::Gauss(S::interval iv, sc_pair sc, omega* big_o, int N, int i, int j){
    double mean = (double) sc.first / (double) sc.second;
    mean = round(mean * resolution) / resolution;
    double sigma = 0;
    for(omega::iterator iter = big_o->begin(); iter != big_o->end(); ++iter){
        double a_key = iter->first;
        int a_value = iter->second;
        double diff = std::abs(mean - a_key);
        double term = (diff * diff) * a_value;
        sigma += term;
    }
    sigma = std::sqrt(sigma / (double) N);
    double l2sigma = sigma <= 1 ? 0 : log2(sigma); // if < 1 could be negative
    double score = N * (l2sigma + log2sqrt2pi + (1/(2 * logbase)));
    return score;
}

double Data::Laplace(S::interval iv, sc_pair sc, omega* big_o, int N){
    double mean = (double) sc.first / (double) sc.second;
    mean = round(mean * resolution) / resolution;
    double sigma = 0;
    for(omega::iterator iter = big_o->begin(); iter != big_o->end(); ++iter){
        double a_key = iter->first;
        int a_value = iter->second;
        double diff = std::abs(mean - a_key);
        double term = diff * a_value;
        sigma += term;
    }
    sigma = (2 * sigma) / N;
    double l2sigma = sigma <= 1 ? 0 : log2(sigma); // if < 1 could be negative
    double score = N * (l2sigma + (1/logbase));
    return score;
}

double Data::SAE(S::interval iv, sc_pair sc, omega* big_o){
    double mean = (double) sc.first / (double) sc.second;
    mean = round(mean * resolution) / resolution;
    double sae = 0;
    for(omega::iterator iter = big_o->begin(); iter != big_o->end(); ++iter){
        double a_key = iter->first;
        int a_value = iter->second;
        double diff = a_value * S::log2N(std::abs(mean - a_key) + 1); // mean == a_key ? 0 :
        sae += diff;
    }
    return sae;
}

void print(omega* omg){
    std::cout << "Start -----\n";
    for(omega::iterator iter = omg->begin(); iter != omg->end(); ++iter){
        double a_key = iter->first;
        int a_value = iter->second;
        std::cout << a_key << " -> " << a_value << std::endl;
    }
    std::cout << "----- End\n";
}

S::weights Data::precomputeWeightsSAEMap(S::segmentation seg){
    //precalculateCSSandCS(seg);
    int beta = (int) seg.size();
    bool classification = cl0.size() > 0;
    S::weights w;
    data_map_table t;
    for(int i = 0; i < beta; i++){
        std::vector< double > dummy(i+1);
        std::vector< data_segment > ds(i+1);
        w.push_back(dummy);
        t.push_back(ds);
    }
    // init diagonal
    for(int i = 0; i < beta; i++){
        t[i][i].set(getEntry(seg[i]));
        if(classification){
            t[i][i].set1(extractDataMap(seg[i], 1));
            t[i][i].set0(extractDataMap(seg[i], 0));
        }else{
            t[i][i].set1(extractDataMap(seg[i], -1));
            t[i][i].set0(0);
        }
        w[i][i] = ModelWeight(seg[i], t[i][i], i, i);
    }
    for(int i = 1; i < beta; i++){
        S::interval actual_iv = seg[i];
        for(int k = i - 1; k >= 0; k--){
            data_segment last;
            last.init();
            last.set1(new omega());
            if(classification){
                last.set0(new omega());
            }else{
                last.set0(0);
            }
            last.add(t[i][(k+1)]);
            last.add(t[k][k]);
            t[i][k] = last;
            actual_iv.joinWith(seg[k]);
            w[i][k] = ModelWeight(actual_iv, t[i][k], i, k);
            // free memory not needed anymore
            if(k == 0){
                delete t[i][k].data1;
                if(classification)
                    delete t[i][k].data0;
            }else if((k+1) != i){
                delete t[i][k+1].data1;
                if(classification)
                    delete t[i][k+1].data0;
            }
        }
    }
    // free memory
    for(int i = 0; i < beta; i++){
        delete t[i][i].data1;
        if(classification)
            delete t[i][i].data0;
    }
    return w;
}

S::weights Data::precomputeWeightsGaussFast(Fraction* f){
    precalculateCSSandCS(f);
    int beta = (int) f->seg->size();
    bool classification = f->binary;
    S::weights w;
    for(int i = 0; i < beta; i++){
        std::vector< double > dummy(i+1);
        w.push_back(dummy);
    }
    // init diagonal
    for(int i = 0; i < beta; i++){
        w[i][i] = GaussApprox(i, i, f);
    }
    for(int i = 1; i < beta; i++){
        for(int k = i - 1; k >= 0; k--){
            w[i][k] = GaussApprox(i, k, f);
        }
    }
    return w;
}

void Data::precalculateCSSandCS(Fraction* f){
    if(!f->weightsSet){
        f->weightsSet = true;
        int beta = (int) f->seg->size();
        double t = 0;
        double tt = 0;
        f->se1->push_back(t);
        f->sse1->push_back(tt);
        for(int i = 0; i < beta; i++){
            omega* current = extractDataMap((*(f->seg))[i], 1);
                for(omega::iterator iter = current->begin(); iter != current->end(); ++iter){
                    double a_key = iter->first;
                    int a_value = iter->second;
                    t += a_key * a_value;
                    tt += (a_key * a_key) * a_value;
                }
            delete current;
            f->se1->push_back(t);
            f->sse1->push_back(tt);
        }
        if(f->binary){
            t = 0;
            tt = 0;
            f->se0->push_back(t);
            f->sse0->push_back(tt);
            for(int i = 0; i < beta; i++){
                omega* current = extractDataMap((*(f->seg))[i], 0);
                for(omega::iterator iter = current->begin(); iter != current->end(); ++iter){
                    double a_key = iter->first;
                    int a_value = iter->second;
                    t += a_key * a_value;
                    tt += (a_key * a_key) * a_value;
                }
                delete current;
                f->se0->push_back(t);
                f->sse0->push_back(tt);
            }
        }
    }
}
