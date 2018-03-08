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

int Data::categoryCount() const{ 
	return categories.size(); 
}

std::set<int> Data::getCategories(){ 
	return categories;
}

void Data::setData(std::vector<std::vector<double> > input,  bool classifierInLastRow, char meth, int s,  bool ignoreLabel){
	method = meth;
	stepSize = s;
	log2sqrt2pi = log2(std::sqrt(2 * 3.14159265358979323846));
	logbase = log(2);
	resolution = 10;
	data=input;
	categories.clear();

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

	for (int i = 0; i < m; i++){
	          if(classifierInLastRow || ignoreLabel){
			int classifier = data[i][n-1];
			categories.insert(classifier);
			if(!ignoreLabel){
				 crows[classifier].push_back(i);
			}else{
				categories.insert(1);
				crows[1].push_back(i);
				}
			data[i].pop_back();
		}else{
			categories.insert(1);
			crows[1].push_back(i);
			} 
		all.push_back(i);
		}
	}



//This function assumes that no sample IDs are provided
void Data::read_file(const char* filename, bool classifierInLastRow, char meth, int s,  bool ignoreLabel){
    method = meth;
    stepSize = s;
    log2sqrt2pi = log2(std::sqrt(2 * 3.14159265358979323846));
    logbase = log(2);
    
    std::ifstream f;
    std::string line;
    categories.clear();
    f.open(filename);
    if (f.is_open()) {                                                                       
        unsigned int counter = 0;
        while (std::getline(f, line)) {
            std::vector<double> new_vec;
            std::istringstream is(line);
            double actual;
            while (is >> actual) {
                new_vec.push_back(actual);
            }
            if(classifierInLastRow || ignoreLabel){
                int classifier = new_vec[new_vec.size() - 1];
                categories.insert(classifier);
                if(!ignoreLabel){
                    crows[classifier].push_back(counter);
                }else{
                    categories.insert(1);
                    crows[1].push_back(counter);
                }
                new_vec.pop_back();
            }else{
                categories.insert(1);
                crows[1].push_back(counter);
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


Data::Data(const char* filename, bool binaryClassifierInLastRow, char m, int s, bool ignoreLabel){
    read_file(filename, binaryClassifierInLastRow, m, s, ignoreLabel);
}

omega Data::extractDataMap(S::interval iv, int label) const{
    std::vector<int> const &rows = crows.at(label);
    omega elementCount;
    for(int i = iv.start; i < iv.end; i++){
        for(unsigned long k = 0; k < rows.size(); k++){
            int j = rows[k];
            elementCount[(data[j][i])] += 1;
        }
    }
    return elementCount;
}

double Data::GaussApprox(int i, int j, Fraction& f) const{
    double score = 0.0;
    std::set<int>::iterator it;
    for(it = categories.begin(); it != categories.end(); ++it){
        int label = *it;
        double ci = GaussApprox(i, j, label, f) / (double) crows.at(label).size();
        score += ci;
    }
    return score;
}

double Data::GaussApprox(int i, int j, int label, Fraction& f) const{
    double sigma = 0.0;
    int ii = (j + 1);
    int jj = (i + 1);
    int diff = (f.seg)[i].end - (f.seg)[j].start;
    std::vector<double> css;
    std::vector<double> cs;
    int N = crows.at(label).size();
    css = f.sse[label];
    cs = f.se[label];
    double myweight = (double) diff * (double) N;
    double sse = (css[jj] - css[ii-1]) - (1/myweight) * ((cs[jj] - cs[ii-1]) * (cs[jj] - cs[ii-1]));
    if(sse <= 0.0)
        return 0.0;
    sigma = std::sqrt(sse / (double) myweight);
    double score = (myweight/2) * log2(sigma*sigma*2.0*M_PI) + ((myweight/2.0)*log(2.0)) + N * log2(resolution);
    score = score < 0.0 ? 0.0 : score;
    return score;
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

S::weights Data::precomputeWeightsGaussFast(Fraction& f) const{
    precalculateCSSandCS(f);
    int beta = (int) f.seg.size();
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

void Data::precalculateCSSandCS(Fraction& f) const{
    if(!f.weightsSet){
        f.weightsSet = true;
        int beta = (int) f.seg.size();
        std::set<int>::iterator it;
        for(it = categories.begin(); it != categories.end(); ++it){
            int label = *it;
            double t = 0;
            double tt = 0;
            f.se[label].push_back(t);
            f.sse[label].push_back(tt);
            for(int i = 0; i < beta; i++){
                omega current = extractDataMap((f.seg)[i], label);
                    for(omega::iterator iter = current.begin(); iter != current.end(); ++iter){
                        double a_key = iter->first;
                        int a_value = iter->second;
                        t += a_key * a_value;
                        tt += (a_key * a_key) * a_value;
                    }
                f.se[label].push_back(t);
                f.sse[label].push_back(tt);
            }
        }
    }
}  
