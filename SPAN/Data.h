#ifndef DATA_H
#define DATA_H

#include <map>
#include <vector>
#include <utility>
#include <set>

#include "Fraction.h"

/**
 * This file includes the method to read the data and several scoring methods.
 **/

typedef std::map< double, int > omega;
typedef std::map<int, std::vector<int> > rowmap;

class Data {
    
public:
    // assign rows to classifier
    rowmap crows;
    std::set<int> categories;
    std::vector<int> all;
    int n,m;
    double min, max;
    int kbins;
    double log2sqrt2pi;
    double logbase;
    int resolution;
    char method;
    S::matrix data;
    Data() : n(0), m(0) {
        S::matrix dummy(0, std::vector<double>(0));
        data = dummy;
    }
    Data(const char*, bool, char, int, bool);
    void read_file(const char* filename, bool, char, int, bool);
    S::weights precomputeWeightsGaussFast(Fraction&) const;
    unsigned int getRowCount();
    unsigned int getColCount();
    void setData(std::vector<std::vector<double> > input, bool, char, int,  bool);
    double& getElement(unsigned int row,unsigned int col);
    std::vector<double>& getRow(unsigned int row);
    int categoryCount() const; 
    std::set<int> getCategories(); 

protected:
    int stepSize;
    double GaussApprox(int,int,Fraction&) const;
    double GaussApprox(int,int,int, Fraction&) const;
    omega extractDataMap(S::interval, int) const;
    void precalculateCSSandCS(Fraction&) const;
};

#endif
