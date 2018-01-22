#ifndef DATA_H
#define DATA_H

#include <map>
#include <vector>
#include <utility>

#include "Fraction.h"

/**
 * This file includes the method to read the data and several scoring methods.
 **/

typedef std::map< double, int > omega;
typedef std::pair<double,double> sc_pair; //sum count pair

struct count_entry {
    int pc0, pc1, ps0, ps1;
    void add(count_entry e){
        pc1 += e.pc1;
        pc0 += e.pc0;
        ps1 += e.ps1;
        ps0 += e.ps0;
    }
    void init(){
        pc1 = 0;
        pc0 = 0;
        ps1 = 0;
        ps0 = 0;
    }
};

struct data_segment {
    omega *data1;
    omega *data0;
    int pc0, pc1, ps0, ps1;
    void add(count_entry e){
        pc1 += e.pc1;
        pc0 += e.pc0;
        ps1 += e.ps1;
        ps0 += e.ps0;
    }
    void set(count_entry e){
        pc1 = e.pc1;
        pc0 = e.pc0;
        ps1 = e.ps1;
        ps0 = e.ps0;
    }
    void add(data_segment e){
        pc1 += e.pc1;
        pc0 += e.pc0;
        ps1 += e.ps1;
        ps0 += e.ps0;
        if(e.data1){
            for(omega::iterator iter = e.data1->begin(); iter != e.data1->end(); ++iter){
                double a_key = iter->first;
                int a_value = iter->second;
                data1->operator[](a_key) += a_value;
            }
        }
        if(e.data0){
            for(omega::iterator iter = e.data0->begin(); iter != e.data0->end(); ++iter){
                double a_key = iter->first;
                int a_value = iter->second;
                data0->operator[](a_key) += a_value;
            }
        }
    }
    void init(){
        pc1 = 0;
        pc0 = 0;
        ps1 = 0;
        ps0 = 0;
    }
    void set0(omega *new_data){
        data0 = new_data;
    }
    void set1(omega *new_data){
        data1 = new_data;
    }
};

typedef std::vector< std::vector< count_entry > > count_table;
typedef std::vector< std::vector< data_segment > > data_map_table;

class Data {
    
public:
    // assign rows to classifier
    std::vector<int> cl0;
    std::vector<int> cl1;
    std::vector<int> all;
    std::vector<double> css0;
    std::vector<double> css1;
    std::vector<double> cs0;
    std::vector<double> cs1;
    int n,m;
    double min, max;
    int kbins;
    double log2sqrt2pi;
    double logbase;
    int resolution;
    bool commulative;
    char method;
    S::costfunptr eval;
    S::matrix data;
    Data() : n(0), m(0) {
        S::matrix dummy(0, std::vector<double>(0));
        data = dummy;
    }
    Data(const char*, bool, bool, char, int, bool);
    
    void read_file(const char* filename, bool, bool, char, int, bool);
    S::weights precomputeWeightsSAEMap(S::segmentation seg);
    S::weights precomputeWeightsGaussFast(Fraction*);
    unsigned int getRowCount();
    unsigned int getColCount();
    void setData(std::vector<std::vector<double> > input, bool binaryClassifierInLastRow, bool comulat, char meth, int s,  bool ignoreLabel);
    double& getElement(unsigned int row,unsigned int col);

protected:
    int stepSize;
    double GaussApprox(int,int,Fraction*);
    double GaussApprox(int,int,int, Fraction*);
    sc_pair getDensities(S::interval, int);
    count_entry getEntry(S::interval);
    double getComplexity(count_entry);
    double ModelWeight(S::interval, data_segment,int,int);
    double Gauss(S::interval, sc_pair, omega*, int, int, int);
    double Laplace(S::interval, sc_pair, omega*, int);
    double SAE(S::interval, sc_pair, omega*);
    omega* extractDataMap(S::interval, int);
    void precalculateCSSandCS(Fraction*);
};

#endif
