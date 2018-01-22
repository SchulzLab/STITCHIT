#ifndef SEGMENT_H
#define SEGMENT_H

#include <vector>
#include <iostream>

#include "defs.h"
#include "Data.h"

/**
 * This class contains the main execution methods for the Span algorithm ('Still' project name).
 **/

struct LessThanByImputity{
    bool operator()(const S::qe& lhs, const S::qe& rhs) const{
        return lhs.gain < rhs.gain;
    }
};

class Binning {
    
protected:
    int length;
    Data data;
    bool verbose;
    
    void printBlockLists(std::vector<S::block>, std::vector<S::block>);
    void applyDPFlexi(int, Fraction*);
    
    std::vector<int> Merge(std::vector<int>,std::vector<int>);
    double modelCost(int,int);
    
public:
    // constructors
    Binning(Data&,bool);
    void runSPAN(int, Fraction*);
    void printSegmentation(S::segmentation);
};

#endif
