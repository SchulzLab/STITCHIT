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
    
    void applyDPFlexi(int, Fraction&, bool) const;
    
    std::vector<int> Merge(std::vector<int>&,std::vector<int>&) const;
    double modelCost(int,int) const;
    
public:
    // constructors
    Binning(Data&,bool);
    void runSPAN(int, Fraction&, bool) const;
    void printSegmentation(S::segmentation) const;
};

#endif
