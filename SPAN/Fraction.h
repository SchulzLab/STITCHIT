#ifndef FRACTION_H
#define FRACTION_H

#include <vector>
#include "defs.h"
#include <map>

#include <set>

typedef std::map<int, std::vector<double> > datamap;
typedef std::map<int,double> simplemap;

class Fraction {
    
public:
    int start, end, stepSize;
    bool weightsSet, merged;
    std::set<int> categories;
    double initial, compressed, initialData;
    int initialBins;
    S::segmentation seg;
    datamap sse;
    datamap se;
    
    /**
     * start must be the first included element and end is excluded [start,end)
     **/
    Fraction() = default;
    Fraction(int,int,int,std::set<int>);
    
    void mergeIn(Fraction&);
    void updateSegmentation(S::segmentation&);
    int size();
    
protected:
    void cutItOff(std::vector<double>&, std::vector<int>&);
    void adjustErrors(std::vector<int>&);
};

#endif
