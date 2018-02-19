#ifndef FRACTION_H
#define FRACTION_H

#include <vector>
#include "defs.h"

//Version 2.0
#include <map>
#include <set>
typedef std::map<int, std::vector<double>* > datamap;
typedef std::map<int,double> simplemap;
//E Version 2.0

class Fraction {
    
public:
    int start, end, stepSize;
    bool weightsSet, merged;
    std::set<int> categories;
    double initial, compressed;
    S::segmentation* seg;
    datamap sse;
    datamap se;
    
    /**
     * start must be the first included element and end is excluded [start,end)
     **/
    Fraction(int,int,int,std::set<int>);
    ~Fraction();
    
    void mergeIn(Fraction*);
    void updateSegmentation(S::segmentation);
    int size();
    
protected:
    void cutItOff(std::vector<double>*, std::vector<int>);
    void adjustErrors(std::vector<int>);
};

#endif
