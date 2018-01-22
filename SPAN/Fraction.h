#ifndef FRACTION_H
#define FRACTION_H

#include <vector>

#include "defs.h"

class Fraction {
    
public:
    int start, end, stepSize;
    bool binary, weightsSet, merged;
    double initial, compressed;
    S::segmentation* seg;
    std::vector<double>* sse0;
    std::vector<double>* sse1;
    std::vector<double>* se0;
    std::vector<double>* se1;

    // functinons
    /**
     * start must be the first included element and end is excluded [start,end)
     **/
    Fraction(int,int,int,bool);
    ~Fraction();
    
    void mergeIn(Fraction*);
    void updateSegmentation(S::segmentation);
    int size();
    
protected:
    void cutItOff(std::vector<double>*, std::vector<int>);
    // pass indices which should be deleted
    void adjustErrors(std::vector<int>);
};

#endif
