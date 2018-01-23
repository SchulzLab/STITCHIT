#ifndef SPAN_H
#define SPAN_H

#include <getopt.h> //options
#include <stdlib.h> //set handler
#include <cstdio>
#include <fstream>
#include <cmath>
#include <string>
#include <omp.h> //Added by fschmidt
#include "Wall_Time.h"
#include "Data.h"
#include "Segment.h"
#include "Fraction.h"
#include "defs.h"

class SPAN {
    
public:
	SPAN(){};
	std::vector<std::pair<unsigned int, unsigned int> > runSpan(Data& d, int s, unsigned int maxCores);

};

#endif
