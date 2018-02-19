#ifndef DEFS_H
#define DEFS_H

#include <vector>
#include <iostream>
#include <cmath>

namespace S {
    typedef std::vector< std::vector<double> > matrix;
    typedef double (*costfunptr)(double,double);
    
    // define return interval, whereas [start, end)
    struct interval {
        int start, end;
        bool inRange(int p) { return (p >= start && p < end); }
        void print(){
            std::cout << "[" << start << ", " << end << ")\n";
        }
        void joinWith(interval iv){
            if(iv.start < start){
                start = iv.start;
            }
            if(iv.end > end){
                end = iv.end;
            }
        }
        int length(){
            return end - start;
        }
    };
    
    // define bin with complexity
    struct bin {
        interval b;
        double complexity;
        int size;
        void init(interval bl, int s, double e){ b = bl; size = s; complexity = e; }
    };
    
    typedef std::vector< interval > segmentation;
    
    struct binning {
        segmentation bins;
        double totalCosts;
    };
    
    // define struct for priority queue
    struct qe {
        bin b;
        double gain;
    };
    
    typedef std::vector< std::vector < double > > weights;
    
    double inline log2N(double z){
        if(z < 1){
            std::cout << "ERROR: log2N not defined for z < 1! --> returned 0\n";
            return 0;
        }else{
            double logstar = log2(z);
            double sum = logstar;
            while(true){
                logstar = log2(logstar);
                if(logstar <= 0){
                    break;
                }else{
                    sum += logstar;
                }
            }
            return sum + log2(2.865064);
        }
    }
    
    double inline log2fac(int n){
        double sum = 0;
        for(int i = 2; i <= n; i++){
            sum += log2(i);
        }
        return sum;
    }
    
    double inline log2nChoosek(int n, int k){
        if(k > n || k == 0){
            return 0;
        }else{
            return log2fac(n) - log2fac(k) - log2fac(n-k);
        }
    }
}

#endif
