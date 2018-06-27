#include "Fraction.h"

#include <cmath>
#include <assert.h>

Fraction::Fraction(int s, int e, int sS, std::set<int> cat){
    assert(s < e);
    start = s;
    end = e;
    stepSize = sS;
    categories = cat;
    std::set<int>::iterator it;
    for(it = categories.begin(); it != categories.end(); ++it){
        int label = *it;
        sse[label];
        se[label];
    }
    // true after first iteration -- don't calculate costs or initial costs again
    weightsSet = false;
    merged = false;
    initialData = 0.0;
    initialBins = 0;
}

void Fraction::mergeIn(Fraction& f){
    merged = true;
    if(f.start < start){
        assert(f.end == start);
        start = f.start;
        // push front
        for(unsigned int i = 0; i < f.seg.size(); i++){
            seg.insert(seg.begin(), (f.seg)[i]);
            // shift because of zero in beginning
            std::set<int>::iterator it;
            for (it = categories.begin(); it != categories.end(); ++it){
                int label = *it;
                sse[label].insert(sse[label].begin() + 1, (f.sse[label])[(i+1)]);
                se[label].insert(se[label].begin() + 1, (f.se[label])[(i+1)]);
            }
        }
        // add increased sum to old ones
        simplemap maxsse;
        simplemap maxse;
        std::set<int>::iterator it;
        for(it = categories.begin(); it != categories.end(); ++it){
            int label = *it;
            maxsse[label] = (f.sse[label])[f.sse[label].size() - 1];
            maxse[label] = (f.se[label])[f.se[label].size() - 1];
        }
        for(unsigned int i = f.sse[1].size(); i < sse[1].size(); i++){
            std::set<int>::iterator it;
            for(it = categories.begin(); it != categories.end(); ++it){
                int label = *it;
                (sse[label])[i] += maxsse[label];
                (se[label])[i] += maxse[label];
            }
        }
    }else{
        assert(end == f.start);
        end = f.end;
        simplemap maxsse;
        simplemap maxse;
        std::set<int>::iterator it;
        for(it = categories.begin(); it != categories.end(); ++it){
            int label = *it;
            maxsse[label] = (sse[label])[sse[label].size() - 1];
            maxse[label] = (se[label])[se[label].size() - 1];
        }
        for(unsigned int i = 0; i < f.seg.size(); i++){
            seg.push_back((f.seg)[i]);
            std::set<int>::iterator it;
            for(it = categories.begin(); it != categories.end(); ++it){
                int label = *it;
                sse[label].push_back((f.sse[label])[i+1] + maxsse[label]);
                se[label].push_back((f.se[label])[i+1] + maxse[label]);
            }
        }
    }
}

void Fraction::updateSegmentation(S::segmentation& s){
    S::segmentation old = seg;
    while(seg.size() > 0){
        seg.pop_back();
    }
    int old_index = 0;
    std::vector<int> toDelete;
    for(unsigned int i = 0; i < s.size(); i++){
        while(old[old_index].end != s[i].end){
            toDelete.push_back(old_index);
            old_index++;
        }
        old_index++;
        seg.push_back(s[i]);
    }
    adjustErrors(toDelete);
}

void Fraction::adjustErrors(std::vector<int>& toDelete){
    if((sse[1].size() - 1) != seg.size()){
        std::set<int>::iterator it;
        for(it = categories.begin(); it != categories.end(); ++it){
            int label = *it;
            cutItOff(sse[label], toDelete);
            cutItOff(se[label], toDelete);
        }
    }
}

void Fraction::cutItOff(std::vector<double>& dum, std::vector<int>& toDelete){
    int i = 0;
    unsigned int del_i = 0;
    for(std::vector<double>::iterator it = dum.begin() + 1; it != dum.end() - 1;){
        if(del_i >= toDelete.size())
            break;
        if(i == toDelete[del_i]){
            it = dum.erase(it);
            del_i++;
        }else{
            ++it;
        }
        i++;
    }
    assert((dum.size() - 1) == seg.size());
}

int Fraction::size(){
    return end - start;
}
