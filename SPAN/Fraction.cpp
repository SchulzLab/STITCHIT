#include "Fraction.h"

#include <cmath>
#include <assert.h>

Fraction::Fraction(int s, int e, int sS, bool b){
    assert(s < e);
    start = s;
    end = e;
    binary = b;
    stepSize = sS;
    sse1 = new std::vector<double>();
    se1 = new std::vector<double>();
    sse0 = new std::vector<double>();
    se0 = new std::vector<double>();
    seg = new S::segmentation();
    // true after first iteration -- don't calculate costs or initial costs again
    weightsSet = false;
    merged = false;
}

Fraction::~Fraction(){
    delete sse1;
    delete sse0;
    delete se1;
    delete se0;
    delete seg;
}

void Fraction::mergeIn(Fraction* f){
    merged = true;
    if(f->start < start){
        assert(f->end == start);
        start = f->start;
        // push front
        for(int i = 0; i < f->seg->size(); i++){
            seg->insert(seg->begin(), (*f->seg)[i]);
            // shift because of zero in beginning
            sse1->insert(sse1->begin() + 1, (*(f->sse1))[(i+1)]);
            se1->insert(se1->begin() + 1, (*(f->se1))[(i+1)]);
            if(binary){
                sse0->insert(sse0->begin() + 1, (*(f->sse0))[(i+1)]);
                se0->insert(se0->begin() + 1, (*(f->se0))[(i+1)]);
            }
        }
        // add increased sum to old ones
        double sse1_max = (*(f->sse1))[f->sse1->size() - 1];
        double se1_max = (*(f->se1))[f->se1->size() - 1];
        double sse0_max = 0.0;
        double se0_max = 0.0;
        if(binary){
            sse0_max = (*(f->sse0))[f->sse0->size() - 1];
            se0_max = (*(f->se0))[f->se0->size() - 1];
        }
        for(int i = f->sse1->size(); i < sse1->size(); i++){
            (*sse1)[i] += sse1_max;
            (*se1)[i] += se1_max;
            if(binary){
                (*sse0)[i] += sse0_max;
                (*se0)[i] += se0_max;
            }
        }
    }else{
        assert(end == f->start);
        end = f->end;
        double sse1_max = (*sse1)[sse1->size() - 1];
        double se1_max = (*se1)[se1->size() - 1];
        double sse0_max = 0.0;
        double se0_max = 0.0;
        if(binary){
            sse0_max = (*sse0)[sse0->size() - 1];
            se0_max = (*se0)[se0->size() - 1];
        }
        for(int i = 0; i < f->seg->size(); i++){
            seg->push_back((*(f->seg))[i]);
            
            sse1->push_back((*(f->sse1))[i+1] + sse1_max);
            se1->push_back((*(f->se1))[i+1] + se1_max);
            if(binary){
                sse0->push_back((*(f->sse0))[i+1] + sse0_max);
                se0->push_back((*(f->se0))[i+1] + se0_max);
            }
        }
    }
}

void Fraction::updateSegmentation(S::segmentation s){
    S::segmentation old = *seg;
    while(seg->size() > 0){
        seg->pop_back();
    }
    int old_index = 0;
    std::vector<int> toDelete;
    for(int i = 0; i < s.size(); i++){
        while(old[old_index].end != s[i].end){
            toDelete.push_back(old_index);
            old_index++;
        }
        old_index++;
        seg->push_back(s[i]);
    }
    adjustErrors(toDelete);
}

void Fraction::adjustErrors(std::vector<int> toDelete){
    if((sse1->size() - 1) != seg->size()){
        cutItOff(sse1, toDelete);
        cutItOff(se1, toDelete);
        if(binary){
            cutItOff(sse0, toDelete);
            cutItOff(se0, toDelete);
        }
    }
}

void Fraction::cutItOff(std::vector<double>* dum, std::vector<int> toDelete){
    int i = 0;
    int del_i = 0;
    for(std::vector<double>::iterator it = dum->begin() + 1; it != dum->end() - 1;){
        if(del_i >= toDelete.size())
            break;
        if(i == toDelete[del_i]){
            it = dum->erase(it);
            del_i++;
        }else{
            ++it;
        }
        i++;
    }
    assert((dum->size() - 1) == seg->size());
}

int Fraction::size(){
    return end - start;
}
