#pragma once
#include <vector>
#include <iostream>
#include <string>
#include "tskit.h"


struct ibd_segment{
    int id1;
    int id2;
    int bp_start;
    int bp_end;
    float cm_start;
    float cm_end;
    float len_cm;
    void display(){
        std::cout << id1 << "\t" << id2 << "\t" << bp_start << "\t" << bp_end << "\t" << cm_start << "\t" << \
        cm_end << "\t" << len_cm << std::endl;
        }
    
    std::string to_string() {
        return std::to_string(id1) + "\t" + std::to_string(id2) + "\t" + std::to_string(bp_start) + "\t" +
                std::to_string(bp_end) + "\t" + std::to_string(cm_start) + "\t" + std::to_string(cm_end) + "\t" +
                std::to_string(len_cm);
    }
    
    ibd_segment(int id1, int id2, int bp_start, int bp_end, float cm_start, float cm_end, float len_cm) :
        id1(id1), id2(id2), bp_start(bp_start), bp_end(bp_end), cm_start(cm_start), cm_end(cm_end), len_cm(len_cm) {}

};

std::vector<std::vector<int>> createMRCATable(const tsk_treeseq_t &ts);
std::vector<std::vector<int>> createLastLeftTable(const tsk_treeseq_t &ts);
