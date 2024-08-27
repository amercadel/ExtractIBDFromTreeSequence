#include <vector>
#include <iostream>
#include "tskit.h"
#include "extract_ibd_segments.hpp"


std::vector<std::vector<int>> createMRCATable(tsk_treeseq_t &ts){
    int n_samples = ts.num_samples;
    std::vector<std::vector<int>> mrca(n_samples, std::vector(n_samples, 0));
    return mrca;
}
std::vector<std::vector<int>> createLastLeftTable(tsk_treeseq_t &ts){
    int n_samples = ts.num_samples;
    std::vector<std::vector<int>> last_left(n_samples, std::vector(n_samples, 0));
    return last_left;
};




