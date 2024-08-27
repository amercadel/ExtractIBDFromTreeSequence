#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <cmath>


std::vector<std::string> splitOnTab(std::string& line);
struct rateMapData{
    std::vector<int> bp_vec;
    std::vector<float> cm_vec;
    std::vector<int> all_sites;
    std::vector<float> interpolated_cm;
    std::vector<float> interpolateVector(std::vector<int> &sites, std::vector<int> &bp_vec, std::vector<float> &cm_vec);
    int last_bp;
    float last_cm;

    void display(int idx){
        std::cout << "bp value: " << bp_vec[idx] << std::endl << "cm value: " << cm_vec[idx] << std::endl;
    };

};
rateMapData readRateMap(char *filename);
