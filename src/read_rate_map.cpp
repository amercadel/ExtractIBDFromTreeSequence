#include <iostream>
#include <fstream>
#include <sstream>

#include "read_rate_map.hpp"




std::vector<std::string> splitOnTab(std::string& line){
    std::vector<std::string> tokens;
    std::istringstream iss(line);
    std::string token;
    while (getline(iss, token, '\t')) {
        tokens.push_back(token);
    }
    return tokens;
}


std::vector<float> rateMapData::interpolateVector(std::vector<int> &sites, std::vector<int> &bp, std::vector<float> &cm) {
    std::vector<float> interpolated_cm;
    int n = sites.size();
    int m = bp.size();
    int j = 0;
    for (int i = 0; i < n; i++){
        while(j < m && bp[j] < sites[i]){
            j++;
        }
        if (j == 0){
            interpolated_cm.push_back(cm[0]);
        }else if (j == m){
            interpolated_cm.push_back(cm[m - 1]);
        }else{
            int x0 = bp[j - 1];
            int x1 = bp[j];
            float y0 = cm[j - 1];
            float y1 = cm[j];
            int x = sites[i];
            float y = ((y0 * (x1 - x)) + (y1 * (x - x0))) / (x1 - x0);
            interpolated_cm.push_back(y);
        }
    }
    return interpolated_cm;
}


rateMapData readRateMap(char* filename){
    std::ifstream inputFile;
    std::string line;
    std::vector<int> bp_vec;
    std::vector<float> cm_vec;
    inputFile.open(filename);
    std::cout << "opening file\n";
    rateMapData res;
    while (getline(inputFile, line)) {
        int bp;
        float cm;
        std::string tmp;
        std::stringstream input_str(line); 
        std::vector<std::string> data = splitOnTab(line);
        bp = std::stoi(data[1]);
        cm = std::stof(data[3]);
        bp_vec.push_back(bp);
        cm_vec.push_back(cm);

    }
    res.bp_vec = bp_vec;
    res.cm_vec = cm_vec;
    res.last_bp = bp_vec[bp_vec.size() - 1];
    res.last_cm = cm_vec[cm_vec.size() - 1];
    for (int i = 0; i <= res.last_bp + 1; i ++){
        res.all_sites.push_back(i);
    }
    res.interpolated_cm = res.interpolateVector(res.all_sites, res.bp_vec, res.cm_vec);
    return res;

}