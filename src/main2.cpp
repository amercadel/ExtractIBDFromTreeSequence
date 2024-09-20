#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <sstream>
#include <fstream>
#include <chrono>
#include <mutex>
#include <thread>

#include "tskit.h"
#include "read_rate_map.hpp"
#include "extract_ibd_segments.hpp"

std::mutex mtx;

std::vector<std::vector<int>> createMRCATable(const tsk_treeseq_t &ts){
    int n_samples = ts.num_samples;
    std::vector<std::vector<int>> mrca(n_samples, std::vector(n_samples, 0));
    return mrca;
}
std::vector<std::vector<int>> createLastLeftTable(const tsk_treeseq_t &ts){
    int n_samples = ts.num_samples;
    std::vector<std::vector<int>> last_left(n_samples, std::vector(n_samples, 0));
    return last_left;
};

#define check_tsk_error(val)                                                            \
    if (val < 0) {                                                                      \
        fprintf(stderr, "line %d: %s", __LINE__, tsk_strerror(val));                    \
        exit(EXIT_FAILURE);                                                             \
    }


void extract_segments(rateMapData &gen_map, const tsk_treeseq_t &ts, size_t start_index, size_t end_index, std::vector<std::vector<int> > &mrca_last, std::vector<std::vector<int> > &last_left, float min_cutoff, std::mutex &mutex){
    std::vector<float> &map = gen_map.interpolated_cm;
    int ret = 0;
    tsk_tree_t tree;
    // check_tsk_error(ret);
    ret = tsk_tree_init(&tree, &ts, 0);
    // check_tsk_error(ret);
    ret = tsk_tree_first(&tree);
    // check_tsk_error(ret);
    int n_samples = ts.num_samples; // equivalent to id_end
    double sequence_length = tsk_treeseq_get_sequence_length(&ts);
    tsk_id_t mrca;
    std::ofstream output_file;
    char file_name[100];
    snprintf(file_name, sizeof(file_name), "output_%zu_%zu.txt", start_index, end_index);
    output_file.open(file_name, std::ofstream::out);
    if (!output_file.is_open()) {
        std::cerr << "Failed to open output file." << std::endl;
    }
    mutex.lock();
    for(int i = 0; i < n_samples; i++){
        for(int j = 0; j < n_samples; j++){
            ret = tsk_tree_get_mrca(&tree, i, j, &mrca);
            check_tsk_error(ret);
            mrca_last[i][j] = mrca;
            last_left[i][j] = tree.interval.left;
        }
    }
    int tree_cnt = 0;
    int n_trees = tsk_treeseq_get_num_trees(&ts);

    int min_tree_subsample = 5000;
    int left, bp_end, bp_start, last_tree_pos;
    float gen_end, gen_start;
    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&tree)){
        if (tree.interval.left - last_tree_pos < min_tree_subsample){
            tree_cnt++;
            continue;
        }
        last_tree_pos = tree.interval.left;
        for(size_t i = start_index; i < end_index; i++){
            for(size_t j = i + 1; j < n_samples; j++){
                ret = tsk_tree_get_mrca(&tree, i, j, &mrca);
                if(mrca_last[i][j] != mrca){
                    left = tree.interval.left;
                    bp_end = left;
                    bp_start = last_left[i][j];
                    if (bp_end > map.size()){
                        gen_end = map[map.size() - 1];
                    }else{
                        gen_end = map[bp_end];
                    }
                    if (bp_start > map.size()){
                        gen_start = map[map.size()];
                    }else{
                        gen_start = map[bp_start];
                    }
                    if (gen_end - gen_start < min_cutoff){
                        last_left[i][j] = left;
                        mrca_last[i][j] = mrca;
                        continue;
                    }
                    if (gen_end - gen_start > min_cutoff){
                        float len = gen_end - gen_start;
                        ibd_segment seg = ibd_segment(i, j, bp_start, bp_end, gen_start, gen_end, len);
                        std::string out = seg.to_string();
                        output_file << out << std::endl;
                        last_left[i][j] = left;
                        mrca_last[i][j] = mrca;
                    }

                }

            }
        }
        tree_cnt++;


    }
    mutex.unlock();
    for (size_t i = start_index; i < end_index; i++){
        for (size_t j = i + 1; j < n_samples; j++){
            bp_start = last_left[i][j];
            
            bp_end = static_cast<int>(sequence_length);
            gen_end = 0;
            gen_start = 0;
            if (bp_end > map.size()){
                gen_end = map[map.size() - 1];
            }else{
                gen_end = map[bp_end];   
            }
            if(bp_start > map.size()){
                gen_start = map[map.size() - 1];

            }else{
                gen_start = map[bp_start];
            }
            
            if (gen_end - gen_start > min_cutoff){
                float len = gen_end - gen_start;
                ibd_segment seg = ibd_segment(i, j, bp_start, bp_end, gen_start, gen_end, len);
                std::string out = seg.to_string();
                output_file << out << std::endl;
            }
        }
    }
    output_file.close();

    tsk_tree_free(&tree);
}
int main(int argc, char* argv[]){
    char *ts_file = argv[1];
    size_t start_index = std::stoul(argv[2]);
    int end_index = std::stoul(argv[3]);
    char *genetic_map_file = argv[4];
    float minimum_cutoff = std::stof(argv[5]);
    auto start = std::chrono::steady_clock::now();
    std::cout << "reading rate_map" << std::endl;
    rateMapData gen_map = readRateMap(genetic_map_file);
    std::cout << "rate map loaded" << std::endl;
    tsk_treeseq_t ts;
    int ret = 0;
    ret = tsk_treeseq_load(&ts, ts_file, 0);
    check_tsk_error(ret);
    std::vector<std::vector<int>> mrca_last = createMRCATable(ts);
    std::vector<std::vector<int>> last_left = createLastLeftTable(ts);
    std::thread threads[5];
    for(size_t i = 0; i < 5; i++){
        size_t s = (i * 60);
        size_t e = ((i + 1) * 60);
        threads[i] = std::thread(extract_segments, std::ref(gen_map), std::ref(ts), s, e, std::ref(mrca_last), std::ref(last_left), minimum_cutoff, std::ref(mtx));
    }

    for (auto& th : threads) th.join();
    tsk_treeseq_free(&ts);
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << std::chrono::duration<double>(diff).count() << " seconds" << std::endl;
    return 0;

}



