#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <sstream>
#include <fstream>
#include <thread>
#include <algorithm>
#include <functional>
#include <mutex>
#include <chrono>

#include "tskit.h"
#include "read_rate_map.hpp"
#include "extract_ibd_segments.hpp"



#define check_tsk_error(val)                                                            \
    if (val < 0) {                                                                      \
        fprintf(stderr, "line %d: %s", __LINE__, tsk_strerror(val));                    \
        exit(EXIT_FAILURE);                                                             \
    }


int extract_segments(rateMapData &gen_map, tsk_treeseq_t &ts, float min_cutoff){
    std::vector<float> map = gen_map.interpolated_cm;
    int ret = 0;
    int min_tree_subsample = 5000;
    tsk_tree_t tree;
    check_tsk_error(ret);
    ret = tsk_tree_init(&tree, &ts, 0);
    check_tsk_error(ret);
    ret = tsk_tree_first(&tree);
    check_tsk_error(ret);
    int n_samples = ts.num_samples;
    double sequence_length = tsk_treeseq_get_sequence_length(&ts);
    int start_index = 0;
    int end_index = n_samples;
    std::vector<std::vector<int>> mrca_last = createMRCATable(ts);
    std::vector<std::vector<int>> last_left = createLastLeftTable(ts);
    tsk_id_t mrca;
    std::ofstream output_file;
    // Define a char array with sufficient size
    char file_name[100];

    // Use sprintf to format the string
    sprintf(file_name, "output_%d_%d.txt", start_index, end_index);
    output_file.open(file_name, std::ofstream::out);
    if (!output_file.is_open()) {
        std::cerr << "Failed to open output file." << std::endl;
        return 1;
    }

    
    for (int i = start_index; i < end_index; i++){
        for(int j = i + 1; j < n_samples; j++){
            ret = tsk_tree_get_mrca(&tree, i, j, &mrca);
            check_tsk_error(ret);
            mrca_last[i][j] = mrca;
            last_left[i][j] = tree.interval.left;
            
        }
    }
    int left;
    int genomic_end;
    int genomic_start;
    float gen_end;
    float gen_start;
    int last_tree_pos = 0;
    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&tree)){
        if (tree.interval.left - last_tree_pos < min_tree_subsample){
            continue;
        }
        last_tree_pos = tree.interval.left;
        for (size_t i = start_index; i < end_index; i++){
            for(size_t j = i + 1; j < n_samples; j++){
                ret = tsk_tree_get_mrca(&tree, i, j, &mrca);
                if (mrca_last[i][j] != mrca){
                    left = tree.interval.left;
                    genomic_end = left;
                    
                    genomic_start = last_left[i][j];
                    gen_end = 0;
                    gen_start = 0;
                    if (genomic_end > map.size()){
                        gen_end = map[map.size() - 1];
                    }else{
                        gen_end = map[genomic_end];
                    }
                    
                    if (genomic_start > map.size()){
                        gen_start = map[map.size()];
                    }else{
                        gen_start = map[genomic_start];
                    }
                    if (gen_end - gen_start < min_cutoff){
                        last_left[i][j] = left;
                        mrca_last[i][j] = mrca;
                        continue;
                    }
                    if (gen_end - gen_start > min_cutoff){
                        float len = gen_end - gen_start;
                        ibd_segment seg = ibd_segment(i, j, genomic_start, genomic_end, gen_start, gen_end, len);
                        std::string out = seg.to_string();
                        output_file << out << std::endl;
                    }
                    last_left[i][j] = left;
                    mrca_last[i][j] = mrca;
                }

            }
        }
    }
    
    for (size_t i = start_index; i < end_index; i++){
        for (size_t j = i + 1; j < n_samples; j++){
            genomic_start = last_left[i][j];
            
            genomic_end = static_cast<int>(sequence_length);
            gen_end = 0;
            gen_start = 0;
            if (genomic_end > map.size()){
                gen_end = map[map.size() - 1];
            }else{
                gen_end = map[genomic_end];   
            }
            if(genomic_start > map.size()){
                gen_start = map[map.size() - 1];

            }else{
                gen_start = map[genomic_start];
            }
            
            if (gen_end - gen_start > min_cutoff){
                float len = gen_end - gen_start;
                ibd_segment seg = ibd_segment(i, j, genomic_start, genomic_end, gen_start, gen_end, len);
                std::string out = seg.to_string();
                output_file << out << std::endl;
            }
                

            }
        }   
    output_file.close();
    return 0;
}

std::pair<int, int> *generate_subsets(int n_cpus, int n_haplotypes){
    std::pair<int, int>* arr =  new std::pair<int, int>[n_cpus];
    
    if (n_haplotypes % n_cpus == 0){
        int haps_per_cpu = n_haplotypes / (n_cpus);
        for (int i = 0; i < n_cpus; i++){
            std::pair<int, int> p(i * haps_per_cpu, (i + 1) * haps_per_cpu);
            arr[i] = p;
        }
    }
    else{
        int haps_per_cpu = n_haplotypes / (n_cpus - 1);
        for (int i = 0; i < n_cpus - 1; i++){
            std::pair<int, int> p(i * haps_per_cpu, (i + 1) * haps_per_cpu);
            arr[i] = p;
        }
        std::pair<int, int> p((n_cpus - 1) * haps_per_cpu, std::min((n_cpus - 1) * haps_per_cpu, n_haplotypes));
        arr[n_cpus - 1] = p;
        }

    return arr;

}


int main(int argc, char* argv[]){
    char *ts_file = argv[1];
    char *genetic_map_file = argv[2];
    float minimum_cutoff = std::stof(argv[3]);
    auto start = std::chrono::steady_clock::now();
    tsk_treeseq_t ts;
    int ret = 0;
    ret = tsk_treeseq_load(&ts, ts_file, 0);
    check_tsk_error(ret);
    int n_samples = ts.num_samples;
    std::cout << "reading rate_map" << std::endl;
    rateMapData gen_map = readRateMap(genetic_map_file);
    std::cout << "rate map loaded" << std::endl;
    auto end = std::chrono::steady_clock::now();
    int ret_val = extract_segments(gen_map, ts, minimum_cutoff);
    // Store the time difference between start and end
    auto diff = end - start;
    std::cout << std::chrono::duration<double>(diff).count() << " seconds" << std::endl;

    
    return ret_val;





    
    



}


