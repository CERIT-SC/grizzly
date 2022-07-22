/**
 * @file config.cpp
 * @author Jan Šamánek (jansamanek@email.cz)
 * @brief See header file
 * @version 1.1
 * @date 2022-06-30
 * 
 * @copyright Copyright (c) 2022
 * No copyright.
 */

#include "config.hpp"
#include <string>
#include <vector>
#include <limits>

size_t kmer_size = 14;
size_t max_n_num = 2;
size_t max_jump_val = 1000;
size_t genome_length = 0;
size_t num_threads = 1;
bool include_softmask_bases = false;
bool perform_tests = false;
bool do_pair_end_alignment = false;
float print_progress_frequence = 0.01;
float minimum_align = 0.0;
string genome_file;
string gene_file;
vector<string> alignment_files;
string database_file;
string load_grid_file;
string save_grid_file;

const string file_prefix("Samanek&Juric-GridMapping");
const int64_t anchor = numeric_limits<int64_t>::max();
