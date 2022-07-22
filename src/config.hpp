/**
 * @file config.hpp
 * @author Jan Šamánek (jansamanek@email.cz)
 * @brief Component of the Grizzly software.
 * Contains global variables, structures, macros and typedefs.
 * @version 1.1
 * @date 2022-06-30
 * 
 * @copyright Copyright (c) 2022
 * This software is open source. There is no license to use this software. Do whatever you want with it.
 * [DISCLAIMER] code may not satisfy clean-code requirements. Scroll down on your own risk.
 */

#pragma once

#include <vector>
#include <string>
#include <chrono>
#include <boost/program_options.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>

using namespace std;
using namespace std::chrono;
namespace po = boost::program_options;

// GLOBAL VARIABLES
extern size_t kmer_size;
extern size_t max_n_num;
extern size_t max_jump_val;
extern size_t genome_length;
extern size_t num_threads;
extern bool include_softmask_bases;
extern bool perform_tests;
extern bool do_pair_end_alignment;
extern float print_progress_frequence;
extern float minimum_align;
extern string genome_file;
extern string gene_file;
extern vector<string> alignment_files;
extern string database_file;
extern string load_grid_file;
extern string save_grid_file;

extern const string file_prefix;
extern const int64_t anchor;



#define ever (;1;)

// Tests character on newline char
#define is_linefeedchar(x) (x == (char) 10)
// Test thread on being last
#define is_last_thread(x) (x == (num_threads-1))
// Test thread on being first
#define is_first_thread(x) (x == 0)


template<typename T>
using Vec2D = vector<vector<T>>;

template<typename T>
using Vec2D_p = vector<vector<T>*>;

template<typename T>
using Vec3D_p = vector<vector<vector<T>*>>;

template<typename T>
using Cell = vector<T>;

template<typename T>
using Grid = vector<Cell<T>>;

/**
 * @brief Simple structure to store the information about genes
 * 
 */
struct Gene_s {
    int64_t start_idx;
    int64_t end_idx;
    string gene;
    string name;
    bool is_empty = true;
    bool is_rcg;

    Gene_s(int64_t s_idx, int64_t e_idx, string g_name, bool rcg)
    {
        start_idx = s_idx;
        end_idx = e_idx;
        name = g_name;
        is_rcg = rcg;
        is_empty = true;
    }

    Gene_s(){}

};

/**
 * @brief Wrapper around Gene_s structures
 * 
 */
struct GenesAnnotation_s
{
    vector<Gene_s> rcgs;
    vector<Gene_s> fwgs;
};

/**
 * @brief Transformed GenesAnnotation_s structure with merged genes
 * 
 */
struct MergedGenes_s
{
    vector<int64_t> start_ids;
    vector<int64_t> end_ids;
    vector<string> names;
    vector<bool> is_rcg;
};

/**
 * @brief Structure for the actual alignment result
 * 
 */
struct Alignment
{
    vector<int64_t> start_ids;
    size_t aligned_num_parts = 1;
    size_t aligned_block_end_idx = 0;
    size_t read_size = 0;

    Alignment(){}

    Alignment(const size_t rs)
    { 
        read_size = rs;
        start_ids.reserve(rs / kmer_size); 
    }
};
