/**
 * @file inferencer.hpp
 * @author Jan Šamánek (jansamanek@email.cz)
 * @brief Component of the Grizzly software.
 * Contains class for aligning fragments to reference genome in grid structure.
 * @version 1.1
 * @date 2022-06-30
 * 
 * @copyright Copyright (c) 2022
 * This software is open source. There is no license to use this software.
 * [DISCLAIMER] code may not satisfy clean-code requirements. Scroll down on your own risk.
 */

#pragma once
#include "config.hpp"
#include "kmer_converter.hpp"


/**
 * @brief Inference class for mapping read on to the reference genome
 * 
 */
class Inferencer
{
   private:
    vector<vector<const vector<int64_t>*>> indices;
    size_t offset = kmer_size;
    KmerConverterManager kcm;
    const Grid<int64_t> &grid;

    size_t search_subsequence_recursive(const int64_t wanted, const size_t lvl);
    size_t determine_seq_len(int64_t wanted, size_t lvl, size_t& end_lvl, int64_t& end_idx, const bool searching_pairend);
    void add(const vector<size_t> &grid_indecies);
    void fast_align(size_t &i, Alignment &alignment);
    void search_pairend(Alignment& first, string::const_iterator bit, string::const_iterator eit);
    void flush();

   public:

    Inferencer(const Grid<int64_t> &g);
    ~Inferencer();

    /**
     * @brief Function for single-end alignment
     * 
     * @param bit Begin of the read to align
     * @param eit End of the read to align
     * @return Alignment 
     */
    Alignment operator ()(string::const_iterator bit, string::const_iterator eit);

    /**
     * @brief Function for pair-end alignment
     * 
     * @param bit1 Begin of the first read.
     * @param eit1 End of the first read.
     * @param bit2 Begin of the second read.
     * @param eit2 End of the second read.
     * @return Alignment 
     */
    Alignment operator ()(string::const_iterator bit1, string::const_iterator eit1, string::const_iterator bit2, string::const_iterator eit2);
};
