/**
 * @file utils.hpp
 * @author Jan Šamánek (jansamanek@email.cz)
 * @brief Component of the Grizzly software.
 * Containes utility functions.
 * @version 1.1
 * @date 2022-06-30
 * @copyright Copyright (c) 2022
 * This software is open source. There is no license to use this software.
 * [DISCLAIMER] code may not satisfy clean-code requirements. Scroll down on your own risk.
 */

#pragma once

#include"config.hpp"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <omp.h>
#include <iterator>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include <boost/program_options.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>

struct Alignment;
struct MergedGenes_s;
struct GenesAnnotation_s;

/**
 * @brief Function to clean last line in the command line
 * 
 */
void clean_console_line();

/**
 * @brief Function to clean the previous line in the command line
 * 
 */
void clean_previous_console_line();

/**
 * @brief Utility function for printing the progress bar.
 * 
 * @param progress Current value.
 * @param size Max value.
 * @param msg Additional message.
 * @param end End char.
 */
void print_progress(const size_t progress, const size_t size, const string msg = "", const string end = "\r");

/**
 * @brief Get the file size from file descriptor.
 * 
 * @param file File descriptor.
 * @return size_t - Number of bytes
 */
size_t get_file_size(istream &file);

/**
 * @brief Transforms sequence to reverse_complement.
 * 
 * @param bit Start of the sequence.
 * @param eit End of the sequence.
 * @return string - Reverse complement
 */
string reverse_complement(string::const_iterator bit, string::const_iterator eit);

/**
 * @brief File opener using mapped file.
 * 
 * @param filename Input file.
 * @return boost::iostreams::mapped_file. 
 */
boost::iostreams::mapped_file open_file(const char* filename);

/**
 * @brief Function to parse the input arguments.
 * 
 * @param argc Number of args.
 * @param argv Arguments.
 */
void parse_args(int argc, const char* argv[]);

/**
 * @brief Check whether the given sequnce belongs to given gene id.
 * 
 * @param genes Gene annotation.
 * @param gene_idx Gene id.
 * @param start_idx Start of the sequence.
 * @param end_idx End of the sequence.
 * @return true - if belongs to the given gene.
 * @return false - otherwise
 */
bool check_gene_belonging(const vector<Gene_s> &genes, size_t gene_idx, int64_t start_idx, int64_t end_idx);

/**
 * @brief Checks if given sequence is mainly soft-masked.
 * 
 * @param bit Start of the sequence.
 * @param eit End of the sequence.
 * @return true - If is mainly masked.
 * @return false - Otherwise.
 */
bool is_mainly_masked(string::const_iterator bit, string::const_iterator eit);

/**
 * @brief Merge genes from gene annotation.
 * 
 * @param ga Gene annotation.
 * @return MergedGenes_s structure
 */
MergedGenes_s merge_genes(const GenesAnnotation_s &ga);

/**
 * @brief Get the gene ids for given alignment
 * 
 * @param mg Merged genes.
 * @param alignment Alignment.
 * @return unordered_set<int64_t> - Gene ids. 
 */
unordered_set<int64_t> get_gene_ids(const MergedGenes_s &mg, const Alignment &alignment);

/**
 * @brief Parallel reduction of matrix. N rows to the first one.
 * 
 * @param m Matrix. (Number of rows = number of threads)
 */
void paralel_reduction(Vec2D<float> &m);

/**
 * @brief Method for distribution multimapped fragments. (Not in use)
 * 
 * @param hard Counter of unique alignments
 * @param soft Counter of multimapped alignments
 * @param total_distribution_sum Total sum
 * @param hard_sum Sum of hard
 * @param soft_sum Sum of soft
 */
void statistical_magic(vector<float> &hard, vector<float> &soft, double& total_distribution_sum, float hard_sum, float soft_sum);

/**
 * @brief Utility function for printing vector
 * 
 * @tparam T typename in the vector
 * @param v vector
 */
template<typename T>
void print_vector(const vector<T> &v)
{
    if (v.empty())
    {
        printf("[]\n");
        return;
    }
    
    printf("[");
    for (size_t i = 0; i < v.size()-1; i++)
    {
        cout << v[i] << ", ";
    }
    cout << v.back() << "]" << endl;
}

/**
 * @brief argsort
 * 
 * @tparam T type in given vector.
 * @param v Vector.
 * @return vector<size_t> - Sorted indices. 
 */
template<typename T>
vector<size_t> argsort(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

/**
 * @brief Shuffle given vector based on given indices
 * 
 * @tparam T vector type.
 * @param data Data vector.
 * @param ids Sorted indices.
 */
template<typename T>
void sort_by_indices(vector<T>& data, const vector<size_t> ids)
{
    vector<T> temp(data.size());

    for(size_t i = 0 ; i < data.size() ; i++){
        temp[i] = data[ids[i]];
    }

    data = temp;
}

/**
 * @brief Takes N bytes from 'i' index and recast it to given type. Moves 'i'.
 * 
 * @tparam T type to recast
 * @param bytes file in char* form
 * @param i index
 * @return const T recasted value 
 */
template<typename T>
const T recast(const char* bytes, size_t &i)
{
    T num;
    copy(&(bytes[i]), (&bytes[i]) + sizeof(num), reinterpret_cast<char*>(&num));
    i += sizeof(num);
    return num;
}