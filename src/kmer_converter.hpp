/**
 * @file kmer_converter.hpp
 * @author Jan Šamánek (jansamanek@email.cz)
 * @brief Component of the Grizzly software.
 * Containes 2 classes for converting kmer to number.
 * @version 1.1
 * @date 2022-06-30
 * 
 * @copyright Copyright (c) 2022
 * This software is open source. There is no license to use this software.
 * [DISCLAIMER] code may not satisfy clean-code requirements. Scroll down on your own risk.
 */

#pragma once
#include "config.hpp"

/**
 * @brief Class for converting kmer to index in the grid
 * This class is managed by class KmerConverterManager and should not be used alone.
 */
class KmerConverter
{
   public:
 
    size_t value = 0;

    /**
     * @brief Construct a new Kmer Converter object
     * 
     */
    KmerConverter();

    /**
     * @brief Construct a new Kmer Converter object
     * 
     * @param kc Kmer converter to copy
     */
    KmerConverter(const KmerConverter& kc);

    /**
     * @brief Add function for adding new nucleotid to bitset
     * 
     * @param nucl new nucleotid {A,C,G,T}
     */
    void add(char nucl);

    /**
     * @brief Converts value to grid index
     * 
     * @return index pointing into a grid
     */
    size_t get_index();

    /**
     * @brief Resets the value of the kmer
     * 
     */
    void reset();
};

/**
 * @brief Manager for converting kmer into a index
 */
class KmerConverterManager
{
   private:
    vector<KmerConverter> converters;

    /**
     * @brief Add function for adding new nucleotid {A,C,G,T,N}
     * Manager does not check MAX_N_NUM assertion. 
     * It's programmers responsibility if necessery.
     * 
     * @param nucl new nucleotid {A,C,G,T,N}
     */
    void add(char nucl);

    /**
     * @brief Converts all versions of loaded kmer into indices
     * 
     * @return Grid indices of all kmer versions
     */
    const vector<size_t> get_indices();

    /**
     * @brief Reset function.
     * 
     */
    void flush();

   public:
    KmerConverterManager();

    /**
     * @brief Function for converting kmer to index
     * 
     * @param kmer kmer
     * @return const vector<size_t> - Return all indices for given kmer
     */
    const vector<size_t> operator() (const string &kmer);

    /**
     * @brief Function for converting sequence to the indices
     * 
     * @param bit Begin of the sequence
     * @param eit End of the sequence
     * @return const vector<size_t> - Return all indices for given kmer
     */
    const vector<size_t> operator() (string::const_iterator bit, string::const_iterator eit);
};