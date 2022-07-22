/**
 * @file kmer_converter.cpp
 * @author Jan Šamánek (jansamanek@email.cz)
 * @brief See the header file.
 * @version 1.1
 * @date 2022-06-30
 * @copyright Copyright (c) 2022
 * This software is open source. There is no license to use this software. Do whatever you want with it.
 * [DISCLAIMER] code may not satisfy clean-code requirements. Scroll down on your own risk.
 */

#include "config.hpp"
#include "utils.hpp"
#include "kmer_converter.hpp"

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>

/**
 * @brief Construct a new Kmer Converter:: Kmer Converter object
 * 
 */
KmerConverter::KmerConverter() {}
    
/**
 * @brief Construct a new Kmer Converter:: Kmer Converter object
 * 
 * @param kc Kmer Converter to copy
 */
KmerConverter::KmerConverter(const KmerConverter& kc)
{
    value = kc.value;
}

/**
 * @brief Add function for adding new nucleotid to bitset
 * 
 * @param nucl new nucleotid {A,C,G,T}
 */
void KmerConverter::add(char nucl)
{
    // Creates bitset by masking nucl number and taking second and third bit
        // A = 65 = 0b01000|00|1 => 0
        // C = 67 = 0b01000|01|1 => 1
        // T = 84 = 0b01010|10|0 => 2
        // G = 71 = 0b01000|11|1 => 3
    value <<= 2;
    value |= (nucl & 0b110) >> 1;
}

/**
 * @brief Converts value to grid index
 * 
 * @return index pointing into a grid
 */
size_t KmerConverter::get_index()
{
    return value;
}


/**
 * @brief Resets the value of the kmer
 * 
 */
void KmerConverter::reset()
{
    value = 0;
}

// KMER MANAGER FUNCTIONS
// --------------------------

/**
 * @brief Construct a new Kmer Converter Manager:: Kmer Converter Manager object
 * 
 */
KmerConverterManager::KmerConverterManager()
{
    converters.push_back(KmerConverter());
}

/**
 * @brief Function for converting kmer to index
 * 
 * @param kmer kmer
 * @return const vector<size_t> - Return all indices for given kmer
 */
const vector<size_t> KmerConverterManager::operator() (const string &kmer)
{
    for (char nucl : kmer)
    {
        this->add(nucl);
    }

    auto grid_indices = this->get_indices();
    this->flush();
    return grid_indices;
}

/**
 * @brief Function for converting sequence to the indices
 * 
 * @param bit Begin of the sequence
 * @param eit End of the sequence
 * @return const vector<size_t> - Return all indices for given kmer
 */
const vector<size_t> KmerConverterManager::operator() (string::const_iterator bit, string::const_iterator eit)
{
    for (auto it = bit; it < eit; it++)
    {
        this->add(*it);
    }

    auto grid_indices = this->get_indices();
    this->flush();
    return grid_indices;
}

void KmerConverterManager::add(char nucl)
{
    nucl = toupper(nucl);
    if (nucl == 'N')
    {
        const size_t size = converters.size();
        for (size_t i = 0; i < size; i++)
        {
            for(char n : {'C', 'G', 'T'})
            {
                KmerConverter kc(converters[i]);
                kc.add(n);
                converters.push_back(kc);
            }
            converters[i].add('A');
        }

    // Nucleotide is not N
    } else
    {
        for (auto& converter : converters)
        {
            converter.add(nucl);
        }
    }
}

const vector<size_t> KmerConverterManager::get_indices()
{
    vector<size_t> indices;
    indices.reserve(converters.size());

    for (auto c : converters)
    {
        indices.push_back(c.get_index());
    }

    return indices;
}

void KmerConverterManager::flush()
{
    converters.resize(1);
    converters[0].reset();
}