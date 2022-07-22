/**
 * @file inferencer.cpp
 * @author Jan Šamánek (jansamanek@email.cz)
 * @brief See the header file.
 * @version 1.1
 * @date 2022-06-30
 * @copyright Copyright (c) 2022
 * This software is open source. There is no license to use this software. Do whatever you want with it.
 * [DISCLAIMER] code may not satisfy clean-code requirements. Scroll down on your own risk.
 */

#include "inferencer.hpp"
#include "utils.hpp"

vector<int64_t>::const_iterator lower_bound_linear(vector<int64_t>::const_iterator bit, vector<int64_t>::const_iterator eit, int64_t wanted)
{
    return (bit < eit and *bit < wanted) ? 
            lower_bound_linear(++bit, eit, wanted) :
            bit;
    // vector<int64_t>::const_iterator it = bit;
    // for (; it < eit; it++)
    // {
    //     if (*it >= wanted)
    //     {
    //         return it;
    //     }
    // }

    // return it;
}

/**
 * @brief Construct a new Inferencer:: Inferencer object
 * 
 * @param g Grid structure.
 */
Inferencer::Inferencer(const Grid<int64_t> &g) : grid(g) 
{}

/**
 * @brief Destroy the Inferencer:: Inferencer object
 * 
 */
Inferencer::~Inferencer() 
{
    this->indices.clear();
}

/**
 * @brief Looking for the sequence of continues numbers in the given vectors. Returns on first miss.
 * 
 * @param wanted Wanted number.
 * @param lvl Index of the vector to search.
 * @return size_t - Returns sequence length
 */
size_t Inferencer::search_subsequence_recursive(const int64_t wanted, const size_t lvl)
{
    if (lvl >= indices.size()) { return 1; }

    for (auto version_indices : indices[lvl])
    {
        // auto it = lower_bound(version_indices->begin(), version_indices->end(), wanted);   
        auto it = (version_indices->size() < 50) ? 
                    lower_bound_linear(version_indices->begin(), version_indices->end(), wanted) :
                    lower_bound(version_indices->begin(), version_indices->end(), wanted);
        if (it != version_indices->end() and (size_t) (*it - wanted) < max_jump_val)
        {
            return 1 + search_subsequence_recursive(wanted + offset, lvl+1);
        }
    }

    return 1;
}

/**
 * @brief Reset for inferencer.
 * 
 */
void Inferencer::flush() 
{ 
    this->indices.clear();
}

/**
 * @brief Function for single-end alignment
 * 
 * @param bit Begin of the read to align
 * @param eit End of the read to align
 * @return Alignment 
 */
Alignment Inferencer::operator ()(string::const_iterator bit, string::const_iterator eit)
{
    Alignment best_alignment(eit-bit);
    if ((size_t) (eit - bit) < 2*kmer_size)
    {
        return best_alignment;
    }
    
    const size_t kmer_n = (eit - bit) / this->offset;
    indices.reserve(kmer_n);

    for (auto it = bit; it <= eit - kmer_size; it += this->offset)
    {        
        auto grid_indices = kcm(it, it + kmer_size);
        this->add(grid_indices);
    }

    // Choose best alignment for read
    for (size_t i = 0; i < kmer_n-1;)
    {
        Alignment alignment(eit - bit - i*offset);
        this->fast_align(i, alignment);
        // printf("    Aligned %lu", alignment.aligned_num_parts);
        // sum_num_parts += alignment.aligned_num_parts;

        if (alignment.aligned_num_parts > best_alignment.aligned_num_parts)
        {
            best_alignment = alignment;        
        }

        if ( best_alignment.aligned_num_parts >= kmer_n - i  or
             best_alignment.aligned_num_parts <= 1 )
        { 
            break; 
        }
    }
    
    this->flush();
    return best_alignment;
}


/**
 * @brief Function for pair-end alignment
 * 
 * @param bit1 Begin of the first read.
 * @param eit1 End of the first read.
 * @param bit2 Begin of the second read.
 * @param eit2 End of the second read.
 * @return Alignment 
 */
Alignment Inferencer::operator ()(string::const_iterator bit1, string::const_iterator eit1, string::const_iterator bit2, string::const_iterator eit2)
{
    if ((size_t) (eit1 - bit1) < 2*kmer_size or (size_t) (eit2 - bit2) < 2*kmer_size)
    {
        return Alignment();
    }

    Alignment result((eit1-bit1) + (eit2-bit2));
    result = this->operator()(bit1, eit1);
    
    // Did not align anywhere
    if (result.start_ids.empty())
    {
        return result;
    }

    this->flush();

    this->search_pairend(result, bit2, eit2);
    this->flush();
    return result;
}

/**
 * @brief Caches the indices from the grid to the class property.
 * 
 * @param grid_indecies All grid indices for current kmer
 */
void Inferencer::add(const vector<size_t> &grid_indecies)
{
    vector<const vector<int64_t>*> version_indices;
    for (size_t grid_index : grid_indecies)
    {
        version_indices.push_back(&grid[grid_index]);
    }
    
    indices.push_back(version_indices);
}

/**
 * @brief Given the alignment of the first half, look for the rest based on the results.
 * 
 * @param first Alignment for the first read
 */
void Inferencer::search_pairend(Alignment& first, string::const_iterator bit, string::const_iterator eit)
{
    size_t kmer_n = (eit - bit) / offset;
    indices.reserve(kmer_n);
    string rc = reverse_complement(bit, eit);
    
    // Push new indices
    for (auto it = rc.begin(); it <= rc.end() - kmer_size; it += this->offset)
    {        
        auto grid_indices = kcm(it, it + kmer_size);
        this->add(grid_indices);
    }

    vector<int64_t> first_ids;
    first_ids.reserve(first.start_ids.size());
    size_t best_seq_len = 0;

    for (const int64_t idx : first.start_ids)
    {
        size_t end_lvl = 0;
        int64_t end_idx = 0;
        // int64_t reversed_idx = (idx < 0) ? idx + genome_length : idx - genome_length;
        const size_t seq_len = this->determine_seq_len(idx, 0, end_lvl, end_idx, true);
        // printf("  idx = %ld, rev_idx = %ld, seq_len = %lu, best_seq_len = %lu\n", idx, reversed_idx, seq_len, best_seq_len);
    
        if (seq_len > best_seq_len)
        {
            first_ids.clear();
            first_ids.push_back(idx);

            best_seq_len = seq_len;

        } else if (seq_len > 0 and seq_len == best_seq_len)
        {
            first_ids.push_back(idx);
        }
    }

    first.aligned_num_parts += best_seq_len;
    first.start_ids = first_ids;
}

/**
 * @brief Determines the sequence length for given number and lvl.
 * 
 * @param wanted Searched number.
 * @param lvl Index of the vector to search
 * @param end_lvl Last lvl where the match was found. (filled via reference)
 * @param end_idx Last index where the match was found. (filled via reference)
 * @return size_t - Return sequence length.
 */
size_t Inferencer::determine_seq_len(int64_t wanted, size_t lvl, size_t& end_lvl, int64_t& end_idx, const bool searching_pairend)
{
    if (lvl >= indices.size()) { return 0; }

    for (auto version_indices : indices[lvl])
    {
        // auto it = lower_bound(version_indices->begin(), version_indices->end(), wanted);
        auto it = (version_indices->size() < 50) ? 
                    lower_bound_linear(version_indices->begin(), version_indices->end(), wanted) :
                    lower_bound(version_indices->begin(), version_indices->end(), wanted);

        if ((it < version_indices->end() and (size_t) (*it - wanted) < max_jump_val) or 
            (searching_pairend and it > version_indices->begin() and (size_t) (wanted - *(it-1)) < max_jump_val))
        {
            end_lvl = lvl;
            end_idx = *it + kmer_size;
            return 1 + this->determine_seq_len(*it + offset, lvl+1, end_lvl, end_idx, false);
        }
    }

    return this->determine_seq_len(wanted + offset, lvl+1, end_lvl, end_idx, searching_pairend);
}

/**
 * @brief Does the alignment of loaded read on to the reference genome.
 * 
 * @param i Last match index. (filled via reference)
 * @param alignment Alignment structure. (filled via reference)
 */
void Inferencer::fast_align(size_t &i, Alignment &alignment)
{
    vector<int64_t> v;
    size_t best_seq_len;

    for (; i < Inferencer::indices.size()-1;)
    {
        best_seq_len = 1;

        for (const auto version_indices : indices[i])
        {
            for (const int64_t idx : *version_indices)
            {
                const size_t seq_len = this->search_subsequence_recursive(idx + offset, i + 1);

                // Found best alignment so far
                if (seq_len > best_seq_len)
                {
                    v.clear();
                    v.push_back(idx);
                    best_seq_len = seq_len;
                
                // Found alignmet of the same quality
                } else if (seq_len == best_seq_len and seq_len > 1)
                {
                    v.push_back(idx);
                }
            }
        }

        alignment.aligned_num_parts = best_seq_len;
        alignment.aligned_block_end_idx = i + best_seq_len;

        // Precise alignment
        if (alignment.aligned_block_end_idx >= indices.size())
        {
            alignment.start_ids = v;
            i = alignment.aligned_block_end_idx;
            return;
        
        // Found some alignment but not precise
        } else if (not v.empty())
        {
            break;
        
        // Did not found any match
        } else
        {
            i++;
        }
    }

    // No match was found
    if (v.empty())
    {
        i++;
        return;
    }

    // Go over all found alignments and determine their full lengths 
    alignment.start_ids.reserve(v.size());
    i = alignment.aligned_block_end_idx;
    size_t i_copy = i;
    size_t full_best_seq_len = best_seq_len;
    for (const int64_t idx : v)
    {
        size_t end_lvl = i_copy;
        int64_t end_idx = 0;
        size_t full_seq_len = best_seq_len + this->determine_seq_len(idx + (i_copy+1)*offset, i_copy+1, end_lvl, end_idx, false);

        if (full_seq_len > full_best_seq_len)
        {
            alignment.start_ids.clear();
            alignment.start_ids.push_back(idx);
            alignment.aligned_num_parts = full_seq_len;
            alignment.aligned_block_end_idx = end_lvl;
            alignment.read_size = end_idx - idx;
            full_best_seq_len = full_seq_len;
            i = end_lvl;
        
        } else if (full_seq_len == full_best_seq_len)
        {
            alignment.start_ids.push_back(idx);
        }
    }
}