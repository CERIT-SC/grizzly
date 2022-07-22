/**
 * @file extend_grid.cpp
 * @author Jan Šamánek (jansamanek@email.cz)
 * @brief See header file
 * @version 1.1
 * @date 2022-06-30
 * 
 * @copyright Copyright (c) 2022
 * No copyright.
 */

#include "config.hpp"
#include "kmer_converter.hpp"
#include "utils.hpp"
#include "extend_grid.hpp"

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>
#include <omp.h>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>

/**
 * @brief Utility function to parse one line in input file
 * 
 * @param line Loaded line
 * @param transcript_offsets Annotation of offsets for genes/transcripts 
 * @param transcript_offset Offset of current transcript. Filled via reference.
 * @param relative_index Index in transcript of current read. Filled via reference.
 * @param old_seq Old sequence.
 * @param new_seqs New (modified) sequences.
 */
void parse_line( const string &line, 
                const unordered_map<string, size_t>& transcript_offsets, 
                size_t &transcript_offset, 
                size_t &relative_index, 
                string &old_seq, 
                vector<string> &new_seqs)
{
    transcript_offset = 0;
    relative_index = 0;
    old_seq = "";
    new_seqs.clear();

    uint32_t token_id = 0;
    istringstream iss(line);
    string seq;
    string token;
    char del;
    unordered_map<string, size_t>::const_iterator it;

    while(getline(iss, token, '\t'))
    {
    istringstream isst(token);
        switch (token_id)
        {
            // First parsing transcript id
            case 0:
                it = transcript_offsets.find(token);
                if ( it == transcript_offsets.end() )
                {
                    fprintf(stderr, "Could not found chromozome name %s in saved names. Loaded chromozomes = %lu\nEnding with err. code 2.\n", token.c_str(), transcript_offsets.size());
                    exit(2);
                }

                transcript_offset = it->second;
                break;

            // Then parsing relative index in transcript
            case 1:
                relative_index = stoi(token)-1;
                break;

            // Then parsing old sequence
            case 2:
                old_seq = token;
                break;

            // Last parsing new_seqs
            case 3:
                if (find(token.begin(), token.end(), ';') != token.end())   { del = ';'; } 
                else                                                        { del = ','; }

                while (getline(isst, seq, del)) { new_seqs.push_back(seq); }
                break;

            default:
                break;
        }

        token_id++;
    }
}


/**
 * @brief Get the all kmers that supposed to be put in to the grid
 * Creates every combination of reference genome and new sequences.
 * 
 * @param sequence Reference fragment to generated kmers from.
 * @param new_seqs All the new sequences
 * @param is_rcg Flag to create reverse complement from reads
 * @return vector<pair<size_t, int64_t>> - pairs of all combinations where numbers in pairs are:
 * first - grid index,
 * second - kmer index
 */
vector<pair<size_t, int64_t>> get_all_kmers( const string& sequence, 
                                             const vector<string>& new_seqs,
                                             const bool is_rcg)
{
    KmerConverterManager kcm;
    vector<pair<size_t, int64_t>> kmers;

    for (string new_seq : new_seqs)
    {
        for (int64_t i = 0; i < (int64_t) (sequence.size() - kmer_size); i++)
        {
            string kmer;
            for (size_t j = i; j < (size_t) i+kmer_size; j++)
            {
                kmer.push_back((j < new_seq.size()) ? 
                                new_seq[j] : 
                                sequence[j]);
            }

            if (is_rcg)
            {
                kmer = reverse_complement(kmer.begin(), kmer.end());
            }

            auto grid_indices = kcm(kmer);
            for (uint32_t grid_index : grid_indices)
            {
                kmers.push_back(make_pair(grid_index, i));
            }
        }        
    }

    return kmers;
}

void extend_grid(Grid<int64_t> &grid, 
                 const GenesAnnotation_s& ga,
                 ifstream &database_fd)
{
    size_t bytes_count = 0;

    size_t total_counter = 0;
    size_t align_properly_counter = 0;
    size_t no_gene_counter = 0;
    size_t align_badly_counter = 0;


    const size_t bytes_n = get_file_size(database_fd);
    const float print_offset = bytes_n * print_progress_frequence;
    float next_print_barrier = 0;

    unordered_map<string, size_t> tr_ids;
    for (size_t i=0; i < ga.fwgs.size(); i++)
    {
        tr_ids[ga.fwgs[i].name] = i;
    }


    auto total_time_start = std::chrono::high_resolution_clock::now();  
    string line;
    while ( database_fd.good() )
    {
        // All these variables are filled using references in the function
        size_t tr_idx = 0;
        size_t relative_index = 0;
        string old_seq;
        vector<string> new_seqs;

        getline(database_fd, line);
        parse_line(line, tr_ids, tr_idx, relative_index, old_seq, new_seqs);

        bytes_count += line.size();
        size_t index = ga.fwgs[tr_idx].start_idx + relative_index;
        string sequence;

        size_t new_seq_max_len = max_element(new_seqs.begin(), new_seqs.end(),
                                    [](const string& a, const string& b) {
                                        return a.length() < b.length();
                                    })->size();

        const bool belongs = check_gene_belonging(ga.fwgs, tr_idx, index, index + new_seq_max_len);
        
        if (not belongs)
        {
            no_gene_counter++;
            total_counter++;
            if (bytes_count >= next_print_barrier)
            {
                next_print_barrier += print_offset;
                print_progress(bytes_count, bytes_n);
            }
            continue;
        }

        auto start_it = ga.fwgs[tr_idx].gene.begin() + relative_index;
        auto end_it = start_it + old_seq.length();

        string loaded_old_seq(start_it, end_it);
        if (old_seq != loaded_old_seq)
        {
            align_badly_counter++;
            total_counter++;
            if (bytes_count >= next_print_barrier)
            {
                next_print_barrier += print_offset;
                print_progress(bytes_count, bytes_n);
            }
            continue;
        }
        
        start_it = ga.fwgs[tr_idx].gene.begin() + relative_index + old_seq.length();
        end_it = (start_it + new_seq_max_len + kmer_size < ga.fwgs[tr_idx].gene.end()) ? 
                            start_it + new_seq_max_len + kmer_size : 
                            ga.fwgs[tr_idx].gene.end();

        string tr_fragment_fw(start_it, end_it);

        start_it = ga.rcgs[tr_idx].gene.begin() + relative_index;
        end_it = (start_it + new_seq_max_len + kmer_size < ga.rcgs[tr_idx].gene.end()) ? 
                            start_it + new_seq_max_len + kmer_size : 
                            ga.rcgs[tr_idx].gene.end();

        string tr_fragment_rc(start_it, end_it);

        auto kmer_pairs = get_all_kmers(tr_fragment_fw, new_seqs, false);
        for (auto kmer_pair : kmer_pairs)
        {
            size_t grid_idx = kmer_pair.first;
            size_t kmer_idx = kmer_pair.second;

            auto it = upper_bound(grid[grid_idx].begin(), grid[grid_idx].end(), kmer_idx);
            grid[grid_idx].insert(it, kmer_idx);
        }

        kmer_pairs = get_all_kmers(tr_fragment_rc, new_seqs, true);
        for (auto kmer_pair : kmer_pairs)
        {
            size_t grid_idx = kmer_pair.first;
            size_t kmer_idx = kmer_pair.second;

            auto it = upper_bound(grid[grid_idx].begin(), grid[grid_idx].end(), kmer_idx);
            grid[grid_idx].insert(it, kmer_idx);
        }

        align_properly_counter++;
        total_counter++;
        if (bytes_count >= next_print_barrier)
        {
            next_print_barrier += print_offset;
            print_progress(bytes_count, bytes_n);
        }
    }
    auto total_time_stop = std::chrono::high_resolution_clock::now();  
    std::chrono::duration<double, std::milli> fp_ms = total_time_stop - total_time_start;

    printf("Extending grid summary:\n");
    printf("Total modifications             = %lu\n", total_counter);
    printf("Proper modifications            = %lu (%.2f)\n", align_properly_counter, ((float) align_properly_counter) / total_counter * 100);
    printf("Bad modifications (skipped)     = %lu (%.2f)\n", align_badly_counter, ((float) align_properly_counter) / total_counter * 100);
    printf("No gene modifications (skipped) = %lu (%.2f)\n", no_gene_counter, ((float) align_properly_counter) / total_counter * 100);
    printf("Total time                      = %.2lfs\n", fp_ms.count() / 1000);
}