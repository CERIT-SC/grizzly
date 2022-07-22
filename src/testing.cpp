/**
 * @file testing.cpp
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
#include "inferencer.hpp"
#include "testing.hpp"
#include "grid_mapping.hpp"

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>
#include <omp.h>

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>

/**
 * @brief Tests grid performance on taken reads from the reference genome single-end
 * 
 * @param grid Grid structure
 * @param ga Gene annotation
 * @param mg Merged gene annotation
 * @param repeat_n Number of cycles
 * @param max_n Max N chars. per kmer
 * @param deform Flag to artificialy deform a read
 * @param min_read_size Minimum read length
 * @param max_read_size Maximum read length
 */
void test_inference_speed_on_genome(const Grid<int64_t>& grid,
                                    GenesAnnotation_s &ga,
                                    MergedGenes_s &mg, 
                                    const uint32_t repeat_n, 
                                    const uint32_t max_n, 
                                    const bool deform, 
                                    const size_t min_read_size,
                                    const size_t max_read_size)
{
    double total_cycle_time = 0.0;
    size_t total_cycle_count = 0;
    size_t total_number_of_valid = 0;
    size_t total_number_of_nonaligned = 0;
    size_t total_number_of_overlaps = 0;
    size_t read_len_sum = 0;
    vector<Gene_s> test_genes;

    // CACHING GENES
    // Load genes to string vector so they do not overlap
    printf("\nCaching genes\n");

    ifstream fd(gene_file);
    ga.fwgs.clear();
    ga.rcgs.clear();
    load_transcriptome(fd, ga);
    mg = merge_genes(ga);
    fd.close();
    
    test_genes.insert(test_genes.end(), ga.rcgs.begin(), ga.rcgs.end());
    test_genes.insert(test_genes.end(), ga.fwgs.begin(), ga.fwgs.end());
    
    clean_console_line();
    clean_previous_console_line();


    float print_offset = repeat_n * print_progress_frequence;
    float next_print_barrier = 0;
    
    printf("Testing performance\n");
    auto total_time_start = std::chrono::high_resolution_clock::now();  
    #pragma omp parallel default(shared) num_threads(num_threads)
    {
        double cycle_time = 0.0;
        Inferencer inferencer(grid);
        uint32_t valid_results = 0;
        uint32_t number_of_overlaps = 0;
        uint32_t number_of_nonaligned = 0;
        const size_t thread_idx = omp_get_thread_num();
        size_t thread_cycle_num = repeat_n / num_threads;
        if (is_last_thread(thread_idx))
        {
            thread_cycle_num += (repeat_n % num_threads);
        }

        for (size_t i = 0; i < thread_cycle_num; i++)
        {
            // printf("i=%lu\n", i);
            uint32_t gene_idx = 0;
            uint32_t n_count = 0;
            uint32_t read_start = 0;
            uint32_t read_end = 0;
            size_t read_size = 0;
            string read;
            Gene_s tg;

            // GENERATING READ
            while (true)
            {
                gene_idx = rand() % test_genes.size();
                read_size = (min_read_size + (rand() % (max_read_size - min_read_size)));
                tg = test_genes[gene_idx];
                const string &gene = tg.gene;

                if (gene.length() < read_size) { read.clear(); continue; }

                read_start = rand() % ((gene.length() - read_size) + 1);
                read_end = read_start + read_size;

                copy(gene.begin() + read_start, gene.begin() + read_end, back_inserter(read));

                bool is_good = true;
                // Check if any of the kmers in the read has over maximum number of allowed Ns
                for (string::iterator it = read.begin(); it < read.end(); it += kmer_size)
                {
                    n_count = count(it, it+kmer_size, 'N');
                    if (n_count > max_n or (not include_softmask_bases and is_mainly_masked(it, it+kmer_size))) 
                    { is_good = false; break; }
                }
                
                if (not is_good) { read.clear(); continue; }
                else { break; }
            }
            // printf("  Generated\n");

            #pragma omp atomic
            read_len_sum += read.size();

            // DEFORMATION
            uint32_t deform_n = 0;
            if (deform && rand() > RAND_MAX/2)
            {
                // printf("  Deforming\n");
                string selection("ACGT");
                // Max 3 deforms
                deform_n = 1+rand() % 3;
                // string orig_read = read;
                for (uint32_t j = 0; j < deform_n; j++)
                {
                    size_t index = rand() % read_size;
                    auto it = find(selection.begin(), selection.end(), read[index]);
                    if (it == selection.end())
                    { read[index] = selection[rand() % selection.size()]; } 
                    else
                    { read[index] = selection[(it - read.begin() + 1) % selection.size()]; }
                    
                    // printf("Deform on index = %lu from %c to %c\n", index, read[index], selection[(index+1) % selection.size()]);
                }
                // cout << orig_read << endl;
                // cout << read << endl;
            }

            // printf("  Deformed\n");

            #pragma omp atomic
            total_cycle_count++;

            // ALIGNMENT
            auto cycle_start = std::chrono::high_resolution_clock::now();  

            // printf("Start idx = %ld\n", (tg.is_rcg) ? tg.start_idx + read_start - genome_length : tg.start_idx + read_start);
            // printf("  Start\n");
            auto alignment = inferencer(read.begin(), read.end());
            // printf("  Stop. Start ids size = %lu\n", alignment.start_ids.size());
            
            auto cycle_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> fp_ms = cycle_end - cycle_start;            
            cycle_time += fp_ms.count();
            
            if (alignment.start_ids.empty())
            {
                // printf("Did not align\n\n");
                number_of_nonaligned++;
                continue;
            }
            // printf("\n");
            
            // printf("  Getting genes\n");
            auto pred_gene_ids = get_gene_ids(mg, alignment);
            // printf("  Done\n");
            
            // bool found = false;
            for (int64_t gene_idx : pred_gene_ids)
            {
                if (mg.names[gene_idx] == tg.name)
                {
                    valid_results++;
                    // found = true;
                    break;
                }
            }

            // if (not found)
            // {
            //     printf("Did not align well.\n");
            // }

            // printf("\n");            
            
            #pragma omp master
            if (total_cycle_count >= next_print_barrier)
            {
                next_print_barrier += print_offset;
                print_progress(total_cycle_count, repeat_n);
            }
        }

        #pragma omp atomic
        total_number_of_valid += valid_results;
        #pragma omp atomic
        total_cycle_time += cycle_time;
        #pragma omp atomic
        total_number_of_overlaps += number_of_overlaps;
        #pragma omp atomic
        total_number_of_nonaligned += number_of_nonaligned;
    }

    auto total_time_stop = std::chrono::high_resolution_clock::now();  
    std::chrono::duration<double, std::milli> fp_ms = total_time_stop - total_time_start;
    double total_duration = fp_ms.count();

    const size_t total_number_of_invalid = total_cycle_count - total_number_of_valid - total_number_of_overlaps - total_number_of_nonaligned;
    const double valid_estimate = total_number_of_valid / (double) total_cycle_count * 100;
    const double invalid_estimate = total_number_of_invalid / (double) total_cycle_count * 100;
    const double overlap_estimate = total_number_of_overlaps / (double) total_cycle_count * 100;
    const double nonaligned_estimate = total_number_of_nonaligned / (double) total_cycle_count * 100;

    clean_previous_console_line();
    printf("Testing Done!\n");
    printf("Total number of aligments        = %lu\n", total_cycle_count);
    printf("Valid aligments                  = %lu (%.2lf %%)\n", total_number_of_valid, valid_estimate);
    printf("Invalid alignments               = %lu (%.2lf %%)\n", total_number_of_invalid, invalid_estimate);
    printf("Overlapped alignments            = %lu (%.2lf %%)\n", total_number_of_overlaps, overlap_estimate);
    printf("Non-aligned alignments           = %lu (%.2lf %%)\n", total_number_of_nonaligned, nonaligned_estimate);
    printf("Average read length              = %.1f chars.\n\n", ((double) read_len_sum) / total_cycle_count);

    printf("Time summary:\n");
    printf("Number of alignments  = %lu\n", total_cycle_count);
    printf("Total cycle time      = %.3lfs (%.3lfms/thread)\n", total_cycle_time / 1000.0, total_cycle_time / num_threads);
    printf("Average cycle time    = %.4lfms/read/thread\n", total_cycle_time / total_cycle_count / num_threads);
    printf("Total time            = %.3lfs (Includes time for generating reads and other)\n\n", total_duration/1000);
}

/**
 * @brief Tests grid performance on taken reads from the reference genome pair-end
 * 
 * @param grid Grid structure
 * @param ga Gene annotation
 * @param mg Merged gene annotation
 * @param repeat_n Number of cycles
 * @param max_n Max N chars. per kmer
 * @param deform Flag to artificialy deform a read
 * @param min_read_size Minimum read length
 * @param max_read_size Maximum read length
 */
void test_inference_speed_on_genome_pairend(const Grid<int64_t>& grid,
                                            GenesAnnotation_s &ga,
                                            MergedGenes_s &mg, 
                                            const uint32_t repeat_n, 
                                            const uint32_t max_n, 
                                            const bool deform, 
                                            const size_t min_read_size,
                                            const size_t max_read_size)
{
    double total_cycle_time = 0.0;
    size_t total_cycle_count = 0;
    size_t total_number_of_valid = 0;
    size_t total_number_of_nonaligned = 0;
    size_t total_number_of_overlaps = 0;
    size_t read_len_sum = 0;
    vector<Gene_s> test_genes;

    // CACHING GENES
    // Load genes to string vector so they do not overlap
    printf("\nCaching genes\n");

    ifstream fd(gene_file);
    ga.fwgs.clear();
    ga.rcgs.clear();
    load_transcriptome(fd, ga);
    mg = merge_genes(ga);
    fd.close();
    
    test_genes.insert(test_genes.end(), ga.rcgs.begin(), ga.rcgs.end());
    test_genes.insert(test_genes.end(), ga.fwgs.begin(), ga.fwgs.end());
    
    clean_console_line();
    clean_previous_console_line();

    float print_offset = repeat_n * print_progress_frequence;
    float next_print_barrier = 0;
    
    printf("Testing performance\n");
    auto total_time_start = std::chrono::high_resolution_clock::now();  
    #pragma omp parallel default(shared) num_threads(num_threads)
    {
        double cycle_time = 0.0;
        Inferencer inferencer(grid);
        uint32_t valid_results = 0;
        uint32_t number_of_overlaps = 0;
        uint32_t number_of_nonaligned = 0;
        const size_t thread_idx = omp_get_thread_num();
        size_t thread_cycle_num = repeat_n / num_threads;
        if (is_last_thread(thread_idx))
        {
            thread_cycle_num += (repeat_n % num_threads);
        }

        for (size_t i = 0; i < thread_cycle_num; i++)
        {
            // printf("i=%lu\n", i);
            uint32_t gene_idx = 0;
            uint32_t n_count = 0;
            uint32_t read1_start = 0;
            uint32_t read1_end = 0;
            uint32_t read2_start = 0;
            uint32_t read2_end = 0;
            size_t min_seq_size = 0;
            size_t read1_size = 0;
            size_t read2_size = 0;
            string read1, read2;
            Gene_s tg;

            // GENERATING READ
            while (true)
            {
                gene_idx = rand() % test_genes.size();
                read1_size = min_read_size + (rand() % (max_read_size - min_read_size));
                read2_size = min_read_size + (rand() % (max_read_size - min_read_size));
                min_seq_size = read1_size + read2_size + max_jump_val;

                tg = test_genes[gene_idx];
                const string &gene = tg.gene;

                if (gene.length() < min_seq_size) { continue; }

                read1_start = rand() % (gene.length() - min_seq_size + 1);
                read1_end = read1_start + read1_size;

                read2_start = read1_start + rand() % (max_jump_val);
                read2_end = read2_start + read2_size;

                copy(gene.begin() + read1_start, gene.begin() + read1_end, back_inserter(read1));
                copy(gene.begin() + read2_start, gene.begin() + read2_end, back_inserter(read2));

                bool is_good = true;
                // Check if any of the kmers in the read has over maximum number of allowed Ns
                for (string::iterator it = read1.begin(); it < read1.end(); it += kmer_size)
                {
                    n_count = count(it, it+kmer_size, 'N');
                    if (n_count > max_n or (not include_softmask_bases and is_mainly_masked(it, it+kmer_size))) 
                    { is_good = false; break; }
                }

                for (string::iterator it = read2.begin(); it < read2.end(); it += kmer_size)
                {
                    n_count = count(it, it+kmer_size, 'N');
                    if (n_count > max_n or (not include_softmask_bases and is_mainly_masked(it, it+kmer_size))) 
                    { is_good = false; break; }
                }
                
                if (not is_good) { read1.clear(); read2.clear(); continue; }
                else { break; }
            }
            // printf("  Generated\n");

            #pragma omp atomic
            read_len_sum += read1.length() + read2.length();

            // DEFORMATION
            uint32_t deform_n = 0;
            if (deform && rand() > RAND_MAX/2)
            {
                // printf("  Deforming\n");
                string selection("ACGT");
                // Max 3 deforms
                deform_n = 1+rand() % 2;
                // string orig_read = read;
                for (uint32_t j = 0; j < deform_n; j++)
                {
                    size_t index = rand() % read1_size;
                    auto it = find(selection.begin(), selection.end(), read1[index]);
                    if (it == selection.end())
                    { read1[index] = selection[rand() % selection.size()]; } 
                    else
                    { read1[index] = selection[(it - read1.begin() + 1) % selection.size()]; }
                    
                    // printf("Deform on index = %lu from %c to %c\n", index, read[index], selection[(index+1) % selection.size()]);
                }

                for (uint32_t j = 0; j < deform_n; j++)
                {
                    size_t index = rand() % read2_size;
                    auto it = find(selection.begin(), selection.end(), read2[index]);
                    if (it == selection.end())
                    { read2[index] = selection[rand() % selection.size()]; } 
                    else
                    { read2[index] = selection[(it - read2.begin() + 1) % selection.size()]; }
                    
                    // printf("Deform on index = %lu from %c to %c\n", index, read[index], selection[(index+1) % selection.size()]);
                }
                // cout << orig_read << endl;
                // cout << read << endl;
            }

            read2 = reverse_complement(read2.begin(), read2.end());

            // printf("Gene idx = %u, gene_start = %ld, gene_end = %ld is_rcg = %s\n", gene_idx, tg.start_idx, tg.end_idx, (tg.is_rcg) ? "true" : "false");
            // if (tg.is_rcg)
            // {
            //     printf("  read1 read_start = %ld (%ld), read_end = %ld (%ld)\n", tg.start_idx + read1_start, tg.start_idx + read1_start - genome_length, tg.start_idx + read1_end, tg.start_idx + read1_end - genome_length);
            //     printf("  read2 read_start = %ld, read_end = %ld\n", tg.start_idx + read2_start, tg.start_idx + read2_end);
                
            // } else
            // {
            //     printf("  read1 read_start = %ld, read_end = %ld\n", tg.start_idx + read1_start, tg.start_idx + read1_end);
            //     printf("  read2 read_start = %ld (%ld), read_end = %ld (%ld)\n", tg.start_idx + read2_start, tg.start_idx + read2_start - genome_length, tg.start_idx + read2_end, tg.start_idx + read2_end - genome_length);
            // }
            
            // printf("  Deformed\n");

            #pragma omp atomic
            total_cycle_count++;

            // ALIGNMENT
            auto cycle_start = std::chrono::high_resolution_clock::now();  

            // printf("Start idx = %ld\n", (tg.is_rcg) ? tg.start_idx + read_start - genome_length : tg.start_idx + read_start);
            // printf("  Start\n");

            // Alignment_pair alignment = test_on_string(grid, mg, read1, read2);
            
            auto alignment = inferencer(read1.begin(), read1.end(), read2.begin(), read2.end());
            // printf("  start_ids ");
            // print_vector(alignment.start_ids);
            // inferencer.search_pairend(alignment, read2.begin(), read2.end());
            // printf("  After paired ");
            // print_vector(alignment.start_ids);
            

            // printf("  kmer_n = %lu, aligned %lu parts\n", read1.size() / kmer_size, alignment.read1.aligned_num_parts);
            
            // printf("  Stop. Start ids size = %lu\n", alignment.start_ids.size());
            
            auto cycle_end = std::chrono::high_resolution_clock::now();
            // printf("\n");
            std::chrono::duration<double, std::milli> fp_ms = cycle_end - cycle_start;            
            cycle_time += fp_ms.count();
            
            if (alignment.start_ids.empty())
            {
                // printf("Did not align\n\n");
                number_of_nonaligned++;
                continue;
            }
            // printf("\n");
            
            // printf("  Getting genes\n");
            auto pred_gene_ids = get_gene_ids(mg, alignment);
            // printf("  Done\n");
            
            // bool found = false;
            for (int64_t gene_idx : pred_gene_ids)
            {
                if (mg.names[gene_idx] == tg.name)
                {
                    valid_results++;
                    // found = true;
                    break;
                }
            }

            // if (not found)
            // {
            //     printf("Did not align well.\n");
            // }

            // printf("\n");            
            
            #pragma omp master
            if (total_cycle_count >= next_print_barrier)
            {
                next_print_barrier += print_offset;
                print_progress(total_cycle_count, repeat_n);
            }

            read1.clear();
            read2.clear();
        }

        #pragma omp atomic
        total_number_of_valid += valid_results;
        #pragma omp atomic
        total_cycle_time += cycle_time;
        #pragma omp atomic
        total_number_of_overlaps += number_of_overlaps;
        #pragma omp atomic
        total_number_of_nonaligned += number_of_nonaligned;
    }

    auto total_time_stop = std::chrono::high_resolution_clock::now();  
    std::chrono::duration<double, std::milli> fp_ms = total_time_stop - total_time_start;
    double total_duration = fp_ms.count();

    const size_t total_number_of_invalid = total_cycle_count - total_number_of_valid - total_number_of_overlaps - total_number_of_nonaligned;
    const double valid_estimate = total_number_of_valid / (double) total_cycle_count * 100;
    const double invalid_estimate = total_number_of_invalid / (double) total_cycle_count * 100;
    const double overlap_estimate = total_number_of_overlaps / (double) total_cycle_count * 100;
    const double nonaligned_estimate = total_number_of_nonaligned / (double) total_cycle_count * 100;

    clean_previous_console_line();
    printf("Testing Done!\n");
    printf("Total number of aligments        = %lu (reads total %lu)\n", total_cycle_count, 2*total_cycle_count);
    printf("Valid aligments                  = %lu (%.2lf %%)\n", total_number_of_valid, valid_estimate);
    printf("Invalid alignments               = %lu (%.2lf %%)\n", total_number_of_invalid, invalid_estimate);
    printf("Overlapped alignments            = %lu (%.2lf %%)\n", total_number_of_overlaps, overlap_estimate);
    printf("Non-aligned alignments           = %lu (%.2lf %%)\n", total_number_of_nonaligned, nonaligned_estimate);
    printf("Average read length              = %.1f chars.\n\n", ((double) read_len_sum) / total_cycle_count);

    printf("Time summary:\n");
    printf("Number of alignments  = %lu\n", total_cycle_count);
    printf("Total cycle time      = %.3lfs (%.3lfms/thread)\n", total_cycle_time / 1000.0, total_cycle_time / num_threads);
    printf("Average cycle time    = %.4lfms/read/thread\n", total_cycle_time / total_cycle_count / num_threads);
    printf("Total time            = %.3lfs (Includes time for generating reads and other)\n\n", total_duration/1000);
}