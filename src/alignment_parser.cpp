/**
 * @file alignment_parser.cpp
 * @author Jan Šamánek (jansamanek@email.cz)
 * @brief See header file
 * @version 1.1
 * @date 2022-06-30
 * 
 * @copyright Copyright (c) 2022
 * No copyright.
 */

#include "alignment_parser.hpp"
#include "config.hpp"
#include "utils.hpp"
#include "inferencer.hpp"

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>
#include <omp.h>

namespace Alignment_parser
{

/**
 * @brief Skipping line in file
 * 
 * @param bytes Mapped file in char* form.
 * @param i index.
 */
void skip_line(const char* bytes, size_t& i, const size_t end)
{
    while(i < end and not is_linefeedchar(bytes[i++]));
}

/**
 * @brief Get the read from char*
 * 
 * @param bytes mapped file in char* form
 * @param i index
 * @param is_fastq flag for parsing purpuses 
 * @return const string - Read
 */
const string get_read(const char* bytes, size_t& i, const bool is_fastq, const size_t end)
{
    string read;
    skip_line(bytes, i, end);

    // Parse actual sequence
    char nucl;
    while (i < end and not is_linefeedchar((nucl = bytes[i++])))
    {
        if (nucl == 'N')
        {
            // Using integer division to get index of the kmer
            // And multipling again to get the offset
            size_t start = (read.length() / kmer_size) * kmer_size;
            // Count N chars in last kmer
            // If number of N exceedes maximum number allowed in one kmer, than stop loading read
            if ((size_t) count(read.begin() + start, read.end(), 'N')+1 > max_n_num)
            {
                skip_line(bytes, i, end);
                break;
            }
        }
        read.push_back(nucl);
    }

    if (is_fastq)
    {
        skip_line(bytes, i, end);
        skip_line(bytes, i, end);
    }

    return read;
}

/**
 * @brief Utility function that handles input and output files
 * 
 * @param input_file Input alignment file
 * @param input_mf Mapped file. Filled via reference
 * @param output_fd Output histogram file descriptor. Filled via reference
 */
void open_files(const string input_file, boost::iostreams::mapped_file &input_mf, ofstream &output_fd, bool &open_fine)
{
    input_mf = open_file(input_file.c_str());

    if (not input_mf.is_open())
    {
        fprintf(stderr, "Could not open input file \"%s\". Skipping.\n", input_file.c_str());
        open_fine = false;
        return;
    }

    string out_filename = input_file + ".histogram";
    output_fd.open(out_filename.c_str());
    if (not output_fd.good())
    {
        fprintf(stderr, "Could not create output file for input \"%s\". Skipping.", out_filename.c_str());
        input_mf.close();
        open_fine = false;
        return;
    }

    open_fine = true;
}


/**
 * @brief Core function that handles alignment itself.
 * 
 * @param mf Mapped input file
 * @param grid Grid structure
 * @param hard_histograms Matrix of gene counters (gene_n X num_threads). Filled via referece.
 * @param total_num_reads Read counter. Filled via reference.
 * @param total_non_aligned Non-aligned counter. Filled via reference.
 * @param total_read_len_sum Read length counter. Filled via reference.
 * @param is_fastq Flag for input parsing purpuses.
 * @param mg Gene annotation.
 */
void parse( boost::iostreams::mapped_file mf, 
            const Grid<int64_t> &grid,
            Vec2D<float> &hard_histograms,
            size_t &total_num_reads,
            size_t &total_non_aligned,
            size_t &total_read_len_sum,
            size_t &total_multimaps,
            // size_t &total_overlapping_reads,
            const bool is_fastq,
            const MergedGenes_s& mg)
{
    vector<size_t> offsets(num_threads);
    
    const size_t file_size = mf.size();
    const char* bytes = mf.const_data();
    const char delim = (is_fastq) ? '@' : '>';
    size_t bytes_counter = 0;

    double next_print_barrier = 0;
    const double print_offset = file_size * print_progress_frequence;

    #pragma omp parallel default(shared) num_threads(num_threads)
    {
        const size_t thread_idx = omp_get_thread_num();
        // printf("Thread %lu. started\n", thread_idx);

        size_t thread_offset = thread_idx * (file_size/num_threads);
        hard_histograms[thread_idx].resize(mg.names.size());
        size_t thread_read_count = 0;
        size_t thread_non_aligned = 0;
        size_t thread_read_len_sum = 0;
        size_t thread_multimaps = 0;
        // size_t thread_overlapping_reads = 0;
        Inferencer inferencer(grid);

        // Every thread set itself at the beginning of the new record in file
        while(bytes[thread_offset] != delim) { thread_offset++; };

        // Check for fastq files if threads really set themselfs at the beginning of the record
        if (not is_first_thread(thread_idx) and is_fastq)
        {
            size_t thread_offset_copy = thread_offset;
            skip_line(bytes, thread_offset_copy, file_size);
            skip_line(bytes, thread_offset_copy, file_size);
            // If not plus sign then delim was not start of head but part of quality score
            if (bytes[thread_offset_copy] != '+')
            {
                skip_line(bytes, thread_offset, file_size);
            }
        }

        // Find first read for pair end
        if (do_pair_end_alignment)
        {
            size_t thread_offset_copy = thread_offset;
            string header1;
            string header2;

            while (not is_linefeedchar(bytes[thread_offset_copy++])) 
            { 
                header1.push_back(bytes[thread_offset_copy-1]);
            }

            skip_line(bytes, thread_offset_copy, file_size);
            if (is_fastq)
            {
                skip_line(bytes, thread_offset_copy, file_size);
                skip_line(bytes, thread_offset_copy, file_size);
            }

            size_t anchor = thread_offset_copy;
            while (not is_linefeedchar(bytes[thread_offset_copy++])) 
            { 
                header2.push_back(bytes[thread_offset_copy-1]);
            }

            if (header1 != header2)
            {
                thread_offset = anchor;
            }
        }

        // string header1;
        // string header2;
        // size_t thread_offset_copy = thread_offset;
        // while (not is_linefeedchar(bytes[thread_offset_copy++]))
        // {
        //     header1.push_back(bytes[thread_offset_copy-1]);
        // }
        // skip_line(bytes, thread_offset_copy);
        // while (not is_linefeedchar(bytes[thread_offset_copy++]))
        // {
        //     header2.push_back(bytes[thread_offset_copy-1]);
        // }
        // printf("%lu.\n%s\n%s\n%s\n\n", thread_idx, header1.c_str(), header2.c_str(), (header1 == header2) ? "is same" : "differs");

        offsets[thread_idx] = thread_offset;
        #pragma omp barrier

        const size_t start = offsets[thread_idx];
        const size_t end = is_last_thread(thread_idx) ? mf.size() : offsets[thread_idx+1];

        // printf("Thread %lu. start %lu, end %lu\n", thread_idx, start, end);

        string read1, read2;
        read1.reserve(150);
        read2.reserve(150);
        unordered_set<int64_t> gene_ids;

        for (size_t i = start; i < end; i++)
        {
            size_t i_copy = i;
            read1.clear();
            read1 = get_read(bytes, i, is_fastq, end);

            gene_ids.clear();

            Alignment alignment;
            size_t kmer_n = 0;
            if (do_pair_end_alignment)
            {
                read2.clear();
                read2 = get_read(bytes, i, is_fastq, end);

                alignment = inferencer(read1.begin(), read1.end(), read2.begin(), read2.end());
                gene_ids = get_gene_ids(mg, alignment);

                kmer_n = read1.length() / kmer_size + read2.length() / kmer_size;
                thread_read_count++;
                // Integer division
                thread_read_len_sum += (read1.length() / kmer_size * kmer_size) + (read2.length() / kmer_size * kmer_size);

            } else
            {
                alignment = inferencer(read1.begin(), read1.end());
                gene_ids = get_gene_ids(mg, alignment);

                kmer_n = read1.length() / kmer_size;
                thread_read_count++;
                // Integer division
                thread_read_len_sum += read1.length() / kmer_size * kmer_size;
            }

            #pragma omp atomic
            bytes_counter += (i - i_copy);

            if (gene_ids.empty() or kmer_n * minimum_align > alignment.aligned_num_parts)
            {
                thread_non_aligned++;
                continue;

            } else if (gene_ids.size() == 1)
            {
                // if (do_pair_end_alignment)
                // {
                //     printf("%lu. %ld %s\n", thread_read_count, alignment.start_ids.front(), alignment.start_ids.front() < 0 ? to_string(alignment.start_ids.front() + genome_length).c_str() : "");
                
                // } else
                // {
                //     printf("%lu. %lu - %ld %s\n", (thread_read_count-1) / 2 + 1, (thread_read_count-1) % 2 + 1, alignment.start_ids.front(), alignment.start_ids.front() < 0 ? to_string(alignment.start_ids.front() + genome_length).c_str() : "");
                // }
                
                hard_histograms[thread_idx][*gene_ids.begin()]++;
            
            } else
            {
                thread_multimaps++;
                const float frac = 1.f/gene_ids.size();
                for (const int64_t gene_idx : gene_ids)
                {
                    hard_histograms[thread_idx][gene_idx] += frac;
                }
            }

            if (bytes_counter >= next_print_barrier)
            {
                #pragma omp atomic
                next_print_barrier += print_offset;
                print_progress(bytes_counter, file_size);
            }
        }

        #pragma omp atomic
        total_num_reads += thread_read_count;

        #pragma omp atomic
        total_non_aligned += thread_non_aligned;

        #pragma omp atomic
        total_read_len_sum += thread_read_len_sum;

        #pragma omp atomic
        total_multimaps += thread_multimaps;
        
        // printf("Thread %lu. finished\n", thread_idx);
        #pragma omp barrier

        // Parallel reduction
        paralel_reduction(hard_histograms);

        #pragma omp barrier
    }
}

void parse( boost::iostreams::mapped_file mf_f,
            boost::iostreams::mapped_file mf_r,
            const Grid<int64_t> &grid,
            Vec2D<float> &hard_histograms,
            size_t &total_num_reads,
            size_t &total_non_aligned,
            size_t &total_read_len_sum,
            size_t &total_multimaps,
            const bool is_fastq,
            const MergedGenes_s& mg)
{
    vector<size_t> offsets(num_threads);

    const size_t file_size = mf_f.size();
    if (file_size != mf_r.size()){
        cout<<"Different size of files!"<<endl;
        cout<<"Ending alignment."<<endl;
        return;
    }
    const char* bytes = mf_f.const_data();
    const char* bytes_r = mf_r.const_data();
    const char delim = (is_fastq) ? '@' : '>';
    size_t bytes_counter = 0;

    double next_print_barrier = 0;
    const double print_offset = file_size * print_progress_frequence;

    #pragma omp parallel default(shared) num_threads(num_threads)
    {
        const size_t thread_idx = omp_get_thread_num();
        // printf("Thread %lu. started\n", thread_idx);

        size_t thread_offset = thread_idx * (file_size/num_threads);
        hard_histograms[thread_idx].resize(mg.names.size());
        size_t thread_read_count = 0;
        size_t thread_non_aligned = 0;
        size_t thread_read_len_sum = 0;
        size_t thread_multimaps = 0;
        // size_t thread_overlapping_reads = 0;
        Inferencer inferencer(grid);

        // Every thread set itself at the beginning of the new record in file
        while(bytes[thread_offset] != delim) { thread_offset++; };

        // Check for fastq files if threads really set themselfs at the beginning of the record
        if (not is_first_thread(thread_idx) and is_fastq)
        {
            size_t thread_offset_copy = thread_offset;
            skip_line(bytes, thread_offset_copy, file_size);
            skip_line(bytes, thread_offset_copy, file_size);
            // If not plus sign then delim was not start of head but part of quality score
            if (bytes[thread_offset_copy] != '+')
            {
                skip_line(bytes, thread_offset_copy, file_size);
            }
        }

        offsets[thread_idx] = thread_offset;
        #pragma omp barrier

        const size_t start = offsets[thread_idx];
        const size_t end = is_last_thread(thread_idx) ? mf_f.size() : offsets[thread_idx+1];

        // printf("Thread %lu. start %lu, end %lu\n", thread_idx, start, end);

        string read1, read2;
        read1.reserve(150);
        read2.reserve(150);
        unordered_set<int64_t> gene_ids;

        for (size_t i = start; i < end; i++)
        {
            size_t i_copy = i;
            size_t i_r = i;
            read1.clear();
            read1 = get_read(bytes, i, is_fastq, end);

            gene_ids.clear();

            Alignment alignment;
            size_t kmer_n = 0;

            read2.clear();
            read2 = get_read(bytes_r, i_r, is_fastq, end);

            alignment = inferencer(read1.begin(), read1.end(), read2.begin(), read2.end());
            gene_ids = get_gene_ids(mg, alignment);

            thread_read_count++;

            kmer_n = read1.length() / kmer_size + read2.length() / kmer_size;
            // Integer division
            thread_read_len_sum += (read1.length() / kmer_size * kmer_size) + (read2.length() / kmer_size * kmer_size);

            #pragma omp atomic
            bytes_counter += (i - i_copy);

            if (gene_ids.empty() or kmer_n * minimum_align > alignment.aligned_num_parts)
            {
                thread_non_aligned++;
                continue;

            } else if (gene_ids.size() == 1)
            {
                hard_histograms[thread_idx][*gene_ids.begin()]++;
            } else
            {
                thread_multimaps++;
                const float frac = 1.f/gene_ids.size();
                for (const int64_t gene_idx : gene_ids)
                {
                    hard_histograms[thread_idx][gene_idx] += frac;
                }
            }

            // #pragma omp master
            if (bytes_counter >= next_print_barrier)
            {
                #pragma omp atomic
                next_print_barrier += print_offset;
                print_progress(bytes_counter, file_size);
            }
        }

        #pragma omp atomic
        total_num_reads += thread_read_count;

        #pragma omp atomic
        total_non_aligned += thread_non_aligned;

        #pragma omp atomic
        total_read_len_sum += thread_read_len_sum;

        #pragma omp atomic
        total_multimaps += thread_multimaps;

        // printf("Thread %lu. finished\n", thread_idx);
        #pragma omp barrier

        // Parallel reduction
        paralel_reduction(hard_histograms);
    }
}

void align_file(string input_file, 
                const Grid<int64_t> &grid, 
                const MergedGenes_s& mg)
{
    boost::iostreams::mapped_file input_mf;
    ofstream output_fd;

    bool is_opened_fine = true;
    Alignment_parser::open_files(input_file, input_mf, output_fd, is_opened_fine);
    if (not is_opened_fine) { return; }

    const char *bytes = input_mf.const_data();
    const bool is_fastq = bytes[0] == '@';
    size_t total_num_reads = 0;
    size_t total_non_aligned = 0;
    size_t total_read_len_sum = 0;
    size_t total_multimaps = 0;
    // size_t total_overlapping_reads = 0;

    printf("File has been determined as \"%s\"\n", (is_fastq) ? "Fastq": "Fasta");

    auto align_start = std::chrono::high_resolution_clock::now();
    Vec2D<float> hard_histograms(num_threads);
    parse(input_mf, grid, hard_histograms, total_num_reads, total_non_aligned, total_read_len_sum, total_multimaps, is_fastq, mg);
    auto align_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = align_end - align_start;
    const double duration = fp_ms.count();

    clean_console_line();
    clean_previous_console_line();
    printf("Alignment summary:\n");
    printf("Total number of %s%s          = %lu %s\n", (do_pair_end_alignment) ? "fragments" : "reads", (do_pair_end_alignment) ? "" : "    ", total_num_reads, (do_pair_end_alignment) ? string("(total reads " + to_string(total_num_reads*2) + ")").c_str() : "");
    printf("Average efective %s length   %s= %.1lf chars.\n", (do_pair_end_alignment) ? "fragment" : "read", (do_pair_end_alignment) ? "" : "    ", ((double) total_read_len_sum / total_num_reads));
    printf("Number of multi-mapped alignments  = %lu (~%.2lf %%)\n", total_multimaps, ((double) total_multimaps) / total_num_reads * 100);
    printf("Number of non-aligned alignments   = %lu (~%.2lf %%)\n", total_non_aligned, ((double) total_non_aligned) / total_num_reads * 100);
    printf("Total time                         = %.3lfs (~%.3lfms/%s/thread)\n\n", duration / 1000.0, duration * num_threads / total_num_reads, (do_pair_end_alignment) ? "fragment" : "read");

    const string output_file = input_file + ".histogram";
    printf("Saving histogram to \"%s\"\n", output_file.c_str());
    double print_offset = mg.names.size() * print_progress_frequence;
    double next_print_barrier = 0;

    output_fd << "gene\tcount" << endl;
    for (size_t i = 0; i < mg.names.size(); i++)
    {
        output_fd << mg.names[i] << "\t" << to_string(hard_histograms[0][i]) << endl;
        if (i > next_print_barrier)
        {
            next_print_barrier += print_offset;
            print_progress(i, mg.names.size());
        }
    }

    input_mf.close();
    output_fd.close();

    clean_console_line();
    clean_previous_console_line();
    clean_previous_console_line();
    printf("Histogram successfully saved to \"%s\"\n\n", output_file.c_str());
    return;
}

void align_file(string input_file,
                string input_file_r,
                const Grid<int64_t> &grid,
                const MergedGenes_s& mg)
{
    boost::iostreams::mapped_file input_mf;
    boost::iostreams::mapped_file input_mf_r(input_file_r);
    ofstream output_fd;

    bool is_opened_fine = true;
    Alignment_parser::open_files(input_file, input_mf, output_fd, is_opened_fine);
    if (not is_opened_fine) { return; }

    if (not input_mf_r.is_open())
    {
        fprintf(stderr, "Could not open input file \"%s\". Skipping.\n", input_file_r.c_str());
        is_opened_fine = false;
    }

    if (not is_opened_fine) { return; }

    const char *bytes = input_mf.const_data();
    const bool is_fastq = bytes[0] == '@';
    size_t total_num_reads = 0;
    size_t total_non_aligned = 0;
    size_t total_read_len_sum = 0;
    size_t total_multimaps = 0;
    // size_t total_overlapping_reads = 0;

    printf("File has been determined as \"%s\"\n", (is_fastq) ? "Fastq": "Fasta");

    auto align_start = std::chrono::high_resolution_clock::now();
    Vec2D<float> hard_histograms(num_threads);
    parse(input_mf, input_mf_r, grid, hard_histograms, total_num_reads, total_non_aligned, total_read_len_sum, total_multimaps, is_fastq, mg);
    auto align_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = align_end - align_start;
    const double duration = fp_ms.count();

    float mean_frag_length = (float)total_read_len_sum / (float)(total_num_reads - total_non_aligned);

    clean_console_line();
    clean_previous_console_line();
    printf("Alignment summary:\n");
    printf("Total number of %s%s          = %lu %s\n", (do_pair_end_alignment) ? "fragments" : "reads", (do_pair_end_alignment) ? "" : "    ", total_num_reads, (do_pair_end_alignment) ? string("(total reads " + to_string(total_num_reads*2) + ")").c_str() : "");
    printf("Average %s length   %s= %.1lf chars.\n", (do_pair_end_alignment) ? "fragment" : "read", (do_pair_end_alignment) ? "" : "    ", mean_frag_length);
    printf("Number of multi-mapped alignments  = %lu (~%.2lf %%)\n", total_multimaps, ((double) total_multimaps) / total_num_reads * 100);
    printf("Number of non-aligned alignments   = %lu (~%.2lf %%)\n", total_non_aligned, ((double) total_non_aligned) / total_num_reads * 100);
    printf("Total time                         = %.3lfs (~%.3lfms/%s/thread)\n\n", duration / 1000.0, duration * num_threads / total_num_reads, (do_pair_end_alignment) ? "fragment" : "read");

    const string output_file = input_file + ".histogram";
    printf("Saving histogram to \"%s\"\n", output_file.c_str());
    double print_offset = mg.names.size() * print_progress_frequence;
    double next_print_barrier = 0;

    output_fd << "gene\tcount" << endl;
    for (size_t i = 0; i < mg.names.size(); i++)
    {
        output_fd << mg.names[i] << "\t" << to_string(hard_histograms[0][i]) << endl;
        if (i > next_print_barrier)
        {
            next_print_barrier += print_offset;
            print_progress(i, mg.names.size());
        }
    }
    input_mf.close();
    input_mf_r.close();
    output_fd.close();

    clean_console_line();
    clean_previous_console_line();
    clean_previous_console_line();
    printf("Histogram successfully saved to \"%s\"\n\n", output_file.c_str());
}


// End of Namespace
}
