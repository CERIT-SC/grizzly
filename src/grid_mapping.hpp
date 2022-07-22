/**
 * @file grid_mapping.hpp
 * @author Jan Šamánek (jansamanek@email.cz)
 * @brief Component of the Grizzly software.
 * File contains functions for saving/loading/creating grid structures.
 * @version 1.1
 * @date 2022-06-30
 * 
 * @copyright Copyright (c) 2022
 * This software is open source. There is no license to use this software.
 * [DISCLAIMER] code may not satisfy clean-code requirements. Scroll down on your own risk.
 */

#pragma once
#include "config.hpp"
#include "utils.hpp"
#include "kmer_converter.hpp"

/**
 * @brief Function for saving the grid structure on to the grid
 * 
 * @param grid Grid structure.
 * @param ga Gene annotation.
 * @param file Output file.
 */
inline void save_grid(Grid<int64_t> *grid, 
                      const GenesAnnotation_s& ga)
{

    ofstream file(save_grid_file);
    if (not file.good())
    {
        fprintf(stderr, "Could not open output file for saving the grid.\n");
        exit(1);
    }
    

    size_t i = 0;
    const size_t grid_size = grid->size();
    // const size_t chr_offsets_size = chr_offsets.size();
    size_t ids_count = 0;
    auto total_time_start = std::chrono::high_resolution_clock::now();  


    // Write metadata
    //Write file prefix
    file.write((char*) &file_prefix.front(), file_prefix.length());
    // Next write saved genome length
    file.write(reinterpret_cast<const char*>(&genome_length), sizeof(genome_length));
    // Next write KMER size of the grid
    file.write(reinterpret_cast<const char*>(&kmer_size), sizeof(kmer_size));
    // Next write MAX N number per kmer for grid
    file.write(reinterpret_cast<const char*>(&max_n_num), sizeof(max_n_num));
    // Next write MAX JUMP for splice juction
    file.write(reinterpret_cast<const char*>(&max_jump_val), sizeof(max_jump_val));

    // GENES ANNOTATION
    // Writing reverse complement genes
    const size_t rcg_genes_size = ga.rcgs.size();
    file.write(reinterpret_cast<const char*>(&rcg_genes_size), sizeof(rcg_genes_size));
    if (rcg_genes_size > 0)
    {
        // Writing start indices 8B * size
        for (const Gene_s& gene : ga.rcgs)
        {
            file.write(reinterpret_cast<const char*>(&gene.start_idx), sizeof(gene.start_idx));        
        }

        // Writing end indices 8B * size
        for (const Gene_s& gene : ga.rcgs)
        {
            file.write(reinterpret_cast<const char*>(&gene.end_idx), sizeof(gene.end_idx));        
        }
        
        // Writing names of genes
        for (const Gene_s& gene : ga.rcgs)
        {
            string name = gene.name;
            name.push_back('\0');
            file.write((char*) &name.front(), name.length());
        }
    }

    // Writing forward genes
    const size_t fwg_genes_size = ga.fwgs.size();
    file.write(reinterpret_cast<const char*>(&fwg_genes_size), sizeof(fwg_genes_size));
    if (fwg_genes_size > 0)
    {
        // Writing start indices 8B * size
        for (const Gene_s& gene : ga.fwgs)
        {
            file.write(reinterpret_cast<const char*>(&gene.start_idx), sizeof(gene.start_idx));        
        }
        
        // Writing end indices 8B * size
        for (const Gene_s& gene : ga.fwgs)
        {
            file.write(reinterpret_cast<const char*>(&gene.end_idx), sizeof(gene.end_idx));        
        }

        // Writing names of genes
        for (const Gene_s& gene : ga.fwgs)
        {
            string name = gene.name;
            name.push_back('\0');
            file.write((char*) &name.front(), name.length());
        }
    }

    // --------- END of METADATA --------- //
    
    const float print_offset = grid_size * print_progress_frequence;
    float next_print_barrier = 0;

    // size_t bin_idx = 0;
    for (const auto& bin : *grid)
    {
        const size_t bin_size = bin.size();
        ids_count += bin_size;

        // First insert size of bin
        file.write(reinterpret_cast<const char*>(&bin_size), sizeof(bin_size));
        // Than insert bin content
        file.write(reinterpret_cast<const char*>(&bin.front()), bin_size * sizeof(bin.front()));
        // Put anchor for parallel loading

        if (++i >= next_print_barrier)
        {
            next_print_barrier += print_offset;
            print_progress(i, grid_size);
        }
    }

    auto total_time_stop = std::chrono::high_resolution_clock::now();  
    std::chrono::duration<double, std::milli> fp_ms = total_time_stop - total_time_start;

    clean_console_line();

    printf("Saved grid stats:\n");
    printf(" - Kmer size            = %lu\n", kmer_size);
    printf(" - Max. N per kmer      = %lu\n", max_n_num);
    printf(" - Max. splice junction = %lu\n", max_jump_val);
    printf(" - Genome length        = %lu\n", genome_length);
    printf(" - Genes annotation sizes:\n");
    printf("   * rev. com. genes = %lu\n", ga.rcgs.size());
    printf("   * forward genes   = %lu\n", ga.fwgs.size());
    printf(" - Total number of cells in the grid = %lu\n", grid->size());
    printf(" - Total number of kmers in the grid = %lu\n", ids_count);
    printf(" - Total time = %.2lfs\n", fp_ms.count() / 1000);
    file.close();
    // printf(" - Chromozome offsets size: %lu\n", chr_offsets.size());
}

/**
 * @brief Function to load the grid structure from the disk
 * 
 * @param file Input file.
 * @param grid Grid structure. (filled via pointer)
 * @param ga Gene annotation. (filled via reference)
 */
inline void load_grid(  const string& file, 
                        Grid<int64_t> *grid, 
                        // unordered_map<string, size_t>& chr_offsets, 
                        GenesAnnotation_s &ga)
{
    auto mf = open_file(file.c_str());
    if (not mf.is_open())
    {
        fprintf(stderr, "Could not open file \"%s\" Ending with -1.\n", file.c_str());
        exit(-1);
    }

    const char* bytes = mf.const_data();

    size_t bin_n = 0;
    size_t total_kmer_count = 0;
    size_t total_bin_count = 0;
    size_t i = 0;

    double print_offset;
    double next_barrier = 0;

    auto total_time_start = std::chrono::high_resolution_clock::now();  

    // ----------- START OF METADATA --------------- //
    // Read file prefix
    string loaded_prefix;
    while (i < file_prefix.length())
    {
        char c = recast<char>(bytes, i);
        loaded_prefix.push_back(c);
    }

    // Check file prefix
    if (loaded_prefix != file_prefix)
    {
        // printf("Faulty prefix = %s\n", loaded_prefix.c_str());
        fprintf(stderr, "File prefix does not match. Are you sure that this is the right file?\n");
        mf.close();
        exit(1);
    }

    // Then load saved genome length
    genome_length = recast<size_t>(bytes, i);
    // Load kmer size
    kmer_size = recast<size_t>(bytes, i);
    // Load max n number
    max_n_num = recast<size_t>(bytes, i);
    // next load max jump splice junction
    max_jump_val = recast<size_t>(bytes, i);

    // Loading gene info
    // First loading rc genes
    size_t rcg_genes_size = recast<size_t>(bytes, i);
    ga.rcgs.resize(rcg_genes_size);
    // Loading start indices
    size_t j = 0;
    size_t end = i + rcg_genes_size*sizeof(int64_t);
    while (i < end)
    {
        const int64_t gene_start_index = recast<int64_t>(bytes, i);
        ga.rcgs[j].start_idx = gene_start_index;
        ga.rcgs[j].is_rcg = true;
        j++;
    }

    // Loading end indices
    j = 0;
    end = i + rcg_genes_size*sizeof(int64_t);
    while (i < end)
    {
        const int64_t gene_end_index = recast<int64_t>(bytes, i);
        ga.rcgs[j].end_idx = gene_end_index;
        j++;
    }

    // Loading names
    char c;
    j = 0;
    while (j < rcg_genes_size)
    {
        c = recast<char>(bytes, i);
        while (c != '\0')
        {
            ga.rcgs[j].name.push_back(c);
            c = recast<char>(bytes, i);
        }
        j++;
    }

    // Now loading fw genes
    size_t fwg_genes_size = recast<size_t>(bytes, i);
    ga.fwgs.resize(fwg_genes_size);
    j = 0;
    end = i + rcg_genes_size*sizeof(int64_t);
    while (i < end)
    {
        const int64_t gene_start_index = recast<int64_t>(bytes, i);
        ga.fwgs[j].start_idx = gene_start_index;
        ga.fwgs[j].is_rcg = false;
        j++;
    }

    // Loading end indices
    j = 0;
    end = i + rcg_genes_size*sizeof(int64_t);
    while (i < end)
    {
        const int64_t gene_end_index = recast<int64_t>(bytes, i);;
        ga.fwgs[j].end_idx = gene_end_index;
        j++;
    }

    // Loading names
    j = 0;
    while (j < fwg_genes_size)
    {
        c = recast<char>(bytes, i);
        while (c != '\0')
        {
            ga.fwgs[j].name.push_back(c);
            c = recast<char>(bytes, i);
        }
        j++;
    }

    //---------------- END OF METADATA --------------------//

    bin_n = (1<<(kmer_size*2));
    grid->clear();
    grid->resize(bin_n);

    print_offset = bin_n * print_progress_frequence;

    size_t bin_idx = 0;
    while(i < mf.size())
    {
        size_t bin_size = recast<size_t>(bytes, i);
        grid->at(bin_idx).reserve(bin_size);
        
        for (size_t j = 0; j < bin_size; j++)
        {
            int64_t index = recast<int64_t>(bytes, i);
            grid->at(bin_idx).push_back(index);
        }

        bin_idx++;
        total_bin_count++;
        total_kmer_count += bin_size;

        if (total_bin_count > next_barrier)
        {
            next_barrier += print_offset;
            print_progress(total_bin_count, bin_n);
        }
    }
    
    auto total_time_stop = std::chrono::high_resolution_clock::now();  
    std::chrono::duration<double, std::milli> fp_ms = total_time_stop - total_time_start;

    clean_console_line();
    printf("Loaded grid stats:\n");
    printf(" - Kmer size            = %lu\n", kmer_size);
    printf(" - Max. N per kmer      = %lu\n", max_n_num);
    printf(" - Max. splice junction = %lu\n", max_jump_val);
    printf(" - Genome length        = %lu\n", genome_length);
    printf(" - Genes annotation sizes:\n");
    printf("   * rev. com. genes = %lu\n", ga.rcgs.size());
    printf("   * forward genes   = %lu\n", ga.fwgs.size());
    printf(" - Total number of cells in the grid = %lu\n", grid->size());
    printf(" - Total number of kmers in the grid = %lu\n", total_kmer_count);
    printf(" - Total time = %.2lfs\n", fp_ms.count() / 1000);
    // printf(" - Chromozome offsets size: %lu\n", chr_offsets.size());
    mf.close();
}

/**
 * @brief Function for creating new grid
 * 
 * @param grid Grid structure. (filled via reference)
 * @param ga Loaded genes annotation (function will clear the gene.gene properties).
 * @return size_t - Returns number of inserted kmers to the grid.
 */
inline size_t fill_grid_transcriptome(  Grid<int64_t>& grid, 
                                        GenesAnnotation_s& ga)
{
    printf("Filling up the grid...\n");
    size_t kmer_count = 0;
    size_t skipped_genes = 0;
    KmerConverterManager kcm;

    const double size = ga.rcgs.size() + ga.fwgs.size();
    const double print_offset = size * print_progress_frequence;
    double next_barrier = 0;
    double progress = 0;
    auto total_time_start = std::chrono::high_resolution_clock::now();  

    for (vector<Gene_s> *genes : {&ga.rcgs, &ga.fwgs})
    {
        for (Gene_s &gene : *genes)
        {
            if (gene.gene.length() < kmer_size) { progress++; skipped_genes++; continue; }

            for (auto it = gene.gene.begin(); it <= gene.gene.end() - kmer_size; ++it)
            {
                if ((size_t) count(it, it + kmer_size, 'N') > max_n_num or
                    (not include_softmask_bases and is_mainly_masked(it, it + kmer_size)))
                { continue; }

                auto grid_indices = kcm(it, it + kmer_size);
                kmer_count += grid_indices.size();

                for (size_t grid_idx : grid_indices)
                {
                    int64_t kmer_idx = gene.start_idx + (it - gene.gene.begin());
                    grid[grid_idx].push_back(gene.is_rcg ? kmer_idx - genome_length : kmer_idx);
                }
            }

            gene.gene.clear();
            gene.gene.shrink_to_fit();

            if (++progress > next_barrier)
            {
                next_barrier += print_offset;
                print_progress(progress, size);
            }
        }   
    }

    auto total_time_stop = std::chrono::high_resolution_clock::now();  
    std::chrono::duration<double, std::milli> fp_ms = total_time_stop - total_time_start;

    clean_console_line();
    clean_previous_console_line();
    printf("Created new grid.\n");
    printf("Total number of skipped genes  = %lu (out of %lu)\n", skipped_genes/2, ga.fwgs.size());
    printf("Total number of inserted kmers = %lu\n", kmer_count);
    printf("Total time = %.2lfs\n", fp_ms.count() / 1000);

    return kmer_count;
}

/**
 * @brief Function for loading the gene annotation from the input file
 * 
 * @param annotation_fd File descriptor for the annotation file.
 * @param ga Gene annotation structure. (filled via reference)
 * @return size_t - Return total number of bases in the file.
 */
inline size_t load_transcriptome(ifstream &annotation_fd, GenesAnnotation_s &ga)
{
    
    size_t last_idx = 1;
    size_t nucl_count = 0;

    string gene;

    while (annotation_fd.good())
    {
        string line;
        getline(annotation_fd, line);
        if (line.empty()) { continue; }
        
        if (line[0] == '>')
        {
            if (not gene.empty())
            {
                ga.fwgs.back().end_idx += gene.length();
                ga.rcgs.back().end_idx += gene.length();

                ga.fwgs.back().is_empty = false;
                ga.rcgs.back().is_empty = false;

                ga.fwgs.back().gene = gene;
                gene = reverse_complement(gene.begin(), gene.end());
                ga.rcgs.back().gene = gene;

                nucl_count += 2*gene.length();
                last_idx = ga.fwgs.back().end_idx + max_jump_val + 1;

                gene.clear();
            }

            stringstream buffer(line);
            string temp;
            vector<string> values;
            while (getline(buffer, temp, ' ')) 
            {
                values.push_back(temp);
            }

            values[0].erase(0, 1);

            int64_t g_start_idx = last_idx;
            int64_t g_end_idx = g_start_idx;
            string transcript_name = values[0];

            ga.fwgs.push_back(Gene_s(g_start_idx, g_end_idx, transcript_name, false));
            ga.rcgs.push_back(Gene_s(g_start_idx, g_end_idx, transcript_name, true));
        
        } else
        {
            gene += line;
        }
    }

    if (ga.fwgs.empty() and ga.rcgs.empty())
    {
        fprintf(stderr, "Did not manage to load single gene annotation. Ending with err. code 3.\n");
        annotation_fd.close();
        exit(3);
    }

    return nucl_count;
}