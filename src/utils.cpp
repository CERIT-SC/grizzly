/**
 * @file utils.cpp
 * @author Jan Šamánek (jansamanek@email.cz)
 * @brief See the header file.
 * @version 1.1
 * @date 2022-06-30
 * @copyright Copyright (c) 2022
 * This software is open source. There is no license to use this software. Do whatever you want with it.
 * [DISCLAIMER] code may not satisfy clean-code requirements. Scroll down on your own risk.
 */

#include "utils.hpp"
#include "config.hpp"
#include "kmer_converter.hpp"
#include "inferencer.hpp"

void clean_console_line()
{
    printf("\33[2K");
}

void clean_previous_console_line()
{
    printf("\033[A");
    clean_console_line();
}

void print_progress(const size_t progress, const size_t size, const string msg, const string end)
{
	float ratio = ((float) progress) / size;
	int eq_count = ratio * 10 - 1;
    eq_count = eq_count >= 0 ? eq_count : 0;
	string eq(eq_count, '=');
	string space((9 - eq_count), ' ');
    string arrow(ratio >= 0.1 ? ">" : " ");
	printf("[%s%s%s] %s%.2f %%%s%s%s", eq.c_str(), arrow.c_str(), space.c_str(), (ratio < 0.1) ? " " : "", ratio * 100, msg.empty() ? " - " : "", msg.c_str(), end.c_str());
}

size_t get_file_size(istream &file)
{
    // Stop eating new lines in binary mode!
    // file.unsetf(std::ios::skipws);

    // get its size:
    streampos fileSize;

    file.seekg(0, std::ios::end);
    fileSize = file.tellg();
    file.seekg(0, std::ios::beg);
    return fileSize;
}

/**
 * @brief Alter character to reverse complement.
 * 
 * @param c Nucleotide.
 * @return char - Reverse complement.
 */
char alter_char(char c)
{
    switch (c){
        case 'a':
            return 't';
        case 'A':
            return 'T';
        case 'c':
            return 'g';
        case 'C':
            return 'G';
        case 'g':
            return 'c';
        case 'G':
            return 'C';
        case 't':
            return 'a';
        case 'T':
            return 'A';
        default:
            return c;
    }
}

string reverse_complement(string::const_iterator bit, string::const_iterator eit)
{
    string rc;
    rc.reserve(eit-bit);
    for (string::const_iterator it = eit-1; it >= bit; it--)
    {
        rc.push_back(alter_char(*it));
    }

    return rc;
}

boost::iostreams::mapped_file open_file(const char* filename)
{
    boost::iostreams::mapped_file mf;
    try
    {
        // Open input file
        boost::iostreams::mapped_file_params params;
        params.path = filename;
        params.flags = boost::iostreams::mapped_file::mapmode::readonly;
        mf.open(params);
    }
    catch(const std::exception& e){};

    return mf;
}

void parse_args(int argc, const char* argv[])
{
    po::options_description desc("Arguments:");
    desc.add_options()
        ("help,h", "produce help message\n")
        // ("genome_file,g", po::value<string>(&genome_file)->default_value(""), "fasta file to create the grid from. *This file is required for creating new grid.")
        ("gene_annot_file,e", po::value<string>(&gene_file)->default_value(""), "gtf file format containing the information about the reference genome. *This file is required for creating new grid.\n")
        ("load_grid_file,l", po::value<string>(&load_grid_file)->default_value(""), "Path to the saved grid structure.\n")
        ("save_grid_file,s", po::value<string>(&save_grid_file)->default_value(""), "Path to save the grid structure.\n")
        ("alignment_files,a", po::value< vector<string> >(&alignment_files), "Fasta/Fastq file paths which contains reads to be aligned. (Can specify multiple.)\n")
        ("database_file,d", po::value<string>(&database_file)->default_value(""), "Database of changes in the reference genome. For file structure please see README details.\n")
        ("kmer_size,k", po::value<size_t>(&kmer_size)->default_value(14), "kmer size. (Overwritten if loading grid.)\n")
        ("max_n,n", po::value<size_t>(&max_n_num)->default_value(2), "Maximum number of of 'N' characters for 1 kmer. (Overwritten if loading grid.)\n")
        ("max_splice_junction,j", po::value<size_t>(&max_jump_val)->default_value(1000), "Maximum gap between 2 kmers during alignment that is considered as splice juction. (Overwritten if loading grid.)\n")
        ("num_threads,t", po::value<size_t>(&num_threads)->default_value(1), "Number of threads. (range <1 - cpu number of threads>)\n")
        ("minimum_align,i", po::value<float>(&minimum_align)->default_value(0.0), "Specify, in procentage <0-1>, minimum alignment needed to count the read alignment as valid.\n")
        ("progress_frequence,f", po::value<float>(&print_progress_frequence)->default_value(0.01), "Frequency of printing progress. (range <0.001-1.0>)\n")
        ("soft_mask_bases,m", po::bool_switch(&include_softmask_bases)->default_value(false), "Include soft-masked bases (lowercase) from genome when creating grid. This can significantly harm the performance. (false, if not specified)\n")
        ("perform_tests,p", po::bool_switch(&perform_tests)->default_value(false), "If no alignment files were specified. Do performance testing of the grid instead. Gene annotation file is required. (false, if not specified)\n")
        ("do_pair_end,r", po::bool_switch(&do_pair_end_alignment)->default_value(false), "Specify whether you want to do single-end or pair-end alignment. (single-end, if not specified)\n");
        ;

    po::positional_options_description p;
    p.add("alignment_files", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
            options(desc).positional(p).run(), vm);
    po::notify(vm);

    bool fault = false;
    if (gene_file.empty() && load_grid_file.empty() and not vm.count("help"))
    {
        cerr << endl << "ERROR: It is required to specify one of the following options:\n - 'gene_annot_file,e'\n - 'load_grid_file,l'\n Otherwise there is nothing to do." << endl << endl;
        // cout << desc << endl;
        fault = true;

    }
    
    if (vm.count("help") or fault)
    {
        cout << desc << endl;
        cout << "Every other argument will be interpreted as a alignment file." << endl;
        cout << "Every histogram for each alignment file will be saved to the same path with '.histogram' suffix." << endl;
        exit(0);
    
    }
}

bool check_gene_belonging(const vector<Gene_s> &genes, size_t gene_idx, int64_t start_idx, int64_t end_idx)
{
    return (genes[gene_idx].start_idx <= start_idx and genes[gene_idx].end_idx >= end_idx);
}

bool is_mainly_masked(string::const_iterator bit, string::const_iterator eit)
{
    uint32_t lower = 0;
    uint32_t upper = 0;

    for (; bit < eit; bit++)
    {
        if (islower(*bit))    { lower++; } 
        else if (*bit != 'N') { upper++; }
    }

    return lower > upper;
}

MergedGenes_s merge_genes(const GenesAnnotation_s &ga)
{
    MergedGenes_s mg;
    const size_t size = ga.fwgs.size();
    mg.start_ids.reserve(size);
    mg.end_ids.reserve(size);
    mg.names.reserve(size);

    for (const Gene_s& gene : ga.fwgs)
    {
        mg.start_ids.push_back(gene.start_idx);
        mg.end_ids.push_back(gene.end_idx);
        mg.names.push_back(gene.name);
    }
    
    // Sort vectors by start indices
    auto sorted_ids = argsort(mg.end_ids);
    sort_by_indices(mg.start_ids, sorted_ids);
    sort_by_indices(mg.end_ids, sorted_ids);
    sort_by_indices(mg.names, sorted_ids);
    // sort(mg.start_ids.begin(), mg.start_ids.end());

    return mg;
}

unordered_set<int64_t> get_gene_ids(const MergedGenes_s &mg, const Alignment &alignment)
{       
    unordered_set<int64_t> ids;

    if (alignment.start_ids.empty())
    {
        return ids;
    }

    for (int64_t read_start : alignment.start_ids)
    {

        // int64_t read_end = read_start + alignment.read_size;

        // if (read_start < 0) { read_start += genome_length; }
        // else                { read_end   += genome_length; }

        read_start += (read_start < 0) ? genome_length : 0;

        size_t read_gene_idx = lower_bound(mg.end_ids.begin(), mg.end_ids.end(), read_start) - mg.end_ids.begin();
        ids.insert(read_gene_idx);

        // if  ((read_end > mg.end_ids[read_gene_idx] or 
        //     (read_gene_idx+1 < mg.end_ids.size() and mg.start_ids[read_gene_idx+1] < read_start)) and
        //     mg.names[read_gene_idx] != mg.names[read_gene_idx+1])
        // {
        //     ids.insert(read_gene_idx+1);
        // }
    }

    return ids;
}

void paralel_reduction(Vec2D<float> &m)
{
    const size_t row_size = m[0].size();
    size_t thread_idx = omp_get_thread_num();

    for (size_t stride = 2; stride/2 < num_threads; stride *= 2)
    {

        // Using integer division
        const size_t row_to_idx = (thread_idx / stride) * stride;
        const size_t row_from_idx = row_to_idx + stride/2;
        const size_t part_idx = thread_idx % stride;
        const size_t part_size = row_size / stride;

        const size_t start = part_idx*part_size;
        const size_t end = (part_idx == (stride-1) or is_last_thread(thread_idx)) ? 
                            row_size : 
                            (part_idx+1) * part_size;
        
        if (row_from_idx < m.size())
        {
            #pragma omp simd
            for (size_t i = start; i < end; i++)
            {
                m[row_to_idx][i] += m[row_from_idx][i];
            }
        }
        #pragma omp barrier
    }
}

void statistical_magic(vector<float> &hard, vector<float> &soft, double& total_distribution_sum, float hard_sum, float soft_sum)
{    
    vector<float> distribution(hard.size(), 0.f);
    float thread_distribution_sum = 0.0;

    size_t thread_idx = omp_get_thread_num();
    size_t thread_block = hard.size() / num_threads;
    size_t start = thread_block * thread_idx;
    size_t end = (is_last_thread(thread_idx)) ? hard.size() : start + thread_block;

    #pragma omp simd reduction(+ : thread_distribution_sum)
    for(size_t i = start; i < end; i++){
        distribution[i] = (hard[i] / hard_sum) * (soft[i]);

        thread_distribution_sum += distribution[i];
    }

    #pragma omp atomic
    total_distribution_sum += thread_distribution_sum;

    #pragma omp barrier

    auto normalization_factor = soft_sum / total_distribution_sum;

    #pragma omp simd
    for(size_t i = start; i < end; i++){
        hard[i] += distribution[i] * normalization_factor; 
    }
}