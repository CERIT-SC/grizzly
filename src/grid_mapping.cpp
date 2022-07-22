/**
 * @file grid_mapping.cpp
 * @author Jan Šamánek (jansamanek@email.cz)
 * @brief File containing main function for grizzly software.
 * @version 1.1
 * @date 2022-06-30
 * 
 * @copyright Copyright (c) 2022
 * This software is open source. There is no license to use this software. Do whatever you want with it.
 * [DISCLAIMER] code may not satisfy clean-code requirements. Scroll down on your own risk.
 */

#include "config.hpp"
#include "grid_mapping.hpp"
#include "testing.hpp"
#include "utils.hpp"
#include "extend_grid.hpp"
#include "alignment_parser.hpp"

using namespace std;

template <typename T>
void correct_variable(T& var, T lower_bound, T upper_bound, const string& var_name)
{
    if (var < lower_bound)
    {
        cout << "[WARNING]: " << var_name << " can be only in range <" << lower_bound << ", " << upper_bound << ">. Using lower bound instead." << endl;
        var = lower_bound;
        return;
    
    } else if (var > upper_bound)
    {
        cout << "[WARNING]: " << var_name << " can be only in range <" << lower_bound << ", " << upper_bound << ">. Using upper bound instead." << endl;
        var = upper_bound;
    }
}

int main(int argc, char const *argv[])
{
    // Setting console settings to not to be buffered.
    // It is for print_progress function
    setvbuf(stdout, NULL, _IONBF, 0);
    
    // First make sure that the given variables are in terms of given boundaries
    parse_args(argc, argv);
    correct_variable(kmer_size, 8UL, 16UL, "Kmer size");
    correct_variable(max_n_num, 0UL, kmer_size, "Max N bases in kmer");
    correct_variable(max_jump_val, 0UL, 1000000UL, "Max splice juction");
    correct_variable(num_threads, 1UL, (size_t) omp_get_max_threads(), "Number of threads");
    correct_variable(print_progress_frequence, 0.001f, 1.f, "Print progress frequency");
    correct_variable(minimum_align, 0.f, 1.f, "Minimum alignment");
    omp_set_num_threads(num_threads);
    perform_tests = perform_tests and alignment_files.empty() and (not gene_file.empty());

    printf("CONFIG. PRINT\n");
    printf("Kmer size             = %lu\n", kmer_size);
    printf("Max. N chars          = %lu\n", max_n_num);
    printf("Max. split juction    = %lu\n", max_jump_val);
    printf("Annotation file       = %s\n", gene_file.c_str());
    printf("Load file             = %s\n", load_grid_file.c_str());
    printf("Save file             = %s\n", save_grid_file.c_str());
    printf("Database file         = %s\n", database_file.c_str());
    printf("Number of threads     = %lu\n", num_threads);
    printf("Minimum alignment     = %.3f %%\n", minimum_align * 100);
    printf("Print progress freq.  = %.3f\n", print_progress_frequence);
    printf("Include softmasked b. = %s\n", include_softmask_bases ? "Yes" : "No");
    printf("Do pair end alignment = %s\n", do_pair_end_alignment ? "Yes" : "No");
    printf("Do performance tests  = %s\n", (perform_tests) ? "Yes" : "No");
    printf("Aligment files: ");
    print_vector(alignment_files);
    printf("-----------------------------------------------------\n\n");

    Grid<int64_t> grid;
    ifstream annotation_fd(gene_file);
    GenesAnnotation_s ga;

    uint32_t cell_n = (1<<(kmer_size*2));
    grid.resize(cell_n);

    // Creating new grid
    if (load_grid_file.empty())
    {
        uint32_t cell_n = (1<<(kmer_size*2));
        grid.resize(cell_n);

        printf("Creating new grid\n");

        if (not annotation_fd.good())
        {
            fprintf(stderr, "Could not load annotation file \"%s\"\n", gene_file.c_str());
            return 1;
        }

        printf("Loading annotation...\n");

        size_t nucl_count = load_transcriptome(annotation_fd, ga);
        genome_length = ga.fwgs.back().end_idx;

        clean_console_line();
        clean_previous_console_line();
        printf("Annotation loaded successfully.\n");
        printf(" - Number of forward genes   = %lu\n", ga.fwgs.size());
        printf(" - Number of rev. com. genes = %lu\n", ga.rcgs.size());
        printf(" - Total number of bases     = %lu\n", nucl_count);

        
        fill_grid_transcriptome(grid, ga);
        
        annotation_fd.close();
        printf("-----------------------------------------------------\n\n");
    
    // Loading grid from file
    } else 
    {
        printf("Loading grid from \"%s\"\n", load_grid_file.c_str());

        load_grid(load_grid_file, &grid, ga);
        
        printf("Grid loaded\n");
        printf("-----------------------------------------------------\n\n");
    }

    printf("Creating gene pool..\n");
    auto mg = merge_genes(ga);
    clean_console_line();
    clean_previous_console_line();

    // Extending grid
    if (not database_file.empty())
    {
        printf("Extending grid by database from \"%s\".\n", database_file.c_str());
        if (gene_file.empty())
        {
            fprintf(stderr, "Gene annotation file is required for extension\n Skipping extending grid\n");

        } else
        {
            printf("Loading annotation\n");
            ifstream fd(gene_file);
            load_transcriptome(fd, ga);
            fd.close();
            clean_previous_console_line();
        
            fd.open(database_file);
            if (not fd.good())
            {
                fprintf(stderr, "Could not open the databse file from \"%s\" Exiting with value 2.\n", database_file.c_str());
                // mf.close();
                exit(2);
            }

            extend_grid(grid, ga, fd);
            fd.close();     
        }
        
        printf("Grid extended\n");
        printf("-----------------------------------------------------\n\n");
    }

    // Saving grid
    if (not save_grid_file.empty())
    {
        printf("Saving grid to \"%s\"\n", save_grid_file.c_str());
        save_grid(&grid, ga);
        printf("Grid saved to file: %s\n", save_grid_file.c_str());
        printf("-----------------------------------------------------\n\n");
    }
    

    // Does alignment
    if (not alignment_files.empty())
    {
        
        printf("Parsing aligments files\n");
        printf("Doing %s alignments\n\n", (do_pair_end_alignment) ? "pair-end" : "single-end");

        if (do_pair_end_alignment) {
            if (alignment_files.size() % 2 != 0){
                printf("Un-equal number of forward and reverse files!\n");
                return 1;
            }
            size_t i_p = 1;
            for (size_t i = 0; i < alignment_files.size(); i += 2) {
                printf("[%lu/%lu] - %s\n", i_p++, alignment_files.size()/2, alignment_files[i].c_str());
                Alignment_parser::align_file(alignment_files[i], alignment_files[i+1], grid, mg);
            }
        }
        else{
            size_t i = 1;
            for (auto align_file : alignment_files)
            {
                printf("[%lu/%lu] - %s\n", i++, alignment_files.size(), align_file.c_str());
                Alignment_parser::align_file(align_file, grid, mg);
            }
        }
        clean_previous_console_line();
        printf("-----------------------------------------------------\n\n");

    } else if (perform_tests)
    {
        srand(time(NULL));
        const uint32_t cycle_n = 10000000;
        const size_t min_read_size = 80;
        const size_t max_read_size = 120;
        const bool deform = true;

        printf("Nothing to align. performing tests instead.\nTesting config. print:\n");
        printf(" - Number of threads = %lu\n", num_threads);
        printf(" - Number of cycles = %u,\n", cycle_n);
        printf(" - Max number of N chars = %lu,\n", max_n_num);
        printf(" - Artificial deform of read = %s,\n", (deform) ? "true" : "false");
        printf(" - Min. read size %lu, Max. read size %lu,\n", min_read_size, max_read_size);
        printf(" - Doing %s alignments\n", (do_pair_end_alignment) ? "pair-end" : "single-end");
        
        if (do_pair_end_alignment)
        {
            test_inference_speed_on_genome_pairend(grid, ga, mg, cycle_n, max_n_num, deform, min_read_size, max_read_size);
        
        } else
        {
            test_inference_speed_on_genome(grid, ga, mg, cycle_n, max_n_num, deform, min_read_size, max_read_size);
        }    
    }
    
    return 0;
}
