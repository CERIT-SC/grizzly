/**
 * @file alignment_parser.hpp
 * @author Jan Šamánek (jansamanek@email.cz)
 * @brief Component of the Grizzly software.
 * This component handles alignment on the given Fastq & Fasta files 
 * and produces histogram of genes in given file.
 * Histogram will be saved with same name as the input file with added ".histogram" suffix.
 * Please note that the alignment is not 100% accurate and histogram might not reflect reality accurately. 
 * See statistics in log for better understanding of precision of the alignments.
 * @version 1.1
 * @date 2022-06-30
 * 
 * @copyright Copyright (c) 2022
 * This software is open source. There is no license to use this software. Do whatever you want with it.
 * [DISCLAIMER] code may not satisfy clean-code requirements. Scroll down on your own risk.
 */

#pragma once
#include "config.hpp"

namespace Alignment_parser
{

/**
 * @brief Function that handles the alignment on given Fastq/Fasta file
 * 
 * @param input_file path to the Fastq/Fasta file
 * @param grid Created/Loaded grid structure
 * @param mg Genes annotation
 */
void align_file(string input_file, 
                const Grid<int64_t> &grid, 
                const MergedGenes_s& mg);

void align_file(string input_file,
                string input_file_r,
                const Grid<int64_t> &grid,
                const MergedGenes_s& mg);

// End of Namespace
}
