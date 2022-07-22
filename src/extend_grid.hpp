/**
 * @file extend_grid.hpp
 * @author Jan Šamánek (jansamanek@email.cz)
 * @brief Component of the Grizzly software.
 * Contains function for inserting modification to the grid.
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
 * @brief Function which parses the input file and applies changes in the transcripts to the grid.
 * 
 * @param grid Grid structure
 * @param ga Genes annotation
 * @param database_fd File descriptor of the input file.
 * Input file needs to have structure
 * Transcript_name  relative_index   old_sequenc    new_seq1,new_seq2,...
 */
void extend_grid(Grid<int64_t> &grid, 
                 const GenesAnnotation_s& ga,
                 ifstream &database_fd);