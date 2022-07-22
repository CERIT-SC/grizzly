/**
 * @file testing.hpp
 * @author Jan Šamánek (jansamanek@email.cz)
 * @brief Component of the Grizzly software.
 * Containes function for testing grid alignment performance.
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
                                    const size_t max_read_size);


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
                                            const size_t max_read_size);
