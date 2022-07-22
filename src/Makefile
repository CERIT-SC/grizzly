CC=g++
out_file=grid_mapping
in_files=alignment_parser.cpp utils.cpp config.cpp inferencer.cpp kmer_converter.cpp testing.cpp extend_grid.cpp grid_mapping.cpp
flags=-Wall -Wextra -pedantic -lboost_system -lboost_iostreams -lboost_program_options -O3 -std=c++17 -fopenmp 


all:
	@$(CC) $(in_files) -o $(out_file) $(flags)
