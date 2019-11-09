#include <iostream>

void noArgumentsPassed() {
	std::cout << "No arguments passed to program." << std::endl;
}

void printProgramVersion() {
	std::cout << "v" << WHITE_MAPPER_VERSION_MAJOR << "." << WHITE_MAPPER_VERSION_MINOR << std::endl;
}

void printProgramHelp() {
	std::cout << "The program takes 2 arguments:" << std::endl << std::endl
		  << "1st argument should be a file in FASTA or FASTQ format" << std::endl
		  << "2nd argument should be a file in FASTA format" << std::endl << std::endl
		  << "-v (or --version) flag for version message" << std::endl
		  << "-h (or --help) for displaying help message" << std::endl << std::endl
		  << "This help message will be changed and updated later." << std::endl;
}

void invalidArgumentsPassed() {
	std::cout << "Invalid arguments passed to program" << std::endl;
}
