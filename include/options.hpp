#include <iostream>

std::set<std::string> valid_file_formats = {"fasta", "fa", "fastq", "fq", "fasta.gz", "fa.gz", "fastq.gz", "fq.gz", "fna.gz"};
std::set<std::string> fasta_file_formats = {"fasta", "fa", "fasta.gz", "fa.gz", "fna.gz"};
std::set<std::string> fastq_file_formats = {"fastq", "fq", "fastq.gz", "fq.gz"};

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

bool checkFileFormats(std::string first_file, std::string second_file,
						std::string& first_file_extension, std::string& second_file_extension) {
	auto first_file_dot_pos = first_file.find('.');
	auto second_file_dot_pos = second_file.find('.');

	if (first_file_dot_pos != std::string::npos && second_file_dot_pos != std::string::npos) {
		first_file_extension = first_file.substr(first_file_dot_pos + 1);
		second_file_extension = second_file.substr(second_file_dot_pos + 1);

		// std::cout << first_file_extension << " " << second_file_extension << std::endl;
 		if (valid_file_formats.count(first_file_extension)
		 	&& fasta_file_formats.count(second_file_extension)) {
				std::cerr << "All good!" << std::endl;
				return true;
		} else {
			std::cerr << "Invalid file formats passed to program!" << std::endl;
			return false;
		}
	}
}

void outputFileStatistics(std::string filename, int size, long long sum, int first_element, int last_element) {
	std::cerr << "File " << "\"" << filename << "\"" << " statistics:" << std::endl;  
				std::cerr << "Number of sequences: " << size << std::endl;
				std::cerr << "Average length of sequence: " << sum / size << std::endl;
				std::cerr << "Minimal length of sequence: " << first_element << std::endl;
				std::cerr << "Maximal length of sequence: " << last_element << std::endl;
}