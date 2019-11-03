// TODO: 
// - write cleaner and simpler command line arguments checker
// - finish writing file extensions checker
// - write file statistics function
// - restructure code to google C++ style
// - move functions to classes where suitable

#include <iostream>
#include <algorithm>
#include "config.h"
#include "options.h"
#include <bioparser/bioparser.hpp>


class FASTAformat {
	public:
		std::string name, sequence;
		int name_length, sequence_length;

		FASTAformat(
			   const char* name, std::uint32_t name_length,
			   const char* sequence, std::uint32_t sequence_length) {
			this->name = name;
			this->name_length = name_length;
			this->sequence = sequence;
			this->sequence_length = sequence_length;
		}
};

class FASTQformat {
	public:
		std::string name, sequence, quality;
		int name_length, sequence_length, quality_length;

		FASTQformat(
			   const char* name, std::uint32_t name_length,
			   const char* sequence, std::uint32_t sequence_length,
			   const char* quality, std::uint32_t quality_length) {
			this->name = name;
			this->name_length = name_length;
			this->sequence = sequence;
			this->sequence_length = sequence_length;
			this->quality = quality_length;
			this->quality_length = quality_length;
		}
};


int main(int argc, char* argv[]) {
	/*
	std::cout << "Test" << std::endl;
	
	std::cout << argc << std::endl;
	*/

	switch (argc) {
		case 1:
			noArgumentsPassed();
			break;

		case 2:
			if (std::string(argv[1]) == "-v" || std::string(argv[1]) == "--version") {
				printProgramVersion();
			} else if (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help") {
				printProgramHelp();
			} else {
				invalidArgumentsPassed();
			}

			break;

		case 3:
			std::string first_file = std::string(argv[1]), second_file = std::string(argv[2]);
			
			// test search file extension, will be changed
			int it = first_file.size() - 1;
			while (it--) {
				if (first_file[it] == '.') break;
			}
			std::string file_extension = first_file.substr(it + 1, first_file.size() - it);
			//std::cout << file_extension << std::endl;
			break;
	}

	std::vector<std::unique_ptr<FASTAformat>> fasta_objects;
	auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FASTAformat> (argv[1]);
	fasta_parser->parse(fasta_objects, -1);

	int mini = 999999;
	int maxi = 0;
	int sum = 0;
	for (int i = 0; i < fasta_objects.size(); i++) {
		mini = std::min(mini, fasta_objects[i]->sequence_length);
		maxi = std::max(maxi, fasta_objects[i]->sequence_length);
		sum += fasta_objects[i]->sequence_length;
	}
	
	double avg = (double)sum / fasta_objects.size();
	
	std::cout << "Number of sequences: " << fasta_objects.size() << std::endl;
	std::cout << "Average length of sequence: " << avg << std::endl;
	std::cout << "Minimal length of sequence: " << mini << std::endl;
	std::cout << "Maximal length of sequence: " << maxi << std::endl;

	return 0;
}
