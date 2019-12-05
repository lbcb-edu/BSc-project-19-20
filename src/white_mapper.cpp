// TODO: 
// - write cleaner and simpler command line arguments checker
// - finish writing file extensions checker
// - write file statistics function
// - restructure code to google C++ style
// - move functions to classes where suitable

#include <iostream>
#include <algorithm>
#include <set>
#include "config.h"
#include "options.hpp"
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
	
			std::string first_file_extension, second_file_extension;
			if (!checkFileFormats(first_file, second_file, first_file_extension, second_file_extension)) {
				return 0;
			}
			
			std::vector<std::unique_ptr<FASTAformat>> fasta_objects;
			std::vector<std::unique_ptr<FASTQformat>> fastq_objects;
			
			if (fasta_file_formats.count(first_file_extension)) {
				auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FASTAformat>(first_file);
				fasta_parser->parse(fasta_objects, -1);
				std::vector<int> sequence_calc;
				long long sum = 0;
				for (int i = 0; i < fasta_objects.size(); i++) {
					sequence_calc.push_back(fasta_objects[i]->sequence_length);
					sum += fasta_objects[i]->sequence_length;
				}
				std::sort(sequence_calc.begin(), sequence_calc.end());
				
				outputFileStatistics(first_file, fasta_objects.size(), sum, sequence_calc[0], sequence_calc[sequence_calc.size() - 1]);
			}
			
			break;
	}

	return 0;
}

