// TODO: 
// - write cleaner and simpler command line arguments checker
// - finish writing file extensions checker
// - write file statistics function
// - restructure code to google C++ style
// - move functions to classes where suitable
// - reorder task exectution when valid files are given

#include <iostream>
#include <algorithm>
#include <set>
#include "config.h"
#include "options.hpp"
#include "../white_alignment/white_alignment.hpp"
#include "../white_minimizers/white_minimizers.hpp"
#include <bioparser/bioparser.hpp>
#include <map>


class FASTAformat {
	public:
		const char* name;
		const char* sequence;
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

bool comparator(std::pair<std::tuple<unsigned int, unsigned int, bool>, int> a, std::pair<std::tuple<unsigned int, unsigned int, bool>, int> b) {
	return a.second < b.second;
}

bool mapComparatorByPosition(std::pair<std::tuple<unsigned int, unsigned int, bool>, int> val1, std::pair<std::tuple<unsigned int, unsigned int, bool>, int> val2) {
    return std::get<1>(val1.first) < std::get<1>(val2.first);
}


int main(int argc, char* argv[]) {
	unsigned int k = 15, w = 5;
	double f = 0.001;

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
            std::vector<std::unique_ptr<FASTAformat>> reference_genome;
			std::vector<std::unique_ptr<FASTQformat>> fastq_objects;
			
			if (fasta_file_formats.count(first_file_extension)) {
				auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FASTAformat>(first_file);
				fasta_parser->parse(fasta_objects, -1);
                auto fasta_parser2 = bioparser::createParser<bioparser::FastaParser, FASTAformat>(second_file);
                fasta_parser2->parse(reference_genome, -1);
				std::vector<int> sequence_calc;
				long long sum = 0;
				for (int i = 0; i < fasta_objects.size(); i++) {
					sequence_calc.push_back(fasta_objects[i]->sequence_length);
					sum += fasta_objects[i]->sequence_length;
				}

				std::sort(sequence_calc.begin(), sequence_calc.end());
				
				outputFileStatistics(first_file, fasta_objects.size(), sum, sequence_calc[0], sequence_calc[sequence_calc.size() - 1]);
			
                std::string Cigar;
                unsigned int t_begin;
                int firstInd = rand() % fasta_objects.size();
                int secondInd = rand() % fasta_objects.size();
                const char* seq1 = fasta_objects[firstInd]->sequence;
                const char* seq2 = fasta_objects[secondInd]->sequence;
                int lenght1 = fasta_objects[firstInd]->sequence_length;
                int length2 = fasta_objects[secondInd]->sequence_length;
                int optimalAlign = white::PairwiseAlignment(seq1, lenght1, seq2, length2, white::AlignmentType::kGlobal, 2, -1, -2, Cigar, t_begin);
                std::cout << Cigar << std::endl;

				
				std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers;


                for (int i = 0; i < fasta_objects.size(); i++) {
					std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_of_sequence = white::minimizers(fasta_objects[i]->sequence,
																										(unsigned int)fasta_objects[i]->sequence_length,
																										k,
																										w);
					for (auto minimizer : minimizers_of_sequence) {
						minimizers.emplace_back(minimizer);
					}
				}

				std::set<std::tuple<unsigned int, unsigned int, bool>> distinct_minimizers;
				std::map<std::tuple<unsigned int, unsigned int, bool>, int> occurences_of_minimizers;

				for (auto minimizer : minimizers) {
					distinct_minimizers.insert(minimizer);
					if (occurences_of_minimizers.count(minimizer)) {
						occurences_of_minimizers[minimizer]++;
					} else {
						occurences_of_minimizers[minimizer] = 1;
					}
				}

				std::vector<std::pair<std::tuple<unsigned int, unsigned int, bool>, int>> occurences_of_minimizers_vector;

				for (auto it : occurences_of_minimizers) {
					occurences_of_minimizers_vector.push_back(it);
				}


				sort(occurences_of_minimizers_vector.begin(), occurences_of_minimizers_vector.end(), comparator);

				std::cout << "Number of distinct minimizers: " << distinct_minimizers.size() << std::endl;

				int singletons = 0;
				for (auto it : occurences_of_minimizers_vector) {
					if (it.second == 1) {
						singletons++;
					}
					else {
						break;
					}
				}

				std::cout << "Fraction of singletons: " << (double)singletons / minimizers.size() << std::endl;

				int skip = minimizers.size() * f;
				std::cout << skip << std::endl;
				std::cout << occurences_of_minimizers_vector.size() << std::endl;

				std::cout << occurences_of_minimizers_vector[1000].second << std::endl;
				std::cout << occurences_of_minimizers_vector[5000].second << std::endl;

				std::cout << "Number of occurences of the most frequent  minimizer when the top "
						 << f <<  " frequent minimizers are not taken in account: "
							<< occurences_of_minimizers_vector[occurences_of_minimizers_vector.size() - skip - 1].second << std::endl;

                std::vector<std::tuple<unsigned int, unsigned int, bool>> reference_minimizers = white::minimizers(reference_genome[1]->sequence,
                                                                                                                   (unsigned int)reference_genome[1]->sequence_length,
                                                                                                                   k,
                                                                                                                   w);
                std::set<std::tuple<unsigned int, unsigned int, bool>> distinct_reference_minimizers;
                std::map<std::tuple<unsigned int, unsigned int, bool>, int> occurences_of_reference_minimizers;
                std::vector<std::pair<std::tuple<unsigned int, unsigned int, bool>, int>> occurences_of_reference_minimizers_vector;

                for (auto eachMini : reference_minimizers) {
                    distinct_reference_minimizers.insert(eachMini);
                    if (occurences_of_reference_minimizers.count(eachMini)) {
                        occurences_of_reference_minimizers[eachMini]++;
                    } else {
                        occurences_of_reference_minimizers[eachMini] = 1;
                    }
                }
                for (auto oc : occurences_of_reference_minimizers) {
                    occurences_of_reference_minimizers_vector.push_back(oc);
                }

//                for (auto eachMin : distinct_reference_minimizers){
//                    std::cout << "min: " << std::get<0>(eachMin) << ", position: " << std::get<1>(eachMin) << std::endl;
//				}

                sort(occurences_of_reference_minimizers_vector.begin(), occurences_of_reference_minimizers_vector.end(), comparator);
                std::cout << "size: " << occurences_of_reference_minimizers_vector.size() << std::endl;

                int not_taken_in_account = distinct_reference_minimizers.size() * f;
                occurences_of_reference_minimizers_vector.erase(occurences_of_reference_minimizers_vector.begin(),
                                                                occurences_of_reference_minimizers_vector.begin() + not_taken_in_account);
                std::cout << "size: " << occurences_of_reference_minimizers_vector.size() << std::endl;

                //sort occurences of reference minimizers vector by position
                sort(occurences_of_reference_minimizers_vector.begin(), occurences_of_reference_minimizers_vector.end(), mapComparatorByPosition);


                for (int i = 0; i < fasta_objects.size(); i++) {
                    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_of_sequence = white::minimizers(fasta_objects[i]->sequence,
                                                                                                                         (unsigned int)fasta_objects[i]->sequence_length,
                                                                                                                         k,
                                                                                                                         w);
                    std::set<std::tuple<unsigned int, unsigned int, bool>> distinct_mins;
                    std::map<std::tuple<unsigned int, unsigned int, bool>, int> occurences_of_minimizers;

                    for (auto minim : minimizers_of_sequence) {
                       distinct_mins.insert(minim);
                        if (occurences_of_minimizers.count(minim)) {
                            occurences_of_minimizers[minim]++;
                        } else {
                            occurences_of_minimizers[minim] = 1;
                        }
                    }
                    std::vector<std::tuple <unsigned int, unsigned int>> matches;
                    for (auto min_ref : occurences_of_reference_minimizers_vector){
                        for (auto min : distinct_mins){
                            if (std::get<0>(min) == std::get<0>(min_ref.first)){
                                matches.push_back(std::make_tuple(std::get<1>(min), std::get<1>(min_ref.first)));
                            }
                        }
                    }
                    std::cout << fasta_objects[i]->name << "\t" <<
                    fasta_objects[i]->sequence_length << "\t" <<
                    "q_start" << "\t" <<
                    "q_end" << "\t" <<
                    "relative_strand" << "\t" <<
                    reference_genome[1]->name << "\t" <<
                    reference_genome[1]->sequence_length << "\t" <<
                    "t_start" << "\t" <<
                    "t_end" << "\t" <<
                    "number_of_matches" << "\t" <<
                    "alignment_length" << "\t" <<
                    "255" << std::endl;
                }

            }

			break;
	}

	return 0;
}

