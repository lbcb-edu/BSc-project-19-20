#include "white_minimizers.hpp"
#include <vector>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <map>


namespace white {

	std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(const char* sequence, unsigned int sequence_length,
																		unsigned int k, unsigned int window_length) {

	std::vector<std::tuple<unsigned int, unsigned int, bool>> ret;

	unsigned int current_kmer = 0;
	unsigned int rev_current_kmer = 0;
	unsigned int next_kmer = 0;
	unsigned int rev_next_kmer = 0;
	unsigned int mask = 0;

	for (int i = 0; i < k; i++) {
		mask = (mask << 2) + 3;
	}

	for (int i = 0; i < k; i++) {
		switch (sequence[i]) {
			case 'A':
				current_kmer = (current_kmer << 2) & mask;
				break;
			case 'C':
				current_kmer = ((current_kmer << 2) | 1) & mask;
				break;
			case 'G':
				current_kmer = ((current_kmer << 2) | 2) & mask;
				break;
			case 'T':
				current_kmer = ((current_kmer << 2) | 3) & mask;
				break;
		}

		switch (sequence[sequence_length - 1 - i]) {
			case 'A':
				rev_current_kmer = ((rev_current_kmer << 2) | 3)& mask;
				break;
			case 'C':
				rev_current_kmer = ((rev_current_kmer << 2) | 2) & mask;
				break;
			case 'G':
				rev_current_kmer = ((rev_current_kmer << 2) | 1) & mask;
				break;
			case 'T':
				rev_current_kmer = (rev_current_kmer << 2) & mask;
				break;

		}
	}

	unsigned int last_minimizer = 0;
	long long int last_pos = -1;
	bool origin = true;

	for (long long int i = 0; i < sequence_length - (window_length - 1 + k) + 1; i++) { 
		switch (sequence[i + k]) {
			case 'A':
				next_kmer = (current_kmer << 2) & mask;
				break;
			case 'C':
				next_kmer = ((current_kmer << 2) | 1) & mask;
				break;
			case 'G':
				next_kmer = ((current_kmer << 2) | 2) & mask;
				break;
			case 'T':
				next_kmer = ((current_kmer << 2) | 3) & mask;
				break;
		}

		switch (sequence[sequence_length - 1 - (i + k)]) {
			case 'A':
				rev_next_kmer = ((rev_current_kmer << 2) | 3) & mask;
				break;
			case 'C':
				rev_next_kmer = ((rev_current_kmer << 2) | 2) & mask;
				break;
			case 'G':
				rev_next_kmer = ((rev_current_kmer << 2) | 1) & mask;
				break;
			case 'T':
				rev_next_kmer = (rev_current_kmer << 2) & mask;
				break;
		}

		if ((origin && last_pos < i) || (!origin && last_pos > (sequence_length - 1) - i - (k - 1))) {
		
			if (current_kmer <= rev_current_kmer) {
				last_minimizer = current_kmer;
				last_pos = i;
				origin = true;
			}
			else {
				last_minimizer = rev_current_kmer;
				last_pos = (sequence_length - 1) - i - (k - 1);
				origin = false;
			}

			unsigned int end = i + k + window_length - 1;

			for (unsigned int j = i + k; j < end; j++) {
				switch (sequence[j]) {
					case 'A':
						current_kmer = (current_kmer << 2) & mask;
						break;
					case 'C':
						current_kmer = ((current_kmer << 2) | 1) & mask;
						break;
					case 'G':
						current_kmer = ((current_kmer << 2) | 2) & mask;
						break;
					case 'T':
						current_kmer = ((current_kmer << 2) | 3) & mask;
						break;
				}

				switch (sequence[sequence_length - 1 - j]) {
					case 'A':
						rev_current_kmer = ((rev_current_kmer << 2) | 3) & mask;
						break;
					case 'C':
						rev_current_kmer = ((rev_current_kmer << 2) | 2) & mask;
						break;
					case 'G':
						rev_current_kmer = ((rev_current_kmer << 2) | 1) & mask;
						break;
					case 'T':
						rev_current_kmer = (rev_current_kmer << 2) & mask;
						break;
				}

				if (last_minimizer >= current_kmer) {
					last_minimizer = current_kmer;	
					last_pos = j - (k - 1);
					origin = true;
				}
				if (last_minimizer >= rev_current_kmer) {
					last_minimizer = rev_current_kmer;
					last_pos = (sequence_length - 1) - j;
					origin = false;
				}
			}

		}
		
		else {
			current_kmer = 0;
			rev_current_kmer = 0;

			unsigned int start = i + window_length - 1; 
			unsigned int end = i + k + window_length - 1; 

			for (unsigned int j = start; j < end; j++) {
				switch (sequence[j]) {
					case 'A':
						current_kmer = (current_kmer << 2) & mask;
						break;
					case 'C':
						current_kmer = ((current_kmer << 2) | 1) & mask;
						break;
					case 'G':
						current_kmer = ((current_kmer << 2) | 2) & mask;
						break;
					case 'T':
						current_kmer = ((current_kmer << 2) | 3) & mask;
						break;

				}

				switch (sequence[sequence_length - 1 - j]) {
					case 'A':
						rev_current_kmer = ((rev_current_kmer << 2) | 3) & mask;
						break;
					case 'C':
						rev_current_kmer = ((rev_current_kmer << 2) | 2) & mask;
						break;
					case 'G':
						rev_current_kmer = ((rev_current_kmer << 2) | 1) & mask;
						break;
					case 'T':
						rev_current_kmer = (rev_current_kmer << 2) & mask;
						break;
				}
			}

			if (last_minimizer >= current_kmer) {
				last_minimizer = current_kmer;
				last_pos = start;
				origin = true;
			}
			if (last_minimizer >= rev_current_kmer) {
				last_minimizer = rev_current_kmer;
				last_pos = (sequence_length - 1) - start - (k - 1);
				origin = false;
			}
		}

		current_kmer = next_kmer;
		rev_current_kmer = rev_next_kmer;

		std::tuple<unsigned int, unsigned int, bool> current_minimizer = std::make_tuple(last_minimizer, (unsigned int) last_pos, origin);

		if (ret.empty())
			ret.emplace_back(current_minimizer);

		else if (ret[ret.size() - 1] != current_minimizer)
			ret.emplace_back(current_minimizer);

	}

	return ret;
	}

}

