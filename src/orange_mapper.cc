#include <unordered_map>
#include <functional>
#include <algorithm>
#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <iterator>
#include <valarray>
#include <cstdlib>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <set>

#include <bioparser/bioparser.hpp>
#include <thread_pool/thread_pool.hpp>

#include "orange_alignment.h"
#include "orange_minimizers.h"
#include "orange_mapper_conf.h"

namespace orange {
namespace mapper {

auto constexpr kStrandGapLim = static_cast<std::uint32_t>(10e3);
auto constexpr kGrupDiffLim = std::int64_t{500};

/**
 * @brief Supported genome file formats
 */
enum FileType { kFasta, kFastq, kUnsupported };

/**
 * Struct array representing long optins
 */
struct option const long_options[] = {
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'v'},
    {"type", required_argument, NULL, 't'},
    {"match", required_argument, 0, 'm'},
    {"mismatch", required_argument, 0, 's'},
    {"gap", required_argument, 0, 'g'},
    {"cigar", no_argument, NULL, 'v'},
    {"kmer", required_argument, NULL, 'k'},
    {"window", required_argument, NULL, 'w'},
    {"frequency", required_argument, NULL, 'f'},
    {"n_threads", required_argument, NULL, 'n'},
    {NULL, 0, NULL, 0}};

/**
 * @brief Holds FASTA/FASTQ sequecnes
 */
struct Sequence {
    /**
     * @brief Construct a new FastAQ object
     *
     * @param name sequence name
     * @param name_len sequence name length
     *
     * @param seq geneom sequence
     * @param seq_len genom sequence length
     */
    Sequence(char const* name, std::uint32_t name_len, char const* seq,
             std::uint32_t seq_len)
        : name_(name, name_len), seq_(seq, seq_len) {}

    std::string name_;
    std::string seq_;
};

/**
 * @brief Represents a portion of a vector
 *
 * @tparam T
 */
template <typename T>
class VectorSlice {
public:
    using vec_type = std::vector<T> const;
    using iterator = typename vec_type::const_iterator;

    VectorSlice(vec_type const& vec, std::size_t l, std::size_t r)
        : begin_{vec.cbegin()}, end_{vec.cbegin()} {
        std::advance(begin_, l);
        std::advance(end_, r);
    }

    VectorSlice(VectorSlice const& other, std::size_t l, std::size_t r)
        : begin_{other.begin()}, end_{other.begin()} {
        std::advance(begin_, l);
        std::advance(end_, r);
    }

    T& operator[](std::size_t pos) { return *(begin_ + pos); }

    T const& at(std::size_t pos) const { return *(begin_ + pos); }

    iterator begin() const { return begin_; }

    iterator end() const { return end_; }

    std::size_t size() const { return end_ - begin_; }

private:
    iterator begin_;
    iterator end_;
};

/**
 * @brief Represents matching of read K-Mer on reference K-Mer index
 *
 * @details <position on query, position on target>
 *
 */
using KMerMatch = std::tuple<std::uint32_t, std::uint32_t>;

using KMerMatches = std::vector<KMerMatch>;

using MatchSlice = VectorSlice<KMerMatch>;

/**
 * @brief FastA sequence alias.
 */
// Doable since FASTA doesn't contain quality information
using FastA = Sequence;

/**
 * @brief FASTQ genome sequences.
 *
 * @details Extends @ref orange::mapper::Sequence
 *      with string quiality information
 */
struct FastQ : public Sequence {
    FastQ(char const* name, std::uint32_t name_len, char const* seq,
          std::uint32_t seq_len, char const* quality_str,
          std::uint32_t quality_str_len)
        : Sequence(name, name_len, seq, seq_len),
          quality_str_(quality_str, quality_str_len) {}

    std::string quality_str_;
};

using MinimizerIndex =
    std::unordered_map<minimizers::KMerVal, minimizers::KMerLocs>;

using SequencePtr = std::unique_ptr<Sequence>;

using FastQPtr = std::unique_ptr<FastQ>;

/**
 * @brief Vector of pointers to parsed @ref orange::mapper::Sequence
 */
using VecSeqPtr = std::vector<SequencePtr>;

/**
 * @brief @ref orange::mapper::VectorSlice of @ref orange::mapper::VecSeqPtr
 */
using SliceSeqPtr = VectorSlice<SequencePtr>;

/**
 * @brief Vector of pointers to parsed @ref orange::mapper::FastQ
 */
using VecFastQPtr = std::vector<FastQPtr>;

/**
 * @brief @ref orange::mapper::VectorSlice of @ref orange::mapper::VecFastQPtr
 *
 */
using SliceFastQPtr = VectorSlice<FastQPtr>;

/**
 * @brief Simple intarface for creating @ref bioparser FASTA parser
 */
auto createFastaParser =
    std::function<std::unique_ptr<bioparser::Parser<FastA>>(
        std::string const& path)>{
        &bioparser::createParser<bioparser::FastaParser, FastA>};

/**
 * @brief Simple interface for creating @ref bioparser FASTQ parser
 */
auto createFastqParser =
    std::function<std::unique_ptr<bioparser::Parser<FastQ>>(
        std::string const& path)>{
        &bioparser::createParser<bioparser::FastqParser, FastQ>};

/**
 * @brief Prints program version to stderr,
 *      exits upon completion
 */
auto printVersion() {
    std::cerr << orange_mapper_VERSION_MAJOR << '.'
              << orange_mapper_VERSION_MINOR << '.'
              << orange_mapper_VERSION_PATCH << '\n';
    std::exit(EXIT_SUCCESS);
}

/**
 * @brief Prints program usage information to stderr,
 *      exits upon completion
 */
auto printHelp() {
    std::cerr
        << "Genome sequence mapper.\n"
        << "Reads are supported in FASTA and FASTQ formats while "
        << "reference is expected to be FASTA.\n\n"
        << "Usage:\n\t orange_mapper <reads> <reference> "
        << "-t <alignment_type> -m <match> -s <mismatch> -g <gap>\n\n"
        << "Options:\n\t-h\thelp\n\t-v\tverions\n\n\t"
           "-t\talignment type\n\t-m\tmatch score\n\t"
           "-s\tmismatch score\n\t-g\tgap score\n\n\t"
           "-k\tk-mer size\n\t-w\twindow length\n\t"
           "-f\tignore top frequent minimizers\n"
           "-\tinclude CIGAR string in PAF output\n"
           "-tn\tspecify the number of threads for mapping (defaults to one)\n";
    std::exit(EXIT_SUCCESS);
}

/**
 * @brief parses alignment type from comman line argument
 *
 * @param type command line argument value
 * @return auto @ref orange::alignment::AlignmentType
 */
auto parseAlignType(char const* type) {
    auto str_val = std::string_view{type};

    if (str_val == "global")
        return alignment::AlignmentType::kGlobal;
    if (str_val == "local")
        return alignment::AlignmentType::kLocal;
    if (str_val == "semi-global")
        return alignment::AlignmentType::kSemiGlobal;

    throw std::invalid_argument(
        "Unknow alignment type.\n"
        "Use orange_mapper -h for further information");
}

/**
 * @brief Parses command line options.
 *
 * @details After option parsing,
 *      all non option elements are placed at the end of the array.
 *
 *      </br>
 *      Associated actions with each option are executed by
 *      this function considering the option was set at launch.
 *
 * @param argc number of command line arguments
 * @param argv commane line arguments
 *
 * @param conf sequence alignemnt config passed over command line options
 *
 * @return int index of first non-option element in argv
 */
auto parseOptions(int argc, char* argv[], alignment::AlignConf& a_conf,
                  minimizers::MinimizerConf& m_conf, std::uint8_t& n_threads) {
    auto opt = int{};
    while ((opt = getopt_long(argc, argv, "hvm:s:g:t:n:c", long_options,
                              NULL)) != -1) {
        switch (opt) {
            case 'v':
                printVersion();
                break;
            case 'h':
                printHelp();
                break;
            case 't':
                a_conf.type_ = parseAlignType(optarg);
                break;
            case 'm':
                a_conf.match_ = std::atoi(optarg);
                break;
            case 's':
                a_conf.mismatch_ = std::atoi(optarg);
                break;
            case 'g':
                a_conf.gap_ = std::atoi(optarg);
                break;
            case 'c':
                a_conf.cigar_ = true;
                break;
            case 'k':
                m_conf.k_ = std::atoi(optarg);
                break;
            case 'w':
                m_conf.win_len_ = std::atoi(optarg);
                break;
            case 'f':
                m_conf.f_ = std::atoi(optarg);
                break;
            case 'n':
                n_threads = std::atoi(optarg);
                break;
            default:
                throw std::invalid_argument(
                    "Uknown option.\n\tUse: 'orange_mapper -h'");
                break;
        }
    }

    return optind;
}

/**
 * @brief Get the @ref orange::mapper::FileType
 *
 * @param file path to file
 * @return corresponding @ref orange::mapper::FileType
 */
auto parseFileType(std::string_view file) {
    auto format_ex = std::invalid_argument(
        "Unsuppored file format/extension!\n"
        "\t Supported: FASTA (reads and reference), FASTQ "
        "(reads)\n\nUse -h for help\n");

    // Find furthest extension
    auto dot = file.rfind('.');
    if (dot == std::string_view::npos)
        throw std::invalid_argument("Missing file extension.");
    else if (file.size() - dot < 2)
        throw format_ex;
    else if (file[dot + 1] == 'g' && file[dot + 2] == 'z') {
        // Ignoring .gz
        file.remove_suffix(file.size() - dot);
        dot = file.rfind('.');
    }

    file.remove_prefix(dot);
    if (file == ".fasta" || file == ".fa")
        return FileType::kFasta;
    else if (file == ".fastq" || file == ".fq")
        return FileType::kFastq;

    throw format_ex;
}

/**
 * @brief Loads FASTA/FASTQ files into
 *
 * @param path_to_file path to file contaning FASTA/FASTQ sequences
 *
 * @return @ref orange::mapper::VecSeqPtr
 */
auto loadFile(std::string const& path_to_file, FileType const& type) {
    auto parse = [](auto parser, auto& objects) { parser->parse(objects, -1); };

    if (type == FileType::kFasta) {
        auto objects = VecSeqPtr{};
        parse(createFastaParser(path_to_file), objects);

        return objects;
    } else {
        auto objects = VecFastQPtr{};
        parse(createFastqParser(path_to_file), objects);

        return VecSeqPtr{std::make_move_iterator(objects.begin()),
                         std::make_move_iterator(objects.end())};
    }
}

/**
 * @brief Prints stats for loaded sequences
 *
 * @details Calculated stats: <br\>
 *          <ol>
 *              <li>Number of sequences</li>
 *              <li>Average sequence length</li>
 *              <li>Minimum sequence length</li>
 *              <li>Maximum sequence length</li>
 *          <\ol>
 *
 * @param vec_seqs reference to loaded set of sequences
 */
auto printStats(std::string_view const& origin, VecSeqPtr const& vec_seq) {
    auto n_seq = vec_seq.size();

    auto avg_len = double{};
    auto min_len = std::uint32_t{std::numeric_limits<std::uint32_t>::max()};
    auto max_len = std::uint32_t{std::numeric_limits<std::uint32_t>::min()};

    std::for_each(begin(vec_seq), end(vec_seq),
                  [&avg_len, &min_len, &max_len](auto& seq_ptr) {
                      auto const& seq_len =
                          static_cast<std::uint32_t>(seq_ptr->seq_.size());

                      avg_len += seq_len;
                      min_len = std::min(min_len, seq_len);
                      max_len = std::max(max_len, seq_len);
                  });

    avg_len /= n_seq;

    std::cerr << '\n'
              << origin << " Stats:\n"
              << "\tNumber or sequences: " << n_seq << '\n'
              << "\tAverage length: " << avg_len << '\n'
              << "\tMinimum length: " << min_len << '\n'
              << "\tMaximum length: " << max_len << '\n';
}

auto printRngAlign(Sequence const& query, Sequence const& target,
                   alignment::AlignConf conf) {
    auto cigar = std::string{};
    auto target_begin = std::uint32_t{};

    auto const align_score = alignment::pairwiseAlignment(
        query.seq_.c_str(), query.seq_.size(), target.seq_.c_str(),
        target.seq_.size(), conf.type_, conf.match_, conf.mismatch_, conf.gap_,
        cigar, target_begin);

    std::cerr << "Random alignment score: " << align_score << '\n'
              << "CIGAR:\n\t" << cigar << "\n\n";
}

auto printMinimizerStats(VecSeqPtr const& reads,
                         minimizers::MinimizerConf const& conf) {
    using MinimzMap = std::unordered_map<minimizers::KMerVal, std::uint64_t>;
    using MinimizerIter = MinimzMap::iterator;

    auto minimz = std::unordered_map<minimizers::KMerVal, std::uint64_t>{};

    for (auto const& it : reads) {
        auto& read = *it.get();
        for (auto const& minim : minimizers::minimizers(
                 read.seq_.c_str(), read.seq_.size(), conf.k_, conf.win_len_)) {
            ++minimz[std::get<0>(minim)];
        }
    }

    std::cerr << "Number of distincs minimizers for reads: " << minimz.size()
              << '\n';

    auto cmp = [](MinimizerIter const& l, MinimizerIter const& r) {
        return l->second > r->second;
    };

    auto n_singletons = std::uint64_t{0};
    auto ignore_set = std::multiset<MinimizerIter, decltype(cmp)>(cmp);
    auto ignore_cnt = static_cast<std::size_t>(minimz.size() * conf.f_ + 1);

    for (auto iter = minimz.begin(); iter != minimz.end(); ++iter) {
        ignore_set.insert(iter);
        if (ignore_set.size() > ignore_cnt)
            ignore_set.erase(std::prev(ignore_set.end()));

        if (iter->second == 1)
            ++n_singletons;
    }

    auto most_freq = (*ignore_set.rbegin())->second;

    ignore_set.erase(std::prev(ignore_set.end()));
    for (auto const& it : ignore_set) minimz.erase(it);
    ignore_set.clear();

    std::cerr << "Fraction: " << 1.0 * n_singletons / most_freq << '\n';
}

/**
 * @brief Creates a minimizer index for the reference genome
 */
MinimizerIndex createRefMinimzIndex(std::string const& ref,
                                    minimizers::MinimizerConf const& conf) {
    auto ref_index = MinimizerIndex{};
    using MinimizerIter = MinimizerIndex::iterator;

    for (auto const& [kmer, pos, org] : minimizers::minimizers(
             ref.c_str(), ref.size(), conf.k_, conf.win_len_)) {
        ref_index[kmer].emplace_back(pos, org);
    }

    auto cmp = [](MinimizerIter const& l, MinimizerIter const& r) {
        return l->second.size() > r->second.size();
    };

    auto ignore_cnt = ref_index.size() * conf.f_;
    auto ignore_set = std::multiset<MinimizerIter, decltype(cmp)>(cmp);

    for (auto iter = ref_index.begin(); iter != ref_index.end(); ++iter) {
        ignore_set.insert(iter);
        if (ignore_set.size() > ignore_cnt)
            ignore_set.erase(std::prev(ignore_set.end()));
    }

    /* clang-format off */
    for (auto const& it : ignore_set) 
        ref_index.erase(it);
    ignore_set.clear();
    /* clang-format on */

    return ref_index;
}

/**
 * @brief
 *
 * @details LIS is done on the reference position, 3rd argument of tuple
 *
 * @param slice
 * @return auto
 */
auto LISAlgo(MatchSlice const& matches) {
    /* Dynamic LIS building */
    std::vector<std::size_t> lis{0u};  //<<< Current max lis sequence
    std::vector<std::size_t> org{
        0u};  //<<< for each element, remembers the positon of seuqence start

    /**
     * @brief Lamba expression used for binary search in LIS
     */
    auto cmp_lambda = [&matches](std::size_t const& a, std::size_t const& b) {
        return std::get<1>(matches.at(a)) < std::get<1>(matches.at(b));
    };

    auto b_search = [&lis, &cmp_lambda](std::size_t const pos) {
        return std::upper_bound(lis.begin(), lis.end(), pos, cmp_lambda);
    };

    for (std::size_t i = 1; i < matches.size(); ++i) {
        auto res = b_search(i);
        if (res == lis.end()) {
            org.push_back(org.back());
            lis.push_back(i);
        } else {
            if (res == lis.begin())
                org[0] = i;
            else {
                auto pos = std::distance(lis.begin(), res);
                org[pos] = org[pos - 1];
            }

            *res = i;
        }
    }

    return std::make_tuple(std::make_pair(std::get<0>(matches.at(org.back())),
                                          std::get<0>(matches.at(lis.back()))),

                           std::make_pair(std::get<1>(matches.at(org.back())),
                                          std::get<1>(matches.at(lis.back()))),
                           lis.size());
}

/**
 * @brief Returns absolute difference of two undinged 32 bit ints
 */
std::uint32_t unsigned_dif(std::uint32_t const& a, std::uint32_t const& b) {
    return a > b ? a - b : b - a;
};

/**
 * @brief signed difference between query and target position
 */
std::int64_t qt_dif(KMerMatch const& a) {
    return static_cast<std::int64_t>(std::get<0>(a)) - std::get<1>(a);
};

auto generatePAF(SequencePtr const& query, SequencePtr const& target,
                 MatchSlice const& matches, char const rel_strand,
                 alignment::AlignConf const& a_conf,
                 minimizers::MinimizerConf const& m_conf) {
    auto [query_se, target_se, n_matches] = [&matches]() {
        using StartEnd = std::pair<std::uint32_t, std::uint32_t>;
        using LISResult = std::tuple<StartEnd, StartEnd, std::size_t>;

        auto const is_split = [&matches](std::size_t const& curr,
                                         std::size_t const& prev) -> bool {
            return unsigned_dif(std::get<0>(matches.at(curr)),
                                std::get<0>(matches.at(prev))) > kStrandGapLim;
        };

        auto prev = 0;
        auto ret = LISResult{{0, 0}, {0, 0}, std::size_t{1}};
        for (std::size_t curr{0}; curr < matches.size(); ++curr) {
            if (curr + 1 == matches.size() || is_split(curr, prev)) {
                auto res = LISAlgo(std::move(MatchSlice(matches, prev, curr)));
                if (std::get<2>(res) > std::get<2>(ret))
                    ret = std::move(res);

                prev = curr;
            }
        }

        return ret;
    }();

    // Ignore insignificat
    if (n_matches <= 4)
        return;

    std::stringstream paf_ss;
    auto gen_seq_data = [&paf_ss](SequencePtr const& seq,
                                  auto start_end) -> void {
        paf_ss << seq.get()->name_ << '\t';
        paf_ss << seq.get()->seq_.size() << '\t';
        paf_ss << start_end.first << '\t';
        paf_ss << start_end.second << '\t';
    };

    if (rel_strand == '-')
        std::swap(query_se.first, query_se.second);

    gen_seq_data(query, query_se);
    paf_ss << rel_strand << '\t';
    gen_seq_data(target, target_se);

    auto const& [query_start, query_end] = query_se;
    auto const& [target_start, target_end] = target_se;

    auto const& que = *query.get();
    auto const& tar = *target.get();

    if (!a_conf.cigar_) {
        paf_ss << n_matches << '\t';
        paf_ss << target_end - target_start << '\t';
        paf_ss << "255\t\n";
    } else {
        std::string cigar{};
        std::uint32_t target_begin{0};

        alignment::pairwiseAlignment(
            que.seq_.c_str() + query_start, query_end - query_start + m_conf.k_,
            tar.seq_.c_str() + target_start,
            target_end - target_start + m_conf.k_, a_conf.type_, a_conf.match_,
            a_conf.mismatch_, a_conf.gap_, cigar, target_begin);

        auto nxt_num_in_cigar = [&cigar](auto& iter) {
            std::string buff{};
            while (iter < cigar.cend() && std::isdigit(*iter))
                buff += *(iter++);
            return std::stoi(buff);
        };

        std::uint32_t cigar_matches{0};
        std::uint32_t alignment_len{0};

        auto iter = cigar.cbegin();
        while (iter < cigar.cend()) {
            auto c_num = nxt_num_in_cigar(iter);
            auto curr_char = *(iter++);
            if (curr_char == 'I' || curr_char == 'M') {
                alignment_len += c_num;
                if (curr_char == 'M')
                    cigar_matches += c_num;
            } else
                alignment_len -= c_num;
        }

        paf_ss << cigar_matches << '\t';
        paf_ss << alignment_len << '\t';
        paf_ss << "\t\n";

        paf_ss << cigar << '\n';
    }

    std::cout << paf_ss.str();
}

std::unordered_map<char,
                   std::function<bool(KMerMatch const&, KMerMatch const&)>>
    sort_dict{
        // Reminder: KMerMatch(query_pos, target_pos)
        {'+',
         [](KMerMatch const& a, KMerMatch const& b) -> bool { return a < b; }},

        {'-',
         [](KMerMatch const& a, KMerMatch const& b) -> bool {
             if (std::get<0>(a) != std::get<0>(b))
                 return std::get<0>(a) > std::get<0>(b);
             return std::get<1>(a) < std::get<1>(b);
         }},

        // Sort diagonals
        {'d',
         [](KMerMatch const& a, KMerMatch const& b) -> bool {
             return std::get<0>(a) - std::get<1>(a) <
                    std::get<0>(b) - std::get<1>(b);
         }}

    };

void mapReads(SliceSeqPtr reads, SequencePtr const& target,
              MinimizerIndex const& target_index,
              alignment::AlignConf const& a_conf,
              minimizers::MinimizerConf const& m_conf) {
    for (auto const& read : reads) {
        // One read -> One PAF mapping
        std::array<KMerMatches, 2>
            matches;  // 1 -> on the same strand, 0 -> differ

        for (auto const& kmer : minimizers::minimizers(
                 read.get()->seq_.c_str(), read.get()->seq_.size(), m_conf.k_,
                 m_conf.win_len_)) {
            auto [kmer_val, kmer_pos, kmer_org] = kmer;
            auto ref_kmer = target_index.find(kmer_val);

            if (ref_kmer != target_index.end()) {
                auto ref_data = ref_kmer->second;
                for (auto [ref_pos, ref_org] : ref_data)
                    matches[kmer_org == ref_org].emplace_back(kmer_pos,
                                                              ref_pos);
            }
        }

        auto const genPAF = [&read, &target, &a_conf, &m_conf](
                                MatchSlice const& matches,
                                char const rel_strand) -> void {
            generatePAF(read, target, matches, rel_strand, a_conf, m_conf);
        };

        for (auto& match_group : matches)
            std::sort(match_group.begin(), match_group.end(), sort_dict['d']);

        auto const group_split = [](KMerMatches const& matches,
                                    std::size_t const& l,
                                    std::size_t const& r) {
            // takes the advantege of ascending sorting
            return qt_dif(matches[r]) - qt_dif(matches[l]) > kGrupDiffLim;
        };

        auto const group_matches = [&genPAF, &group_split](
                                       KMerMatches& matches,
                                       char const rel_strand) {
            std::size_t group_start{0};
            for (std::size_t group_end{0}; group_end <= matches.size();
                 ++group_end) {
                if (group_end == matches.size() ||
                    group_split(matches, group_start, group_end)) {
                    std::sort(matches.begin() + group_start,
                              matches.begin() + group_end,
                              sort_dict[rel_strand]);

                    genPAF(MatchSlice(matches, group_start, group_end),
                           rel_strand);
                    group_start = group_end;
                }
            }
        };

        if (!matches[1].empty())
            group_matches(matches[1], '+');
        if (!matches[0].empty())
            group_matches(matches[1], '-');
    }
}

/**
 * @brief Splits mapping procedure in thread blocks
 *
 * @return auto
 */
void threadMapping(VecSeqPtr const& reads, SequencePtr const& reference,
                   minimizers::MinimizerConf const& m_conf,
                   alignment::AlignConf const& a_conf,
                   std::uint8_t const n_threads) {
    auto thread_pool = thread_pool::createThreadPool(n_threads);
    auto futures = std::vector<std::future<void>>{};

    futures.reserve(std::ceil(reads.size() / n_threads));

    auto jmp = static_cast<std::size_t>(reads.size() / n_threads);
    auto ref_index = createRefMinimzIndex(reference.get()->seq_, m_conf);

    std::size_t block_start{0};
    for (; block_start + jmp < reads.size(); block_start += jmp) {
        auto slice = SliceSeqPtr(reads, block_start, block_start + jmp);
        futures.push_back(thread_pool->submit(
            mapReads, slice, std::ref(reference), std::ref(ref_index),
            std::ref(a_conf), std::ref(m_conf)));
    }

    // Last slice
    auto const slice = SliceSeqPtr(reads, block_start, reads.size());
    futures.push_back(thread_pool->submit(mapReads, slice, std::ref(reference),
                                          std::ref(ref_index), std::ref(a_conf),
                                          std::ref(m_conf)));

    std::for_each(futures.begin(), futures.end(),
                  [](auto const& f) { f.wait(); });
}

}  // namespace mapper
}  // namespace orange

/**
 * @brief Main program entry point
 *
 * @param argc number of command line arguments
 * @param argv command line arguments (options and/or sources)
 *
 * @return int execution success status
 */
int main(int argc, char* argv[]) {
    using namespace orange;

    try {
        auto n_threads = std::uint8_t{1};
        // Alignment configuration
        auto a_conf = alignment::AlignConf{};
        // Minimizer configuration
        auto m_conf = minimizers::MinimizerConf{};

        auto arg_index =
            mapper::parseOptions(argc, argv, a_conf, m_conf, n_threads);

        // Check number of passed arguments
        if (argc - arg_index < 2) {
            std::cerr << "Invalid number of arguments.\n"
                         "\tRequired two input files:\n"
                         "\t\t1. Set of fragments (FASTA/FASTQ)\n"
                         "\t\t2. Reference genome (FASTA)\n\n"
                         "Followed up with alignment algorithm specification "
                         "{global, local, semi-global}\n"
                         "and match, mismatch, gap scoring values for the "
                         "coresponding type.\n"
                         "\tEg. orange_mapper escherichia_coli_r7_reads.fastq"
                         " escherichia_coli_reference.fasta -t global -m 1 -s "
                         "-1 -g -1\n";
            return EXIT_FAILURE;
        }

        // Extracting resource paths
        auto path_to_reads = std::string{argv[arg_index]};
        auto path_to_ref = std::string{argv[arg_index + 1]};

        // Determine scan formats
        auto reads_type = mapper::parseFileType(path_to_reads);
        auto refs_type = mapper::parseFileType(path_to_ref);

        // Report start of file loading
        // std::cout << "Started loading files\n";

        // Load and parse data from files
        auto const reads = mapper::loadFile(path_to_reads, reads_type);
        auto const refs = mapper::loadFile(path_to_ref, refs_type);

        // Report end of file loading
        // std::cout << "Finsihed loading files\n\n";

        // For fast I/O
        std::ios_base::sync_with_stdio(false);
        std::cin.tie(nullptr);

        for (auto const& ref : refs)
            mapper::threadMapping(reads, ref, m_conf, a_conf, n_threads);

    } catch (std::exception const& e) {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
