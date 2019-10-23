#pragma once

#include <bioparser/bioparser.hpp>
#include <memory>
#include <type_traits>
#include <vector>

#include "sequence.hpp"

namespace mapper {

template <typename T>
using Sequences = ::std::vector<::std::unique_ptr<T>>;

template <typename DataFormat>
::std::vector<::std::unique_ptr<Sequence>> parse(
    const ::std::string& input_file) {
  constexpr bool is_fasta = std::is_same_v<DataFormat, Sequence>;

  Sequences<DataFormat> sequences;
  ::std::unique_ptr<::bioparser::Parser<DataFormat>> parser;

  if constexpr (is_fasta) {
    parser.reset(
        ::bioparser::createParser<::bioparser::FastaParser, DataFormat>(
            input_file)
            .release());
  } else {
    parser.reset(
        ::bioparser::createParser<::bioparser::FastqParser, DataFormat>(
            input_file)
            .release());
  }

  while (parser->parse(sequences, 500 * 1024 * 1024))
    ;

  if constexpr (is_fasta) {
    return sequences;
  } else {
    Sequences<Sequence> retval(sequences.size());
    for (int i = 0; i < retval.size(); ++i)
      retval[i].reset(sequences[i].release());
    return retval;
  }
}

}  // namespace mapper