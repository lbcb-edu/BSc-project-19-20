#pragma once

#include <algorithm>
#include <bioparser/bioparser.hpp>
#include <memory>
#include <type_traits>
#include <vector>

#include "sequence.hpp"

namespace mapper {

namespace detail {

template <typename DataFormat>
auto GetParser(const ::std::string& input_file) {
  if constexpr (::std::is_same_v<DataFormat, Sequence>)
    return ::bioparser::createParser<::bioparser::FastaParser, DataFormat>(
        input_file);
  else
    return ::bioparser::createParser<::bioparser::FastqParser, DataFormat>(
        input_file);
}

}  // namespace detail

template <typename T>
using Sequences = ::std::vector<::std::unique_ptr<T>>;

template <typename DataFormat>
Sequences<Sequence> Parse(const ::std::string& input_file) {
  Sequences<DataFormat> sequences;
  detail::GetParser<DataFormat>(input_file)->parse(sequences, -1);

  if constexpr (std::is_same_v<DataFormat, Sequence>) {
    return sequences;
  } else {
    Sequences<Sequence> retval(sequences.size());
    ::std::move(sequences.begin(), sequences.end(), retval.begin());
    return retval;
  }
}

}  // namespace mapper