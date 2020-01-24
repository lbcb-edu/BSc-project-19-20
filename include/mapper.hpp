#pragma once

#include <algorithm>
#include <bioparser/bioparser.hpp>
#include <iterator>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include "sequence.hpp"

namespace mapper {

template <typename T>
using Sequences = ::std::vector<::std::unique_ptr<T>>;

template <template <typename> typename Parser>
Sequences<Sequence> Parse(const ::std::string& input_file) {
  Sequences<Sequence> sequences;
  ::bioparser::createParser<Parser, Sequence>(input_file)->parse(sequences, -1);
  return sequences;
}

}  // namespace mapper