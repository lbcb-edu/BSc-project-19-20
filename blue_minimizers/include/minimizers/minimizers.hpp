#pragma once

#include <algorithm>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include <common/strong_type.hpp>

namespace blue {

using KMerInfo = ::std::tuple<unsigned, unsigned, bool>;

using SequenceLength = StrongType<unsigned, struct SequenceLengthTag>;
using KType = StrongType<unsigned, struct KTypeTag>;
using WindowLength = StrongType<unsigned, struct WindowLengthTag>;

::std::vector<KMerInfo> minimizers(const char* sequence,
                                   SequenceLength sequence_length, KType k,
                                   WindowLength window_length);

}  // namespace blue