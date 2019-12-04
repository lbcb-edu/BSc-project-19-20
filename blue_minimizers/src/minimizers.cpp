#include "minimizers/minimizers.hpp"

namespace blue {

::std::vector<KMerInfo> minimizers(const char* sequence,
                                   SequenceLength sequence_length, KType k,
                                   WindowLength window_length) {
  ::std::vector<KMerInfo> ret;
  ret.emplace_back(k.get(), window_length.get(), true);

  return ret;
}

}  // namespace blue