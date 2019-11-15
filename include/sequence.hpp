#pragma once

#include <string>
#include <utility>

namespace mapper {

struct Sequence {
  ::std::string name, sequence;
  Sequence(const char* name, ::std::uint32_t name_len, const char* sequence,
           ::std::uint32_t seq_len)
      : name(name, name_len), sequence(sequence, seq_len) {}
};

struct QSequence : Sequence {
  ::std::string quality;
  QSequence(const char* name, ::std::uint32_t name_len, const char* sequence,
            ::std::uint32_t seq_len, const char* quality, ::std::uint32_t q_len)
      : Sequence(name, name_len, sequence, seq_len), quality(quality, q_len) {}
};

}  // namespace mapper