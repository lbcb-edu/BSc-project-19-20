#pragma once

#include <string>
#include <utility>

namespace mapper {

struct Sequence {
  ::std::string name, sequence;
  ::std::string quality = "255";

  Sequence(const char* name, ::std::uint32_t name_len, const char* sequence,
           ::std::uint32_t seq_len)
      : name(name, name_len), sequence(sequence, seq_len) {}

  Sequence(const char* name, ::std::uint32_t name_len, const char* sequence,
           ::std::uint32_t seq_len, const char* quality, ::std::uint32_t q_len)
      : name(name, name_len),
        sequence(sequence, seq_len),
        quality(quality, q_len) {}
};

}  // namespace mapper