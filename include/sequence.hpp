#pragma once

#include <string>
#include <utility>

namespace mapper {

struct Sequence {
  ::std::string name, sequence;
  Sequence(const char* name, ::std::uint32_t, const char* sequence,
           ::std::uint32_t)
      : name(name), sequence(name) {}
};

struct QSequence : Sequence {
  ::std::string quality;
  QSequence(const char* name, ::std::uint32_t i1, const char* sequence,
            ::std::uint32_t i2, const char* quality, ::std::uint32_t)
      : Sequence(name, i1, sequence, i2), quality(quality) {}
};

}  // namespace mapper