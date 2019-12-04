#pragma once

#include <utility>

namespace blue {

template <typename T, typename Tag>
class StrongType {
 private:
  T value_;

 public:
  explicit StrongType(const T& value) : value_(value) {}
  explicit StrongType(T&& value) : value_(::std::move(value)) {}

  T& get() {
    return value_;
  }

  const T& get() const {
    return value_;
  }
};

}  // namespace algn