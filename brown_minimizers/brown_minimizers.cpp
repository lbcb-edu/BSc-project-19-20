#include "brown_minimizers.hpp"
#include <math.h>
#include <deque>
#include <iostream>
#include <tuple>
#include <vector>

using namespace std;

struct S {
  unsigned int b : 2;
};

namespace brown {

vector<tuple<unsigned int, unsigned int, bool>> minimizers(
    const char* sequence, unsigned int sequence_length, unsigned int k,
    unsigned int window_length) {
  vector<S> codedSequence;
  vector<tuple<unsigned int, unsigned int, bool>> mins;
  for (int i = 0; i < sequence_length; i++) {
    if (sequence[i] == 'A') {
      S s = {0};
      codedSequence.push_back(s);
    } else if (sequence[i] == 'C') {
      S s = {1};
      codedSequence.push_back(s);
    } else if (sequence[i] == 'G') {
      S s = {2};
      codedSequence.push_back(s);
    } else if (sequence[i] == 'T') {
      S s = {3};
      codedSequence.push_back(s);
    }
  }
  deque<unsigned int> window;
  deque<unsigned int> sums;
  unsigned int smallest = -1;
  unsigned int smallPos = -1;
  unsigned int previous = -1;
  unsigned int pos = 0;
  for (int i = 0; i < sequence_length; i++) {
    if ((i + 1) <= window_length + k - 1) {
      window.push_back((unsigned int)codedSequence[i].b);
      if ((i + 1) >= k) {
        unsigned int sum = 0;
        for (int j = 0; j < k; j++) {
          sum += pow(10, k - j - 1) * window[j + i + 1 - k];
        }
        sums.push_back(sum);
        if (smallest == -1 || sum < smallest) {
          smallest = sum;
          pos = i + 1 - k;
        }
        if (i + 1 == window_length + k - 1) {
          mins.push_back(make_tuple(smallest, pos, true));
          previous = smallest;
        }
      }
    } else {
      window.pop_front();
      window.push_back((unsigned int)codedSequence[i].b);
      sums.pop_front();
      unsigned int sum = 0;

      for (int j = 0; j < k; j++) {
        sum += pow(10, k - j - 1) * window[j + window_length - 1];
      }
      sums.push_back(sum);
      smallest = -1;
      for (int j = 0; j < sums.size(); j++) {
        if (smallest == -1 || sums[j] < smallest) {
          smallest = sums[j];
          pos = i + 1 - k;
        }
      }
      if (smallest != previous) {
        mins.push_back(make_tuple(smallest, pos, true));
        previous = smallest;
      }
    }
  }

  return mins;
}
}  // namespace brown