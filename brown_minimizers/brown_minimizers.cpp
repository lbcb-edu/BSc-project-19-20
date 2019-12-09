#include "brown_minimizers.hpp"
#include <vector>
#include <tuple>
#include <iostream>
#include <deque>
#include <math.h>

using namespace std;

struct S {
 unsigned int b : 2;
};

namespace brown {

    vector<tuple<unsigned int, unsigned int, bool>> minimizers(const char* sequence, unsigned int sequence_length,
                                                                    unsigned int k,
                                                                    unsigned int window_length) {
        vector<S> codedSequence;
        vector<tuple<unsigned int, unsigned int, bool>> mins;
        for (int i = 0; i < sequence_length; i++)
        {
            if(sequence[i] == 'A') {
                S s = {0};
                codedSequence.push_back(s);
            }
            else if(sequence[i] == 'C') {
                S s = {1};
                codedSequence.push_back(s);
            }
            else if(sequence[i] == 'G') {
                S s = {2};
                codedSequence.push_back(s);
            }
            else if(sequence[i] == 'T') {
                S s = {3};
                codedSequence.push_back(s);
            }
        }
        deque<int> window;
        int smallest = -1;
        int cnt = 0;
        int pos = 0;
        for (int i = 0; i < sequence_length; i++) {
            if((i + 1) <= window_length + k - 1) {
                window.push_back((int)codedSequence[i].b);
                //cout << window[i] << "\n";
                if((i + 1) >= k) {
                    int sum = 0;
                    for (int j = 0; j < k; j++)
                    {
                        sum += pow(10, k-j-1)*window[j + cnt];
                    }
                    //cout << sum << " ";
                    if(smallest == -1 || sum < smallest) {
                        smallest = sum;
                        pos = i + 1 - k;
                    } 
                    cnt++;
                }
                if(cnt == window_length) mins.push_back(make_tuple(smallest, pos, true));
            } else {
                window.pop_front();
                window.push_back((int)codedSequence[i].b);
            }
        }
        // for (int i = 0; i < mins.size(); i++)
        // {
        //     cout << get<0>(mins[i]) << " " << get<1>(mins[i]) << " " << get<2>(mins[i]) << "\n"; 
        // }
        
        return mins;
        
        }
}