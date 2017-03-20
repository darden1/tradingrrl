#ifndef CPPTESTSUM_H
#define CPPTESTSUM_H

#include <vector>
using namespace std;

namespace cpp_test_sum{
    class cppTestSum{
        public:
            int n_iter;
            vector<double> sum;
            // Constructor
            cppTestSum(int n_iter);
            // Destructor
            ~cppTestSum();
            // Menber function
            void calc_sum();
    };
}
#endif