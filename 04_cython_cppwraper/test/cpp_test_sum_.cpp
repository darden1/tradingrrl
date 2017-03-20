#include <vector>
#include "cpp_test_sum_.h"

using namespace std;

namespace cpp_test_sum{
    // Constructor
    cppTestSum::cppTestSum(int n_iter)
    :
    n_iter(n_iter),
    sum(n_iter, 0.0)
    {}
    // Destructor
    cppTestSum::~cppTestSum(){}
    // Menber function
    void cppTestSum::calc_sum(){
        for(int i=1; i<n_iter; ++i){
            sum[i] = sum[i-1] + i;
        }
    }
}
