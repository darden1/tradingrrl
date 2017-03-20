from libcpp.vector cimport vector

cdef extern from "cpp_test_sum_.h" namespace "cpp_test_sum":
    cdef cppclass cppTestSum:
        cppTestSum(int n_iter) except +
        int n_iter
        vector[double] sum
        
        void calc_sum()

cdef class TestSum(object):
    cdef cppTestSum* thisptr

    def __cinit__(self,int n_iter):
        self.thisptr = new cppTestSum(n_iter)
        
    def __dealloc(self):
        del self.thisptr

    def calc_sum(self):
        return self.thisptr.calc_sum()

    property n_iter:
        def __get__(self): return self.thisptr.n_iter
        def __set__(self, n_iter): self.thisptr.n_iter = n_iter

    property sum:
        def __get__(self): return self.thisptr.sum
        def __set__(self, sum): self.thisptr.sum = sum
