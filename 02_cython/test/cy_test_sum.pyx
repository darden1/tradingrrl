import numpy as np
cimport numpy as np
cimport cython

ctypedef np.float64_t DOUBLE_t

cdef class TestSum(object):

    cdef public int n_iter
    cdef public np.ndarray sum
    
    def __init__(self, n_iter):
        self.n_iter = n_iter
        self.sum = np.zeros(n_iter, dtype=np.float64)
        
    def calc_sum(self):
        cdef int i
        cdef np.ndarray[DOUBLE_t, ndim=1, mode="c"] sum
        sum = self.sum
        for i in range(1, self.n_iter):
            sum[i] = sum[i-1] + i
        self.sum = sum
