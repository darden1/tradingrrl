import time
import py_test_sum
import cy_test_sum

n_iter = 10000000

#--- Python sum
py_ts = py_test_sum.TestSum(n_iter)

tic = time.clock()
py_ts.calc_sum()
toc = time.clock()

print("Python Sum: " + str(py_ts.sum[-1])  + ". Elapsed time: " + str(toc-tic) + " sec.")

#--- Cython sum
cy_ts = cy_test_sum.TestSum(n_iter)

tic = time.clock()
cy_ts.calc_sum()
toc = time.clock()

print("Cython Sum: " + str(cy_ts.sum[-1])  + ". Elapsed time: " + str(toc-tic) + " sec.")