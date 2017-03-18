import numpy as np

class TestSum(object):

    def __init__(self, n_iter):
        self.n_iter = n_iter
        self.sum = np.zeros(n_iter)
        
    def calc_sum(self):
        for i in range(1, self.n_iter):
            self.sum[i] = self.sum[i-1] + i
