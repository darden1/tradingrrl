# -*- coding: utf-8 -*-
#cython: boundscheck=False
#cython: wraparound=False
#cython: nonecheck=False

import numpy as np
cimport numpy as np
cimport cython

ctypedef np.float64_t DOUBLE_t
ctypedef np.int_t INT_t

import time
import pickle
import pandas as pd
from datetime import datetime as dt
import matplotlib.pyplot as plt

cdef class TradingRRL(object):
    cdef public int T
    cdef public int M
    cdef public int init_t
    cdef public double mu
    cdef public double sigma
    cdef public double rho
    cdef public np.ndarray all_t
    cdef public np.ndarray all_p
    cdef public np.ndarray t
    cdef public np.ndarray p
    cdef public np.ndarray r
    cdef public np.ndarray x
    cdef public np.ndarray F
    cdef public np.ndarray R
    cdef public np.ndarray w
    cdef public np.ndarray w_opt
    cdef public np.ndarray epoch_S
    cdef public int n_epoch
    cdef public int progress_period
    cdef public double q_threshold
    cdef public np.ndarray sumR
    cdef public np.ndarray sumR2
    cdef public double A
    cdef public double B
    cdef public double S
    cdef public double S_opt
    cdef public double dSdA
    cdef public double dSdB
    cdef public double dAdR
    cdef public np.ndarray dBdR
    cdef public np.ndarray dRdF
    cdef public np.ndarray dRdFp
    cdef public np.ndarray dFpdw
    cdef public np.ndarray dFdw
    cdef public np.ndarray dSdw

    def __init__(self, T=1000, M=200, init_t=10000, mu=10000, sigma=0.04, rho=1.0, n_epoch=10000):
        self.T = T
        self.M = M
        self.init_t = init_t
        self.mu = mu
        self.sigma = sigma
        self.rho = rho
        self.all_t = None
        self.all_p = None
        self.t = None
        self.p = None
        self.r = None
        self.x = np.zeros([T,M+2], dtype=np.float64)
        self.F = np.zeros(T+1, dtype=np.float64)
        self.R = np.zeros(T, dtype=np.float64)
        self.w = np.ones(M+2, dtype=np.float64)
        self.w_opt = np.ones(M+2, dtype=np.float64)
        self.epoch_S = np.empty(0, dtype=np.float64)
        self.n_epoch = n_epoch
        self.progress_period = 100
        self.q_threshold = 0.7
        self.sumR = np.zeros(T, dtype=np.float64)
        self.sumR2 = np.zeros(T, dtype=np.float64)
        self.A = 0.0
        self.B = 0.0
        self.S = 0.0
        self.S_opt = 0.0
        self.dSdA = 0.0
        self.dSdB = 0.0
        self.dAdR = 0.0
        self.dBdR = np.zeros(T, dtype=np.float64)
        self.dRdF = np.zeros(T, dtype=np.float64)
        self.dRdFp = np.zeros(T, dtype=np.float64)
        self.dFpdw = np.zeros(M+2, dtype=np.float64)
        self.dFdw = np.zeros(M+2, dtype=np.float64)
        self.dSdw = np.zeros(M+2, dtype=np.float64)

    def load_csv(self, fname):
        tmp = pd.read_csv(fname, header=None)
        tmp_tstr = tmp[0] +" " + tmp[1]
        tmp_t = [dt.strptime(tmp_tstr[i], '%Y.%m.%d %H:%M') for i in range(len(tmp_tstr))]
        tmp_p = list(tmp[5])
        self.all_t = np.array(tmp_t[::-1])
        self.all_p = np.array(tmp_p[::-1])

    def quant(self, f):
        fc = f.copy()
        fc[np.where(np.abs(fc) < self.q_threshold)] = 0
        return np.sign(fc)

    def set_t_p_r(self):
        self.t = self.all_t[self.init_t:self.init_t+self.T+self.M+1]
        self.p = self.all_p[self.init_t:self.init_t+self.T+self.M+1]
        self.r = -np.diff(self.p)

    def set_x_F(self):
        cdef int i
        cdef int j
        cdef np.ndarray[DOUBLE_t, ndim=1, mode="c"] xt
        cdef np.ndarray[DOUBLE_t, ndim=2, mode="c"] x
        cdef np.ndarray[DOUBLE_t, ndim=1, mode="c"] r
        cdef np.ndarray[DOUBLE_t, ndim=1, mode="c"] F
        cdef np.ndarray[DOUBLE_t, ndim=1, mode="c"] w
        x = self.x
        r = self.r
        F = self.F
        w = self.w
        for i in xrange(self.T-1, -1, -1):
            """
            x[i] = np.zeros(self.M+2, dtype=np.float64)
            x[i][0] = 1.0
            x[i][self.M+2-1] = F[i+1]
            for j in xrange(1, self.M+2-1, 1):
                x[i][j] = r[i+j-1]
            """
            xt = np.zeros(self.M+2, dtype=np.float64)
            xt[0] = 1.0
            xt[self.M+2-1] = F[i+1]
            for j in xrange(1, self.M+2-1, 1):
                xt[j] = r[i+j-1]
            x[i] = xt
            F[i] = np.tanh(np.dot(w, x[i]))
            
        self.x = x
        self.F = F

    def calc_R(self):
        self.R = self.mu*(self.F[1:]*self.r[:self.T] - self.sigma*np.abs(-np.diff(self.F)))

    def calc_sumR(self):
        self.sumR = np.cumsum(self.R[::-1])[::-1].copy(order='C')
        self.sumR2 = np.cumsum((self.R**2)[::-1])[::-1].copy(order='C')

    def calc_dSdw(self):
        cdef int i
        cdef np.ndarray[DOUBLE_t, ndim=2, mode="c"] x
        cdef np.ndarray[DOUBLE_t, ndim=1, mode="c"] r
        cdef np.ndarray[DOUBLE_t, ndim=1, mode="c"] F
        cdef np.ndarray[DOUBLE_t, ndim=1, mode="c"] R
        cdef np.ndarray[DOUBLE_t, ndim=1, mode="c"] w
        cdef np.ndarray[DOUBLE_t, ndim=1, mode="c"] sumR
        cdef np.ndarray[DOUBLE_t, ndim=1, mode="c"] sumR2
        cdef double A
        cdef double B
        cdef double S
        cdef double dSdA
        cdef double dSdB
        cdef double dAdR
        cdef np.ndarray[DOUBLE_t, ndim=1, mode="c"] dBdR
        cdef np.ndarray[DOUBLE_t, ndim=1, mode="c"] dRdF
        cdef np.ndarray[DOUBLE_t, ndim=1, mode="c"] dRdFp
        cdef np.ndarray[DOUBLE_t, ndim=1, mode="c"] dFpdw
        cdef np.ndarray[DOUBLE_t, ndim=1, mode="c"] dFdw
        cdef np.ndarray[DOUBLE_t, ndim=1, mode="c"] dSdw

        self.set_x_F()
        self.calc_R()
        self.calc_sumR()
        x = self.x
        r = self.r
        F = self.F
        R = self.R
        w = self.w
        sumR = self.sumR
        sumR2 = self.sumR2

        A = sumR[0]/self.T
        B = sumR2[0]/self.T
        S = A/np.sqrt(B-A**2)
        dSdA = S * (1 + S**2) / A 
        dSdB = -S**3 / 2 / A**2
        dAdR = 1.0/self.T
        dBdR = 2.0/self.T*R
        dRdF = -self.mu*self.sigma*np.sign(-np.diff(F))
        dRdFp = self.mu*r[:self.T] + self.mu*self.sigma*np.sign(-np.diff(F))

        dFpdw = np.zeros(self.M+2, dtype=np.float64)
        dFdw = np.zeros(self.M+2, dtype=np.float64)
        dSdw = np.zeros(self.M+2, dtype=np.float64)
        for i in xrange(self.T-1, -1, -1):
            if i!=self.T-1:
                dFpdw = dFdw.copy()
            dFdw  = (1 - F[i]**2)*(x[i] + w[self.M+2-1]*dFpdw)
            dSdw += (dSdA*dAdR + dSdB*dBdR[i])*(dRdF[i]*dFdw + dRdFp[i]*dFpdw)

        self.sumR = sumR
        self.sumR2 = sumR2
        self.A = A
        self.B = B
        self.S = S
        self.dSdA = dSdA
        self.dSdB = dSdB
        self.dAdR = dAdR
        self.dBdR = dBdR
        self.dRdF = dRdF
        self.dRdFp = dRdFp
        self.dFpdw = dFpdw
        self.dFdw = dFdw
        self.dSdw = dSdw

    def update_w(self):
        self.w += self.rho*self.dSdw

    def fit(self):
        cdef double tic
        cdef double toc
        cdef int pre_epoch_times
        cdef int e_index

        pre_epoch_times = len(self.epoch_S)

        self.calc_dSdw()
        print("Epoch loop start. Initial sharp's ratio is " + str(self.S) + ".")
        self.S_opt = self.S
        
        tic = time.clock()
        for e_index in range(self.n_epoch):
            self.calc_dSdw()
            if self.S > self.S_opt:
                self.S_opt = self.S
                self.w_opt = self.w.copy()
            self.epoch_S = np.append(self.epoch_S, self.S)
            self.update_w()
            if e_index % self.progress_period  == self.progress_period-1:
                toc = time.clock()
                print("Epoch: " + str(e_index + pre_epoch_times + 1) + "/" + str(self.n_epoch + pre_epoch_times) +". Shape's ratio: " + str(self.S) + ". Elapsed time: " + str(toc-tic) + " sec.")
        toc = time.clock()
        print("Epoch: " + str(e_index + pre_epoch_times + 1) + "/" + str(self.n_epoch + pre_epoch_times) +". Shape's ratio: " + str(self.S) + ". Elapsed time: " + str(toc-tic) + " sec.")
        self.w = self.w_opt.copy()
        self.calc_dSdw()
        print("Epoch loop end. Optimized sharp's ratio is " + str(self.S_opt) + ".")

    def save_weight(self):
        pd.DataFrame(self.w).to_csv("w.csv", header=False, index=False)
        pd.DataFrame(self.epoch_S).to_csv("epoch_S.csv", header=False, index=False)
        
    def load_weight(self):
        tmp = pd.read_csv("w.csv", header=None)
        self.w = tmp.T.values[0]


def plot_hist(n_tick, R):
    rnge = max(R)-min(R)
    tick = rnge/n_tick
    tick_min = [min(R)-tick*0.5 +i*tick for i in range(n_tick)]
    tick_max = [min(R)+tick*0.5 +i*tick for i in range(n_tick)]
    tick_center = [min(R) +i*tick for i in range(n_tick)]
    tick_val=[0.0]*n_tick
    for i in range(n_tick ):
        tick_val[i]=len(set(np.where(tick_min[i]<np.array(R))[0].tolist()).intersection(np.where(np.array(R)<=tick_max[i])[0]))
    plt.bar(tick_center,tick_val, width=tick)
    plt.grid()
    plt.show()
