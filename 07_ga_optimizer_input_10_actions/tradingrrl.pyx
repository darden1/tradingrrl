from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "tradingrrl_.h" namespace "tradingrrl":
    cdef cppclass cppTradingRRL:
        cppTradingRRL(int T, int M, int init_t, double  mu, double sigma, double rho, int n_epoch) except +
        int T
        int M
        int init_t
        double  mu
        double sigma
        double rho
        int n_epoch
        int progress_period
        double q_threshold
        vector[string] all_t
        vector[double] all_p
        vector[string] t
        vector[double] p
        vector[double] r
        vector[ vector[double] ] x
        vector[double] F
        vector[double] R
        vector[double] w
        vector[double] w_opt
        vector[double] epoch_S
        vector[double] sumR
        vector[double] sumR2
        double A
        double B
        double S
        double S_opt
        double dSdA
        double dSdB
        double dAdR
        vector[double] dBdR
        vector[double] dRdF
        vector[double] dRdFp
        vector[double] dFpdw
        vector[double] dFdw
        vector[double] dSdw

        void load_csv(string fname)
        int quant(double f, double threshold)
        double sign(double f)
        void set_t_p_r()
        void set_x_F()
        void calc_R()
        void calc_sumR()
        double calc_S(vector[double] _w);
        void calc_dSdw()
        void update_w()
        void fit()
        void save_weight()
        void load_weight()
      
cdef class TradingRRL(object):
    cdef cppTradingRRL* crrl

    def __cinit__(self,int T=1000, int M=200, int init_t=10000, double  mu=10000, double sigma=0.04, double rho=1.0, int n_epoch=10000):
        self.crrl = new cppTradingRRL(T, M, init_t, mu, sigma, rho, n_epoch)
        
    def __dealloc(self):
        del self.crrl

    def load_csv(self, str fname):
        return self.crrl.load_csv(fname)

    def quant(self, double f, double threshold):
        return self.crrl.quant(f, threshold)

    def sign(self, double f):
        return self.crrl.sign(f)

    def set_t_p_r(self):
        return self.crrl.set_t_p_r()

    def set_x_F(self):
        return self.crrl.set_x_F()

    def calc_R(self):
        return self.crrl.calc_R()

    def calc_sumR(self):
        return self.crrl.calc_sumR()

    def calc_S(self, vector[double] _w):
        return self.crrl.calc_S(_w)

    def calc_dSdw(self):
        return self.crrl.calc_dSdw()

    def update_w(self):
        return self.crrl.update_w()

    def fit(self):
        return self.crrl.fit()

    def save_weight(self):
        return self.crrl.save_weight()

    def load_weight(self):
        return self.crrl.load_weight()

    property T:
        def __get__(self): return self.crrl.T
        def __set__(self, T): self.crrl.T = T

    property M:
        def __get__(self): return self.crrl.M
        def __set__(self, M): self.crrl.M = M

    property init_t:
        def __get__(self): return self.crrl.init_t
        def __set__(self, init_t): self.crrl.init_t = init_t

    property mu:
        def __get__(self): return self.crrl.mu
        def __set__(self, mu): self.crrl.mu = mu

    property sigma:
        def __get__(self): return self.crrl.sigma
        def __set__(self, sigma): self.crrl.sigma = sigma

    property rho:
        def __get__(self): return self.crrl.rho
        def __set__(self, rho): self.crrl.rho = rho

    property n_epoch:
        def __get__(self): return self.crrl.n_epoch
        def __set__(self, n_epoch): self.crrl.n_epoch = n_epoch

    property progress_period:
        def __get__(self): return self.crrl.progress_period
        def __set__(self, progress_period): self.crrl.progress_period = progress_period

    property q_threshold:
        def __get__(self): return self.crrl.q_threshold
        def __set__(self, q_threshold): self.crrl.q_threshold = q_threshold

    property all_t:
        def __get__(self): return self.crrl.all_t
        def __set__(self, all_t): self.crrl.all_t = all_t


    property all_p:
        def __get__(self): return self.crrl.all_p
        def __set__(self, all_p): self.crrl.all_p = all_p

    property t:
        def __get__(self): return self.crrl.t
        def __set__(self, t): self.crrl.t = t

    property p:
        def __get__(self): return self.crrl.p
        def __set__(self, p): self.crrl.p = p

    property r:
        def __get__(self): return self.crrl.r
        def __set__(self, r): self.crrl.r = r

    property x:
        def __get__(self): return self.crrl.x
        def __set__(self, x): self.crrl.x = x

    property F:
        def __get__(self): return self.crrl.F
        def __set__(self, F): self.crrl.F = F

    property R:
        def __get__(self): return self.crrl.R
        def __set__(self, R): self.crrl.R = R

    property w:
        def __get__(self): return self.crrl.w
        def __set__(self, w): self.crrl.w = w

    property w_opt:
        def __get__(self): return self.crrl.w_opt
        def __set__(self, w_opt): self.crrl.w_opt = w_opt

    property epoch_S:
        def __get__(self): return self.crrl.epoch_S
        def __set__(self, epoch_S): self.crrl.epoch_S = epoch_S

    property sumR:
        def __get__(self): return self.crrl.sumR
        def __set__(self, sumR): self.crrl.sumR = sumR

    property sumR2:
        def __get__(self): return self.crrl.sumR2
        def __set__(self, sumR2): self.crrl.sumR2 = sumR2

    property A:
        def __get__(self): return self.crrl.A
        def __set__(self, A): self.crrl.A = A

    property B:
        def __get__(self): return self.crrl.B
        def __set__(self, B): self.crrl.B = B

    property S:
        def __get__(self): return self.crrl.S
        def __set__(self, S): self.crrl.S = S

    property S_opt:
        def __get__(self): return self.crrl.S_opt
        def __set__(self, S_opt): self.crrl.S_opt = S_opt

    property dSdA:
        def __get__(self): return self.crrl.dSdA
        def __set__(self, dSdA): self.crrl.dSdA = dSdA

    property dSdB:
        def __get__(self): return self.crrl.dSdB
        def __set__(self, dSdB): self.crrl.dSdB = dSdB

    property dAdR:
        def __get__(self): return self.crrl.dAdR
        def __set__(self, dAdR): self.crrl.dAdR = dAdR

    property dBdR:
        def __get__(self): return self.crrl.dBdR
        def __set__(self, dBdR): self.crrl.dBdR = dBdR

    property dRdF:
        def __get__(self): return self.crrl.dRdF
        def __set__(self, dRdF): self.crrl.dRdF = dRdF

    property dRdFp:
        def __get__(self): return self.crrl.dRdFp
        def __set__(self, dRdFp): self.crrl.dRdFp = dRdFp

    property dFpdw:
        def __get__(self): return self.crrl.dFpdw
        def __set__(self, dFpdw): self.crrl.dFpdw = dFpdw

    property dFdw:
        def __get__(self): return self.crrl.dFdw
        def __set__(self, dFdw): self.crrl.dFdw = dFdw

    property dSdw:
        def __get__(self): return self.crrl.dSdw
        def __set__(self, dSdw): self.crrl.dSdw = dSdw