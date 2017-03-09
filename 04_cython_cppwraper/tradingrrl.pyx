from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "tradingrrl_.h" namespace "tradingrrl":
    cdef cppclass cppTradingRRL:
        cppTradingRRL(int T, int M, int init_t, double  mu, double sigma, double rho, int n_epoch)
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
        vector[double] epoch_S
        vector[double] sumR
        vector[double] sumR2
        double A
        double B
        double S
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
        void set_all_t_p(vector[string] _all_pt, vector[double] _all_p)
        void set_w(vector[double] _w)
        void set_t_p_r()
        void set_x_F()
        void calc_R()
        void calc_sumR()
        void calc_dSdw()
        void update_w()
        void fit()
        void save_weight()
        void load_weight()
      
cdef class TradingRRL(object):
    cdef cppTradingRRL* crrl
    
    cdef public int T
    cdef public int M
    cdef public int init_t
    cdef public double mu
    cdef public double sigma
    cdef public double rho
    cdef public int n_epoch
    cdef public int progress_period
    cdef public double q_threshold
    cdef public vector[string] all_t
    cdef public vector[double] all_p
    cdef public vector[string] t
    cdef public vector[double] p
    cdef public vector[double] r
    cdef public vector[ vector[double] ] x
    cdef public vector[double] F
    cdef public vector[double] R
    cdef public vector[double] w
    cdef public vector[double] epoch_S
    cdef public vector[double] sumR
    cdef public vector[double] sumR2
    cdef public double A
    cdef public double B
    cdef public double S
    cdef public double dSdA
    cdef public double dSdB
    cdef public double dAdR
    cdef public vector[double] dBdR
    cdef public vector[double] dRdF
    cdef public vector[double] dRdFp
    cdef public vector[double] dFpdw
    cdef public vector[double] dFdw
    cdef public vector[double] dSdw
   
    def __cinit__(self,int T=1000, int M=200, int init_t=10000, double  mu=10000, double sigma=0.04, double rho=1.0, int n_epoch=10000):
        self.crrl = new cppTradingRRL(T, M, init_t, mu, sigma, rho, n_epoch)
        self.T = T
        self.M = M
        self.init_t = init_t
        self.mu = mu
        self.sigma = sigma
        self.rho = rho
        
    def __dealloc(self):
        del self.crrl
   
    def update_val(self):
        self.w = self.crrl.w
        self.all_t = self.crrl.all_t
        self.all_p = self.crrl.all_p
        self.t = self.crrl.t 
        self.p = self.crrl.p
        self.r = self.crrl.r
        self.x = self.crrl.x
        self.F = self.crrl.F
        self.R = self.crrl.R
        self.w = self.crrl.w 
        self.epoch_S = self.crrl.epoch_S
        self.n_epoch = self.crrl.n_epoch
        self.progress_period = self.crrl.progress_period 
        self.q_threshold = self.crrl.q_threshold
        self.sumR = self.crrl.sumR
        self.sumR2 =self.crrl.sumR2 
        self.A = self.crrl.A
        self.B = self.crrl.B
        self.S = self.crrl.S
        self.dSdA = self.crrl.dSdA
        self.dSdB = self.crrl.dSdB
        self.dAdR = self.crrl.dAdR
        self.dBdR = self.crrl.dBdR
        self.dRdF = self.crrl.dRdF
        self.dRdFp = self.crrl.dRdFp
        self.dFpdw = self.crrl.dFpdw
        self.dFdw = self.crrl.dFdw
        self.dSdw = self.crrl.dSdw

    def load_csv(self, str fname):
        return self.crrl.load_csv(fname)

    def quant(self, double f, double threshold):
        return self.crrl.quant(f, threshold)

    def sign(self, double f):
        return self.crrl.sign(f)

    def set_all_t_p(self, vector[string] _all_t, vector[double] _all_p):
        return self.crrl.set_all_t_p(_all_t, _all_p)

    def set_w(self, vector[double] _w):
        return self.crrl.set_w(_w)

    def set_t_p_r(self):
        return self.crrl.set_t_p_r()

    def set_x_F(self):
        return self.crrl.set_x_F()

    def calc_R(self):
        return self.crrl.calc_R()

    def calc_sumR(self):
        return self.crrl.calc_sumR()

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