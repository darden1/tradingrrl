#ifndef TRADINGRRL_H
#define TRADINGRRL_H

#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <math.h>
#include <ctime>

using namespace std;

namespace tradingrrl{
    class cppTradingRRL{
        public:
            int T;
            int M;
            int init_t;
            double mu;
            double sigma;
            double rho;
            int n_epoch;
            int progress_period;
            double q_threshold;
            vector<string> all_t;
            vector<double> all_p;
            vector<string> t;
            vector<double> p;
            vector<double> r;
            vector< vector<double> > x;
            vector<double> F;
            vector<double> R;
            vector<double> w;
            vector<double> w_opt;
            vector<double> epoch_S;
            vector<double> sumR;
            vector<double> sumR2;
            double A;
            double B;
            double S;
            double S_opt;
            double dSdA;
            double dSdB;
            double dAdR;
            vector<double> dBdR;
            vector<double> dRdF;
            vector<double> dRdFp;
            vector<double> dFpdw;
            vector<double> dFdw;
            vector<double> dSdw;
    
        private:
            vector<string> split(string& input, char delimiter);

        public:
            // Constructor
            cppTradingRRL(int T, int M, int init_t, double  mu, double sigma, double rho, int n_epoch);
            // Destructor
            ~cppTradingRRL();
            // Menber functions
            void load_csv(string fname);
            int quant(double f, double threshold);
            double sign(double f);
            void set_t_p_r();
            void set_x_F();
            void calc_R();
            void calc_sumR();
            void calc_dSdw();
            void update_w();
            void fit();
            void save_weight();
            void load_weight();
    };
}
#endif