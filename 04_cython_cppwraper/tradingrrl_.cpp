#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <math.h>
#include <ctime>
#include "tradingrrl_.h"

using namespace std;

namespace tradingrrl{
    // Constructor
    cppTradingRRL::cppTradingRRL(int T, int M, int init_t, double  mu, double sigma, double rho, int n_epoch)
    :
    T(T),
    M(M),
    init_t(init_t),
    mu(mu),
    sigma(sigma),
    rho(rho),
    n_epoch(n_epoch),
    progress_period(100),
    q_threshold(0.7),
    x(T, vector<double>(M+2, 0.0)),
    F(T+1, 0.0),
    R(T, 0.0),
    w(M+2, 1.0),
    w_opt(M+2, 1.0),
    sumR(T, 0.0),
    sumR2(T, 0.0),
    A(0.0),
    B(0.0),
    S(0.0),
    S_opt(0.0),
    dSdA(0.0),
    dSdB(0.0),
    dAdR(0.0),
    dBdR(T, 0.0),
    dRdF(T, 0.0),
    dRdFp(T, 0.0),
    dFpdw(M+2, 0.0),
    dFdw(M+2, 0.0),
    dSdw(M+2, 0.0)
    {}
    
    // Destructor
    cppTradingRRL::~cppTradingRRL(){}
    
    // Menber functions
    vector<string> cppTradingRRL::split(string& input, char delimiter){
        istringstream stream(input);
        string field;
        vector<string> result;
        while (getline(stream, field, delimiter)) {
            result.push_back(field);
        }
        return result;
    }
    
    void cppTradingRRL::load_csv(string fname){
        ifstream ifs(fname.c_str());
        vector<string> tmp_t;
        vector<double> tmp_p;
        string tmp_tstr;
        double tmp_pdbl;
        int read_count =0;

        // --- read data
        string line;
        while (getline(ifs, line)){
            vector<string> strvec = split(line, ',');
            stringstream ss;
            tmp_tstr = strvec.at(0) + " " + strvec.at(1);
            ss << strvec.at(5);
            ss >> tmp_pdbl;
            tmp_t.push_back(tmp_tstr);
            tmp_p.push_back(tmp_pdbl);
            read_count +=1;
        }
        // --- Reverse chart data
        for( int i=0 ; i<read_count ; ++i ){
            all_t.push_back(tmp_t[read_count-1-i]);
            all_p.push_back(tmp_p[read_count-1-i]);
        }
    }
    
    int cppTradingRRL::quant(double f, double threshold){
        if(fabs(f)< threshold){
            return 0;
        }
        else if(f>0.0){
            return 1;
        }
        else{
            return -1;
        }
    }

    double cppTradingRRL::sign(double f){
        if(f==0.0){
            return 0.0;
        }
        else if(f>0.0){
            return 1.0;
        }
        else{
            return -1.0;
        }
    }

    void cppTradingRRL::set_t_p_r(){
        int i;
        for(i=0 ; i<T+M; ++i ){
            t.push_back(all_t[init_t+i]);
            p.push_back(all_p[init_t+i]);
            r.push_back(all_p[init_t+i] - all_p[init_t+i+1]);
        }
    }
    
    void cppTradingRRL::set_x_F(){
        for(int i=T-1 ; i>=0 ; --i ){
            x[i][0] = 1.0;
            x[i][M+2-1] = F[i+1];
            for(int j=1 ; j<M+2-1 ; ++j ){
                x[i][j] = r[i+(j-1)];
            }
            double wdotx = 0.0;
            for(int k=0 ; k<M+2 ; ++k ){
                wdotx += w[k]*x[i][k];
            }
            F[i] = tanh(wdotx);
        }
    }
    
    void cppTradingRRL::calc_R(){
        int i;
        for(i=T-1 ; i>=0 ; --i ){
            R[i]  = mu*(F[i+1]*r[i] - sigma*fabs(F[i] - F[i+1]));
        }
    }

    void cppTradingRRL::calc_sumR(){
        sumR[T-1]  = R[T-1];
        sumR2[T-1] = R[T-1]*R[T-1];
        for(int i=T-2 ; i>=0 ; --i ){
            sumR[i]  = sumR[i+1]  + R[i];
            sumR2[i] = sumR2[i+1] + R[i]*R[i];
        }
    }

    void cppTradingRRL::calc_dSdw(){
        set_x_F();
        calc_R();
        calc_sumR();
    
        A = sumR[0]/T;
        B = sumR2[0]/T;
        S = A/sqrt(B-A*A);
        dSdA = S*(1+S*S)/A;
        dSdB = -S*S*S/2/A/A;
        dAdR = 1.0/T;   
        
        for(int j=0 ; j<M+2; ++j ){
            dFpdw[j] = 0.0;
            dFdw[j]  = 0.0;
            dSdw[j]  = 0.0;
        }
        for(int i=T-1 ; i>=0 ; --i ){
            dBdR[i]  = 2.0/T*R[i];
            dRdF[i]  = - mu*sigma*sign(F[i] - F[i+1]);
            dRdFp[i] =   mu*r[i] + mu*sigma*sign(F[i] - F[i+1]);
            for(int j=0 ; j<M+2; ++j ){
                if(i!=T-1){
                    dFpdw[j] = dFdw[j];
                }
                dFdw[j]=(1.0 - F[i]*F[i])*(x[i][j] + w[M+2-1]*dFpdw[j]);
                dSdw[j]+=(dSdA*dAdR + dSdB*dBdR[i])*(dRdF[i]*dFdw[j] + dRdFp[i]*dFpdw[j]);
            }
        }
    }

    void cppTradingRRL::update_w(){
        for(int j=0 ; j<M+2; ++j ){
            w[j] += rho * dSdw[j];
        }
    }
    
    void cppTradingRRL::fit(){
        time_t tic, toc;
        int pre_epoch_times;
        int e_index;
        
        pre_epoch_times = epoch_S.size();

        calc_dSdw();
        cout<<"Epoch loop start. Initial sharp's ratio is "<<S<<"."<<endl;
        S_opt = S;
        
        tic = time(0);
        for(e_index=0 ; e_index<n_epoch; ++e_index ){
            calc_dSdw();
            if(S > S_opt){
                S_opt = S;
                w_opt = w;
            }
            epoch_S.push_back(S);
            update_w();
            if(e_index%progress_period == progress_period-1){
                toc = time(0); 
                cout << "Epoch: " << e_index + pre_epoch_times + 1 << " / " <<n_epoch + pre_epoch_times <<". Shape's ratio: "<<S<<". Elapsed time: "<< difftime(toc, tic) <<" sec."<< endl;
            }
        }
        toc = time(0);
        cout << "Epoch: " << e_index + pre_epoch_times << " / " <<n_epoch + pre_epoch_times <<". Shape's ratio: "<<S<<". Elapsed time: "<< difftime(toc, tic) <<" sec."<< endl;
        w = w_opt;
        calc_dSdw();
        cout << "Epoch loop end. Optimized sharp's ratio is "<<S_opt<<"."<< endl;
    }
    
    void cppTradingRRL::save_weight(){
        ofstream ofs1("w.csv"); 
        for(int j=0 ; j<M+2; ++j ){
            ofs1<<w[j]<<endl;
        }
        ofstream ofs2("epoch_S.csv");
        for(int e_index=0; e_index < (int)epoch_S.size(); ++e_index)
        {
            ofs2<<epoch_S[e_index]<<endl;
        }
        /*
        ofstream ofs3("r.csv");
        for(int j=0 ; j<T+M; ++j ){
            ofs3<<r[j]<<endl;
        }
        */
    }
    
    void cppTradingRRL::load_weight(){
        ifstream ifs("w.csv");
        vector<double> tmp_w;
        double tmp_wdbl;
        int read_count =0;
    
        // --- read data
        string line;
        while (getline(ifs, line)){
            stringstream ss;
            ss << line;
            ss >> tmp_wdbl;
            tmp_w.push_back(tmp_wdbl);
            read_count +=1;
        }
        w = tmp_w;
    }
}