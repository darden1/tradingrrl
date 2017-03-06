#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <math.h>
#include <cmath>
#include <ctime>

using namespace std;

class ChartData{
    public:
        vector<string> all_t;
        vector<double> all_p;

    private:
        vector<string> split(string& input, char delimiter)
        {
            istringstream stream(input);
            string field;
            vector<string> result;
            while (getline(stream, field, delimiter)) {
                result.push_back(field);
            }
            return result;
        }
        
    public:
        void load_csv(const char* fname)
        {
            ifstream ifs(fname);
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
};

class TradingRRL{
    public:
        int T;
        int M;
        int init_t;
        double q_threshold;
        double mu;
        double sigma;
        double alpha;
        int n_epoch;
        int progress_period;
        vector<string> t;
        vector<double> p;
        vector<double> r;
        vector< vector<double> > x;
        vector<double> F;
        vector<double> R;
        vector<double> w;
        vector<double> epoch_S;
        vector<double> sumR;
        vector<double> sumR2;
        double A;
        double B;
        double S;
        double dSdA;
        double dSdB;
        double dAdR;
        vector<double> dBdR;
        vector<double> dRdF;
        vector<double> dRdFp;
        vector<double> dFpdw;
        vector<double> dFdw;
        vector<double> dSdw;
        
        // ---  Constructor
        TradingRRL(int T=1000, int M=200, int init_t=10000, double q_threshold=0.7, double  mu=10000, double sigma=0.04, double alpha=1.0, int n_epoch=10000)
        :
        T(T),
        M(M),
        init_t(init_t),
        q_threshold(q_threshold),
        mu(mu),
        sigma(sigma),
        alpha(alpha),
        n_epoch(n_epoch),
        progress_period(100),
        x(T, vector<double>(M+2, 0.0)),
        F(T+1, 0.0),
        R(T, 0.0),
        w(M+2, 1.0),
        sumR(T+1, 0.0),
        sumR2(T+1, 0.0),
        A(0.0),
        B(0.0),
        S(0.0),
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
        
        int quant(double f, double threshold){
            if(abs(f)< threshold){
                return 0;
            }
            else if(f>0.0){
                return 1;
            }
            else{
                return -1;
            }
        }

        double sign(double f){
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
        
        void set_t_p_r(ChartData cd){
            int i;
            for(i=0 ; i<T+M; ++i ){
                t.push_back(cd.all_t[init_t+i]);
                p.push_back(cd.all_p[init_t+i]);
                r.push_back(cd.all_p[init_t+i] - cd.all_p[init_t+i+1]);
            }
        }
        
        void set_x_F(){
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
        
        void calc_R(){
            int i;
            for(i=T-1 ; i>=0 ; --i ){
                R[i]  = mu*(F[i+1]*r[i] - sigma*abs(F[i] - F[i+1]));
            }
        }
    
        void calc_sumR(){
            sumR[T-1]  = R[T-1];
            sumR2[T-1] = R[T-1]*R[T-1];
            for(int i=T-2 ; i>=0 ; --i ){
                sumR[i]  = sumR[i+1]  + R[i];
                sumR2[i] = sumR2[i+1] + R[i]*R[i];
            }
        }

        void calc_dSdw(){
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

        void update_w(){
            for(int j=0 ; j<M+2; ++j ){
                w[j] += alpha * dSdw[j];
            }
        }
        
        void fit(){
            time_t tic, toc;
            int pre_epoch_times;
            int e_index;
            
            pre_epoch_times = epoch_S.size();

            calc_dSdw();
            cout<<"Epoch loop start. Initial sharp's ratio is "<<S<<"."<<endl;
            
            tic = time(0);
            for(e_index=0 ; e_index<n_epoch; ++e_index ){
                calc_dSdw();
                update_w();
                epoch_S.push_back(S);
                if(e_index%progress_period == progress_period-1){
                    toc = time(0); 
                    cout << "Epoch: " << e_index + pre_epoch_times + 1 << " / " <<n_epoch + pre_epoch_times <<". Shape's ratio: "<<S<<". Elapsed time: "<< difftime(toc, tic) <<" sec."<< endl;
                }
            }
            toc = time(0); 
            cout << "Epoch loop end. Optimized sharp's ratio is "<<S<<"."<< endl;
            cout << "Epoch: " << e_index + pre_epoch_times << " / " <<n_epoch + pre_epoch_times <<". Shape's ratio: "<<S<<". Elapsed time: "<< difftime(toc, tic) <<" sec."<< endl;
        }
        
        void save_weights(){
            ofstream ofs1("w.csv"); 
            for(int j=0 ; j<M+2; ++j ){
                ofs1<<w[j]<<endl;
            }
            ofstream ofs2("epoch_S.csv");
            for(int e_index=0; e_index<epoch_S.size(); ++e_index)
            {
                ofs2<<epoch_S[e_index]<<endl;
            }
            ofstream ofs3("r.csv");
            for(int j=0 ; j<T+M; ++j ){
                ofs3<<r[j]<<endl;
            }
        }

};
    
int main()
{
    const char* fname = "../data/USDJPY30.csv";
    int init_t = 6000;
    int T = 1000;
    int M = 200;
    double q_threshold = 0.7;
    double mu = 10000;
    double sigma = 0.04;
    double alpha = 1.0;
    int n_epoch = 10000;
    
    ChartData cd;
    cd.load_csv(fname);
    
    TradingRRL rrl(T, M, init_t, q_threshold, mu, sigma, alpha, n_epoch);
    rrl.set_t_p_r(cd);
    rrl.fit();
    rrl.save_weights();

}