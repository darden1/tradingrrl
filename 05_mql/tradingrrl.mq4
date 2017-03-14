//+------------------------------------------------------------------+
//|                                                   tradingrrl.mq4 |
//|                        Copyright 2014, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "Copyright 2014, MetaQuotes Software Corp."
#property link      "http://www.mql5.com"
#property version   "1.03"
#property strict

#define MAGIC 20170214

extern int    slippage        = 2;
extern double lots            = 0.1;
extern bool   stop_flag       = False;
extern bool   limit_flag      = False;
extern double stop_points     = 1000;
extern double limit_points    = 1000;

extern int    T               = 150;
extern int    M               = 30;
extern double mu              = 10000;
extern double sigma           = 0.04;
extern double rho             = 2.0;
extern double q_threshold     = 0.7;
extern int    n_epoch         = 10000;
extern int    n_tick_update_w = 100;
extern bool   write_log       = True;

int    n_tick;
bool   init_flag;
double Fp;

//+------------------------------------------------------------------+
//| RRL agent                                                        |
//+------------------------------------------------------------------+
class TradingRRL{
    public:
        /*
        int T;
        int M;
        double mu;
        double sigma;
        double rho;
        int n_epoch;
        double q_threshold;
        */
        int progress_period;

        double r[];
        double x[];
        double F[];
        double R[];
        double w[];
        double w_opt[];
        double epoch_S[];
        double sumR[];
        double sumR2[];
        double A;
        double B;
        double S;
        double S_opt;
        double dSdA;
        double dSdB;
        double dAdR;
        double dBdR[];
        double dRdF[];
        double dRdFp[];
        double dFpdw[];
        double dFdw[];
        double dSdw[];
        
        // Constructor
        //TradingRRL(int T, int M, double mu, double sigma, double rho, int n_epoch, double q_threshold);
        TradingRRL();
        // Destructor
        ~TradingRRL();
        // Menber functions
        int quant(double f);
        double sign(double f);
        double dot(double &a[], double &b[]);
        double tanh(double f);
        void set_r();
        void set_x(int t_index, double Fp);
        double calc_F(int t_index, double Fp);
        void set_F();
        void calc_R();
        void calc_sumR();
        void calc_dSdw();
        void update_w();
        void fit();
        void save_weight();
        void load_weight();
};

// Constructor
//TradingRRL::TradingRRL(int T, int M, double mu, double sigma, double rho, int n_epoch, double q_threshold){
TradingRRL::TradingRRL(){
   /*
   // Member variables can't be initialized, when the names of the member variables are the same with the ones of the global variables.
   T = T;
   M = M;
   mu =mu;
   sigma = sigma;
   rho = rho;
   n_epoch = n_epoch;
   q_threshold = q_threshold;
   */
   progress_period = 100;
   A = 0.0;
   B = 0.0;
   S = 0.0;
   S_opt = 0.0;
   dSdA = 0.0;
   dSdB = 0.0;
   dAdR = 0.0;

   ArrayResize(r, T+M);
   ArrayResize(x, M+2);
   ArrayResize(F, T+1);
   ArrayResize(R, T);
   ArrayResize(w, M+2);
   ArrayResize(w_opt, M+2);
   ArrayResize(sumR, T+1);
   ArrayResize(sumR2, T+1);
   ArrayResize(dBdR, T);
   ArrayResize(dRdF, T+1);
   ArrayResize(dRdFp, T+1);
   ArrayResize(dFdw, M+2);
   ArrayResize(dFpdw, M+2);
   ArrayResize(dSdw, M+2);
   ArrayResize(epoch_S, n_epoch);
   
   ArrayInitialize(x, 0.0);
   ArrayInitialize(F, 0.0);
   ArrayInitialize(R, 0.0);
   ArrayInitialize(w, 1.0);
   ArrayInitialize(w_opt, 1.0);
   ArrayInitialize(sumR, 0.0);
   ArrayInitialize(sumR2, 0.0);
   ArrayInitialize(dBdR, 0.0);
   ArrayInitialize(dRdF, 0.0);
   ArrayInitialize(dRdFp, 0.0);
   ArrayInitialize(dFdw, 0.0);
   ArrayInitialize(dFpdw, 0.0);
   ArrayInitialize(dSdw, 0.0);
   ArrayInitialize(epoch_S, 0.0);
}

// Destructor
TradingRRL::~TradingRRL(){}

// Menber functions
int TradingRRL::quant(double f){
    if(MathAbs(f) < q_threshold){
        return 0;
    }
    else if(f>0.0){
        return 1;
    }
    else{
        return -1;
    }
}

double TradingRRL::sign(double f){
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

double TradingRRL::dot(double &a[], double &b[]){
   double adotb = 0.0;
   for(int i=0; i<ArraySize(a); i++){
      adotb += a[i]*b[i];
   }
   return adotb;
}

double TradingRRL::tanh(double f){
   if(MathAbs(f)>20.0){
      if(f>0.0){
         return 1.0;
      }
      else{
         return -1.0;
      }
   }
   return (MathExp(f)-MathExp(-f))/(MathExp(f)+MathExp(-f));
}

void TradingRRL::set_r(){
   for(int i=0; i<T+M; ++i){
      r[i] = Close[i] - Close[i+1];
   }
}

void TradingRRL::set_x(int t_index, double _Fp){
   x[0] = 1.0;
   x[M+2-1] = _Fp;
   for(int j=1 ; j<M+2-1 ; ++j ){
       x[j] = r[t_index+(j-1)];
   }
}

double TradingRRL::calc_F(int t_index, double _Fp){
   set_x(t_index, _Fp);
   return tanh(dot(w,x));
}

void TradingRRL::set_F(){
    for(int i=T-1 ; i>=0 ; --i ){
        F[i] = calc_F(i, F[i+1]);
    }
}

void TradingRRL::calc_R(){
    int i;
    for(i=T-1 ; i>=0 ; --i ){
        R[i]  = mu*(F[i+1]*r[i] - sigma*MathAbs(F[i] - F[i+1]));
    }
}

void TradingRRL::calc_sumR(){
    sumR[T-1]  = R[T-1];
    sumR2[T-1] = R[T-1]*R[T-1];
    for(int i=T-2 ; i>=0 ; --i ){
        sumR[i]  = sumR[i+1]  + R[i];
        sumR2[i] = sumR2[i+1] + R[i]*R[i];
    }
}

void TradingRRL::calc_dSdw(){
    set_F();
    calc_R();
    calc_sumR();
    A    = sumR[0]/T;
    B    = sumR2[0]/T;
    S    = A/sqrt(B-A*A);
    dSdA = S*(1+S*S)/A;
    dSdB = -S*S*S/2/A/A;
    dAdR = 1.0/T;   
    ArrayInitialize(dFdw, 0.0);
    ArrayInitialize(dFpdw, 0.0);
    ArrayInitialize(dSdw, 0.0);
    for(int i=T-1 ; i>=0 ; --i ){
        dBdR[i]  = 2.0/T*R[i];
        dRdF[i]  = - mu*sigma*sign(F[i] - F[i+1]);
        dRdFp[i] =   mu*r[i] + mu*sigma*sign(F[i] - F[i+1]);
        set_x(i, F[i+1]);
        for(int j=0 ; j<M+2; ++j ){
            if(i!=T-1) dFpdw[j] = dFdw[j];
            dFdw[j]  = (1.0 - F[i]*F[i])*(x[j] + w[M+2-1]*dFpdw[j]);
            dSdw[j] += (dSdA*dAdR + dSdB*dBdR[i])*(dRdF[i]*dFdw[j] + dRdFp[i]*dFpdw[j]);
        }
    }
}

void TradingRRL::update_w(){
    for(int j=0 ; j<M+2; ++j ){
        w[j] += rho * dSdw[j];
    }
}

void TradingRRL::fit(){
    double tic, toc;
    int e_index;

    ArrayInitialize(w, 1.0);
    ArrayInitialize(w_opt, 1.0);
    calc_dSdw();
    Print("Epoch loop start. Initial sharp's ratio is " + DoubleToString(S) +".");
    S_opt = S;
    
    tic = GetTickCount()*1.0e-3;
    for(e_index=0; e_index<n_epoch; ++e_index){
        calc_dSdw();
        if(S > S_opt){
            S_opt = S;
            ArrayCopy(w_opt, w);
        }
        epoch_S[e_index] = S;
        update_w();
        if(e_index%progress_period == progress_period-1){
            toc = GetTickCount()*1.0e-3;
            Print("Epoch: "+ IntegerToString(e_index + 1) + " / " + IntegerToString(n_epoch) + ". Shape's ratio: "+ DoubleToString(S) + ". Elapsed time: " + DoubleToString(toc-tic, 4) + " sec.");
        }
    }
    toc = GetTickCount()*1.0e-3;
    Print("Epoch: "+ IntegerToString(e_index + 1) + " / " + IntegerToString(n_epoch) + ". Shape's ratio: "+ DoubleToString(S) + ". Elapsed time: " + DoubleToString(toc-tic, 4) + " sec.");
    ArrayCopy(w, w_opt);
    calc_dSdw();
    Print("Epoch loop end. Optimized sharp's ratio is " + DoubleToString(S_opt) +".");
}

void TradingRRL::save_weight(){
   int handle;
   handle= FileOpen("w.csv", FILE_WRITE|FILE_CSV);
   if(handle>0){
      for(int i=0; i<ArraySize(w); ++i){
         FileWrite(handle, w[i]);
      }
      FileClose(handle);
   }
   handle= FileOpen("epoch_S.csv", FILE_WRITE|FILE_CSV);
   if(handle>0){
      for(int i=0; i<ArraySize(epoch_S); ++i){
         FileWrite(handle, epoch_S[i]);
      }
      FileClose(handle);
   }
   /*
   handle= FileOpen("r.csv", FILE_WRITE|FILE_CSV);
   if(handle>0){
      for(int i=0; i<ArraySize(r); ++i){
         FileWrite(handle, r[i]);
      }
      FileClose(handle);
   }
   */
}

//--- Create TradingRRL object as global.
TradingRRL *rrl;

//+------------------------------------------------------------------+
//| Expert initialization function                                   |
//+------------------------------------------------------------------+
int OnInit(){
   n_tick    = 0;
   init_flag = True;
   Fp        = 0.0;
   rrl = new TradingRRL();
   return(INIT_SUCCEEDED);
}
//+------------------------------------------------------------------+
//| Expert deinitialization function                                 |
//+------------------------------------------------------------------+
void OnDeinit(const int reason){
   delete rrl;   
}

//+------------------------------------------------------------------+
//| Market functions                                                 |
//+------------------------------------------------------------------+
void closeBuyPos(){  
   int res;
   for(int i=OrdersTotal()-1; i>=0; i--){
      if(OrderSelect(i,SELECT_BY_POS,MODE_TRADES)==false) continue;
      if(OrderMagicNumber()!=MAGIC || OrderSymbol()!=Symbol()) continue;
      if(OrderType()==OP_BUY){  
         res = OrderClose(OrderTicket(),OrderLots(),Bid,slippage,OrangeRed);
         continue;
      }
   }
}

void closeSellPos(){  
   int res;
   for(int i=OrdersTotal()-1; i>=0; i--){
      if(OrderSelect(i,SELECT_BY_POS,MODE_TRADES)==false) continue;
      if(OrderMagicNumber()!=MAGIC || OrderSymbol()!=Symbol()) continue;
      if(OrderType()==OP_SELL){
         res = OrderClose(OrderTicket(),OrderLots(),Ask,slippage,OrangeRed);
         continue;
      }
   }
}

int countBuyPos(){  
   int n_buy_pos = 0;
   for(int i=OrdersTotal()-1; i>=0; i--){
      if(OrderSelect(i,SELECT_BY_POS,MODE_TRADES)==false) continue;
      if(OrderMagicNumber()!=MAGIC || OrderSymbol()!=Symbol()) continue;
      if(OrderType()==OP_BUY){  
         n_buy_pos++;
         continue;
      }
   }
   return n_buy_pos;
}

int countSellPos(){  
   int n_sell_pos = 0;
   for(int i=OrdersTotal()-1; i>=0; i--){
      if(OrderSelect(i,SELECT_BY_POS,MODE_TRADES)==false) continue;
      if(OrderMagicNumber()!=MAGIC || OrderSymbol()!=Symbol()) continue;
      if(OrderType()==OP_SELL){
         n_sell_pos++;
         continue;
      }
   }
   return n_sell_pos;
}

void orderStopLimit(bool _stop_flag, double _stop_points, bool _limit_flag, double _limit_points){  
   int res;
   for(int i=0; i<OrdersTotal(); i++){
      if(OrderSelect(i,SELECT_BY_POS,MODE_TRADES)==false) continue;
      if(OrderMagicNumber()!=MAGIC || OrderSymbol()!=Symbol()) continue;
      //--- For long position.
      if(OrderType()==OP_BUY){  
         if(OrderStopLoss()==0){
            res = OrderModify(OrderTicket(), OrderOpenPrice(),
                              (Bid - (_stop_points  - 2) * Point) * _stop_flag,
                              (Bid + (_limit_points - 2) * Point) * _limit_flag,
                              0, ForestGreen);
            continue;
         }
      }
      //--- For short position.
      if(OrderType()==OP_SELL){
         if(OrderStopLoss()==0){
            res = OrderModify(OrderTicket(), OrderOpenPrice(),
                              (Ask + (_stop_points  - 2) * Point) * _stop_flag,
                              (Ask - (_limit_points - 2) * Point) * _limit_flag,
                              0, ForestGreen);
            continue;
         }
      }
   }
}


//+------------------------------------------------------------------+
//| Expert tick function                                             |
//+------------------------------------------------------------------+
void OnTick(){

   // --- Place stop and limit.
   if(stop_flag || limit_flag){
      orderStopLimit(stop_flag, stop_points, limit_flag, limit_points);
   }
   
   //---Evaluate only the first tick.
   if(Volume[0]>1 || IsTradeAllowed()==false) return;
   
   //--- Wait until the number of bars become T+M.
   n_tick +=1;
   if(init_flag){
      if(Bars <= T+M) return;
      init_flag = False;
      n_tick = 0;
   }

   //--- Training agent.
   if(n_tick % n_tick_update_w == 0){
      Print("Update weights");
      rrl.set_r();
      rrl.fit();
      if(write_log) rrl.save_weight();
   }
   
   //--- The agent decide action with optimized weight.
   double F;
   int qF;
   rrl.set_r();
   F  = rrl.calc_F(0, Fp);
   qF = rrl.quant(F);
   Fp = F;

   //--- Place an order following agent action.
   int res, n_buy_pos, n_sell_pos;
   
   n_buy_pos  = countBuyPos();
   n_sell_pos = countSellPos();
   
   Print("F: " + DoubleToString(F));
   Comment("F: " + DoubleToString(F));
   
   //--- long
   if(qF == 1){ 
      if(n_buy_pos == 0){
         closeSellPos();
         res = OrderSend(Symbol(),OP_BUY,lots,Ask,slippage,0,0,"",MAGIC,0,Blue);
      }
   }
   //--- short
   else if(qF == -1){
      if(n_sell_pos == 0){
         closeBuyPos();
         res = OrderSend(Symbol(),OP_SELL,lots,Bid,slippage,0,0,"",MAGIC,0,Red);
      }
   }
   //--- neutral
   else{
      closeBuyPos();
      closeSellPos();
   }
   
}