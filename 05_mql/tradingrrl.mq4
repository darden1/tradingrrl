//+------------------------------------------------------------------+
//|                                                   tradingrrl.mq4 |
//|                        Copyright 2014, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "Copyright 2014, MetaQuotes Software Corp."
#property link      "http://www.mql5.com"
#property version   "1.00"
#property strict

#define MAGIC 20170214

extern int    slippage = 2;
extern double lots = 0.1;
extern bool stop_flag = False;
extern bool limit_flag = False;
extern double stop_points = 1000;
extern double limit_points = 1000;

extern int T = 150;
extern int M = 30;
extern double q_threshold = 0.7;
extern double mu = 10000;
extern double sigma = 0.04;
extern double alpha = 2.0;
extern int n_epoch = 10000;
extern int n_tick_update_w = 150;

extern bool write_log = True;

int n_tick;
bool init_flag;
double pre_F;
double r[], w[];

//+------------------------------------------------------------------+
//| Expert initialization function                                   |
//+------------------------------------------------------------------+
int OnInit()
  {
//---
   n_tick = 0;
   init_flag =True;
   pre_F = 0.0;
   ArrayResize(r, T+M);
   ArrayResize(w, M+2);
   init_w(w);
   
//---
   return(INIT_SUCCEEDED);
  }
//+------------------------------------------------------------------+
//| Expert deinitialization function                                 |
//+------------------------------------------------------------------+
void OnDeinit(const int reason)
  {
//---
   
  }
//+------------------------------------------------------------------+
//| My functions                                                     |
//+------------------------------------------------------------------+
void closeBuyPos(){  
   int res;
   for(int i=OrdersTotal()-1; i>=0; i--)
   {
      if(OrderSelect(i,SELECT_BY_POS,MODE_TRADES)==false)continue;
      if(OrderMagicNumber()!=MAGIC || OrderSymbol()!=Symbol())continue;
      if(OrderType()==OP_BUY)
         {  
               res = OrderClose(OrderTicket(),OrderLots(),Bid,slippage,OrangeRed);
               continue;
         }
   }
}

void closeSellPos(){  
   int res;
   for(int i=OrdersTotal()-1; i>=0; i--)
   {
      if(OrderSelect(i,SELECT_BY_POS,MODE_TRADES)==false)continue;
      if(OrderMagicNumber()!=MAGIC || OrderSymbol()!=Symbol())continue;
      if(OrderType()==OP_SELL)
         {
               res = OrderClose(OrderTicket(),OrderLots(),Ask,slippage,OrangeRed);
               continue;
         }
   }
}

int countBuyPos(){  
   int n_buy_pos = 0;
   
   for(int i=OrdersTotal()-1; i>=0; i--)
   {
      if(OrderSelect(i,SELECT_BY_POS,MODE_TRADES)==false)continue;
      if(OrderMagicNumber()!=MAGIC || OrderSymbol()!=Symbol())continue;
      if(OrderType()==OP_BUY)
         {  
               n_buy_pos++;
               continue;
         }
   }
   return n_buy_pos;
}

int countSellPos(){  
   int n_sell_pos = 0;
   for(int i=OrdersTotal()-1; i>=0; i--)
   {
      if(OrderSelect(i,SELECT_BY_POS,MODE_TRADES)==false)continue;
      if(OrderMagicNumber()!=MAGIC || OrderSymbol()!=Symbol())continue;
      if(OrderType()==OP_SELL)
         {
               n_sell_pos++;
               continue;
         }
   }
   return n_sell_pos;
}


void orderStopLimit(bool stop_flag, double stop_points, bool limit_flag, double limit_points)
{  
   int res;
   for(int i=0; i<OrdersTotal(); i++)
   {
      if(OrderSelect(i,SELECT_BY_POS,MODE_TRADES)==false)continue;
      if(OrderMagicNumber()!=MAGIC || OrderSymbol()!=Symbol())continue;
      //---Position is buy
      if(OrderType()==OP_BUY)
      {  
         if(OrderStopLoss()==0)
         {
            res = OrderModify(OrderTicket(), OrderOpenPrice(),
                              (Bid-(stop_points-2)*Point)*stop_flag,
                              (Bid+(limit_points-2)*Point)*limit_flag,
                              0, ForestGreen);
            continue;
         }

      }
      //---Position is sell
      if(OrderType()==OP_SELL)
      {
         if(OrderStopLoss()==0)
         {
            res = OrderModify(OrderTicket(), OrderOpenPrice(),
                              (Ask+(stop_points-2)*Point)*stop_flag,
                              (Ask-(limit_points-2)*Point)*limit_flag,
                              0, ForestGreen);
            continue;
         }
      }
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

int quantizer(double f){
   if(MathAbs(f) < q_threshold){
      return 0;
   }
   else if(f > 0.0)   {
      return 1;
   }
   else   {
      return -1;
   }
}

double dot(double &w[], double &x[]){
   double wdotx = 0.0;
   for(int i=0; i<M+2; i++)
   {
      wdotx += w[i]*x[i];
      //Print("w:" + w[i] +"/x:" + x[i]);
   }
   
   return wdotx;
}

double tanh(double f){
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

void set_r(double &r[]){
   for(int i=0; i<T+M; i++)
   {
      r[i] = Close[i] - Close[i+1];
   }
}

void init_w(double &w[]){
   for(int i=0; i<M+2; i++)
   {
      w[i] = 1.0;
   }
}

void set_x(int t_index, double &x[], double &r[], double l_pre_F){
   x[0] = 1.0;
   x[M+2-1] = l_pre_F;
   for( int j=1 ; j<M+2-1 ; ++j ){
       x[j] = r[t_index+(j-1)];
   }

}


double get_F(int t_index, double &r[], double &w[], double l_pre_F){
   double x[];
   double wdotx, f;
   ArrayResize(x, M+2);
   set_x(t_index, x, r, l_pre_F);
   wdotx = dot(w,x);
   f = tanh(wdotx);
   return f;
}

double update_w(double &r[], double &w[]){
   double wdotx, A, B, S, dSdA, dSdB, dAdR;
   double x[], F[], R[], sum_R[], sum_R2[], dBdR[], dRdF[], dRdFm[], pre_dFdw[], dFdw[], dSdw[];
   ArrayResize(x, M+2);
   ArrayResize(F, T+1);
   ArrayResize(R, T);
   ArrayResize(sum_R, T+1);
   ArrayResize(sum_R2, T+1);
   ArrayResize(dBdR, T);
   ArrayResize(dRdF, T+1);
   ArrayResize(dRdFm, T+1);
   ArrayResize(dFdw, M+2);
   ArrayResize(pre_dFdw, M+2);
   ArrayResize(dSdw, M+2);

   sum_R[T] =0.0;
   sum_R2[T] =0.0;
   F[T] = 0.0;
   
   for( int i=T-1 ; i>=0 ; --i ){
      F[i] = get_F(i, r, w, F[i+1]);
      R[i]  = mu*(F[i+1]*r[i] - sigma*MathAbs(F[i] - F[i+1]));
      sum_R[i] = sum_R[i+1] + R[i];
      sum_R2[i] = sum_R2[i+1] + R[i]*R[i];
      dBdR[i] = 2.0/T*R[i];
      dRdF[i]    = - mu*sigma*sign(F[i] - F[i+1]);
      dRdFm[i]   =   mu*r[i] + mu*sigma*sign(F[i] - F[i+1]);
   }
        
   A = sum_R[0]/T;
   B = sum_R2[0]/T;
   S = A/sqrt(B-A*A);
   
   dSdA = (1+S*S)/sqrt(B-A*A);
   dSdB = -S*S/2/A/sqrt(B-A*A);
   dAdR = 1.0/T;   
   
   
   for(int i=T-1 ; i>=0 ; --i ){
      if(i==T-1){
          for( int j=0 ; j<M+2; ++j ){
              pre_dFdw[j] = 0.0;
              dSdw[j] = 0.0;
          }
      }
      else{
          for( int j=0 ; j<M+2; ++j ){
              pre_dFdw[j] = dFdw[j];
          }
      }
      set_x(i, x, r, F[i+1]);
      for( int j=0 ; j<M+2; ++j ){
          dFdw[j]=(1 - F[i]*F[i])*(x[j] + w[M+2-1]*pre_dFdw[j]);
          dSdw[j]+=(dSdA*dAdR + dSdB*dBdR[i])*(dRdF[i]*dFdw[j] + dRdFm[i]*pre_dFdw[j]);
      }
   }
   
   // --- Update weight w
   for( int j=0 ; j<M+2; ++j ){
      w[j]+=alpha*dSdw[j];
   }
   return S;
 
}

int file_output(string fname, double &array[]){
   int handle= FileOpen(fname, FILE_WRITE|FILE_CSV);
   if(handle>0)
   {
      for(int i=0; i<ArraySize(array); ++i)
      FileWrite(handle, array[i]);
      FileClose(handle);
   }
   return handle;
}
//+------------------------------------------------------------------+
//| Expert tick function                                             |
//+------------------------------------------------------------------+
void OnTick()
  {
//---
   // --- stop and limit
   if(stop_flag || limit_flag){
      orderStopLimit(stop_flag, stop_points, limit_flag, limit_points);
   }
   
   if(Volume[0]>1 || IsTradeAllowed()==false) return;
   
   n_tick +=1;
   if(init_flag){
      if(Bars <= T+M){
         //---Create Bars;
         //Print("CloseSize:" + ArraySize(Close));
         //Print("Bars:" +Bars);
         return;
      }
      init_flag = False;
      n_tick = 0;
   }
   
   
  
   //------ 
   int res, qF, n_buy_pos, n_sell_pos;
   double F;
   set_r(r);
   
   //---
   if(n_tick%n_tick_update_w==0){
      Print("Update weights");
      double S;
      double epoch_S[];
      ArrayResize(epoch_S, n_epoch);
      
      init_w(w);
      for(int e_index=1 ; e_index<=n_epoch; ++e_index){
         S = update_w(r, w);
         epoch_S[e_index-1] = S;
         if(e_index%100==0){
            Print("Epoch: " + IntegerToString(e_index) + "/" + IntegerToString(n_epoch));
         }
      }
      if(write_log){
         file_output("w_opt.csv", w);
         file_output("r.csv", r);
         file_output("epoch_S.csv", epoch_S);
      }
   }
   
   
   F = get_F(0, r, w, pre_F);
   qF = quantizer(F);
   pre_F = F;
   
   n_buy_pos  = countBuyPos();
   n_sell_pos = countSellPos();

   if(qF == 1)
   {
      if(n_buy_pos == 0)
      {
         closeSellPos();
         res = OrderSend(Symbol(),OP_BUY,lots,Ask,slippage,0,0,"",MAGIC,0,Blue);
      }
   }
   else if(qF == -1)
   {
      if(n_sell_pos == 0)
      {
         closeBuyPos();
         res = OrderSend(Symbol(),OP_SELL,lots,Bid,slippage,0,0,"",MAGIC,0,Red);
      }
   }
   else
   {
      closeBuyPos();
      closeSellPos();
   }
  }
//+------------------------------------------------------------------+
