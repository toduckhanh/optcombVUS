#include <Rcpp.h>
using namespace Rcpp;

inline double ind_auc(double a, double b){
  if(a < b) return 1;
  else if(a == b) return 0.5;
  else return 0;
}

static inline double sigmoid(double x){
  double y = 1/(1 + exp(-x));
  return y;
}

static inline double probit(double x){
  NumericVector xs(1, x);
  double y = as<double>(wrap(pnorm(xs, 0.0, 1.0)));
  return y;
}

//' @export
// [[Rcpp::export]]
double full_AUC(NumericVector tt_0, NumericVector tt_1) {
  int n0 = tt_0.size();
  int n1 = tt_1.size();
  double res = 0.0;
  for(int i = 0; i < n0; ++i){
    for(int j = 0; j < n1; ++j){
      res += ind_auc(tt_0[i], tt_1[j]);
    }
  }
  return res/n0/n1;
}

//' @export
// [[Rcpp::export]]
double auc_C_smooth(NumericVector tt1, NumericVector tt2, double h, char link){
  int n1 = tt1.size(), n2 = tt2.size();
  double num = 0.0;
  switch(link){
  case 's' :
    for(int i = 0; i < n1; ++i){
      for(int j = 0; j < n2; ++j){
        num += sigmoid((tt2[j] - tt1[i])/h);
      }
    }
    break;
  case 'p':
    for(int i = 0; i < n1; ++i){
      for(int j = 0; j < n2; ++j){
        num += probit((tt2[j] - tt1[i])/h);
      }
    }
    break;
  }
  return num/n1/n2;
}


