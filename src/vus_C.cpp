#include <Rcpp.h>
using namespace Rcpp;

static inline double indvus(double a, double b, double c) {
  if(a > b || c < b) return 0;
  else{
    double d = 0;
    if(a == b) d += 1;
    if(c == b) d += 1;
    return (8 - 3*d)/(8 + 2*d);
  }
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
double vusC_full(NumericVector tt1, NumericVector tt2, NumericVector tt3){
  int nn1 = tt1.size(), nn2 = tt2.size(), nn3 = tt3.size();
  double num = 0.0;
  for(int i = 0; i < nn1; i++){
    for(int j = 0; j < nn2; j++){
      for(int k = 0; k < nn3; k++){
        num += indvus(tt1[i], tt2[j], tt3[k]);
      }
    }
  }
  double out = num/(nn1*nn2*nn3);
  return out;
}

//' @export
// [[Rcpp::export]]
double vusC_sigmoid(NumericVector tt1, NumericVector tt2, NumericVector tt3, 
                    NumericVector sig_n) {
  int nn1 = tt1.size(), nn2 = tt2.size(), nn3 = tt3.size();
  double num = 0.0;
  for(int i = 0; i < nn1; i++){
    for(int j = 0; j < nn2; j++){
      for(int k = 0; k < nn3; k++){
        num += sigmoid((tt2[j] - tt1[i])/sig_n[0])*sigmoid((tt3[k] - tt2[j])/sig_n[1]);
      }
    }
  }
  double out = num/(nn1*nn2*nn3);
  return out;
}

//' @export
// [[Rcpp::export]]
double vusC_probit(NumericVector tt1, NumericVector tt2, NumericVector tt3, 
                   NumericVector sig_n) {
  int nn1 = tt1.size(), nn2 = tt2.size(), nn3 = tt3.size();
  double num = 0.0;
  for(int i = 0; i < nn1; i++){
    for(int j = 0; j < nn2; j++){
      for(int k = 0; k < nn3; k++){
        num += probit((tt2[j] - tt1[i])/sig_n[0])*probit((tt3[k] - tt2[j])/sig_n[1]);
      }
    }
  }
  double out = num/(nn1*nn2*nn3);
  return out;
}
