#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector lcrec(NumericVector y,
                                 NumericVector logsigma2,
                                 double alpha0,
                                 double alpha1,
                                 double beta1) {
  NumericVector y2 = Rcpp::pow(y,2);
  int n = y2.size();
  for(int i = 1; i < n ; ++i) {
    
    logsigma2[i] = alpha0 + alpha1 * log(y2[i - 1]) + beta1 * logsigma2[i - 1];
  }
  return logsigma2;
  
}
// [[Rcpp::export]]  
DataFrame derivreclm(NumericVector leta2,
                                    NumericVector logsigma2,
                                    NumericVector beta1,
                                    NumericVector lm2) {
    int n = logsigma2.size();
    NumericVector betavector (n);
    
    NumericVector dalpha0x (n);
    NumericVector dalpha1x (n);
    NumericVector dbeta1x (n);
    NumericVector dlm2x (n);
    
    for(int i = 0; i < n ; ++i) {
      
      betavector[i] = Rcpp::pow(beta1, i)[1];
      
      
      
      
      
      
      dalpha0x[i] = Rcpp::sum(betavector);
      dalpha1x[i] = Rcpp::sum(Rcpp::rev(betavector[Rcpp::Range(0,i)]) * leta2[Rcpp::Range(0,i)]);
      dbeta1x[i] = Rcpp::sum(Rcpp::rev(betavector[Rcpp::Range(0,i)]) * logsigma2[Rcpp::Range(0,i)]);
      dlm2x[i] = Rcpp::sum(Rcpp::rev(betavector[Rcpp::Range(0,i)]) * lm2[Rcpp::Range(0,i)]);                 
      
    }
    
    DataFrame df = DataFrame::create(Named("da0") = dalpha0x,
                                     Named("da1") = dalpha1x,
                                     Named("db1") = dbeta1x,
                                     Named("dlm") = dlm2x);
    return df;
    
  }
// [[Rcpp::export]]  
DataFrame derivreclev(NumericVector leta2,
                                     NumericVector logsigma2,
                                     NumericVector beta1) {
    int n = logsigma2.size();
    NumericVector betavector (n);
    
    NumericVector dalpha0x (n);
    NumericVector dalpha1x (n);
    NumericVector dbeta1x (n);
    
    for(int i = 0; i < n ; ++i) {
      
      betavector[i] = Rcpp::pow(beta1, i)[1];
      
      
      
      
      
      
      dalpha0x[i] = Rcpp::sum(betavector);
      dalpha1x[i] = Rcpp::sum(Rcpp::rev(betavector[Rcpp::Range(0,i)]) * leta2[Rcpp::Range(0,i)]);
      dbeta1x[i] = Rcpp::sum(Rcpp::rev(betavector[Rcpp::Range(0,i)]) * logsigma2[Rcpp::Range(0,i)]);
      
      
    }
    
    DataFrame df = DataFrame::create(Named("da0") = dalpha0x,
                                     Named("da1") = dalpha1x,
                                     Named("db1") = dbeta1x);
    return df;
    
  }
  
// [[Rcpp::export]]    
NumericVector cstoch(NumericVector h,
                                    NumericVector eta,
                                    double sbeta,
                                    double sigma,
                                    double mu,
                                    double l) {
    
    for(int i = 1; i < l  ; ++i) {
      
      h[i] = mu + sbeta * h[i - 1] + sigma * eta[i];
    }
    
    return h;
    
  }
  