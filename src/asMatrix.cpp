#include <Rcpp.h>    
using namespace Rcpp;    

// [[Rcpp::export]]    
NumericMatrix asMatrix(NumericVector rp,
                       NumericVector cp,
                       NumericVector z,
                       int nrows,
                       int ncols){
  
  int k = z.size() ;    
  
  NumericMatrix mat(nrows, ncols);    
  
  for (int i = 0; i < k; i++){    
    mat(rp[i],cp[i]) = z[i];    
  }    
  
  return mat;    
}
