#include <RcppArmadillo.h>
#include <unordered_map>
#include <vector>
#include <string>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat MergeRowsForDense(const arma::mat& matrix, const std::vector<std::string>& groupings) {
  
  int nrow = matrix.n_rows;
  int ncol = matrix.n_cols;
  std::unordered_set<std::string> unique_groupings(groupings.begin(), groupings.end());
  int n = unique_groupings.size();
  
  if (nrow != static_cast<int>(groupings.size())) stop("groupings must be the length of nrow(matrix)");
  
  std::unordered_map<std::string, int> rowIndices(n); //reserve space when declaring
  
  arma::mat result(n, ncol, arma::fill::zeros); // reserve space for result matrix
  
  int index = 0;
  for (int i = 0; i < nrow; ++i) {
    const std::string& row = groupings[i]; // avoid copying string
    auto it = rowIndices.find(row);
    if (it != rowIndices.end()) {
      result.row(it->second) += matrix.row(i);
    } else {
      rowIndices.emplace(row, index);
      result.row(index) = matrix.row(i);
      ++index;
    }
  }
  
  return result;
}

// [[Rcpp::export]]
arma::sp_mat MergeRowsForSparse(const arma::sp_mat& matrix, const std::vector<std::string>& groupings) {
  // Create a mapping matrix
  std::unordered_map<std::string, int> unique_groupings;
  for (const auto &group : groupings) {
    if (unique_groupings.find(group) == unique_groupings.end()) {
      unique_groupings[group] = unique_groupings.size();
    }
  }
  arma::sp_mat mapping(groupings.size(), unique_groupings.size());
  for (size_t i = 0; i < groupings.size(); ++i) {
    mapping(i, unique_groupings[groupings[i]]) = 1;
  }
  
  // Multiply mapping by matrix
  arma::sp_mat result = mapping.t() * matrix;
  
  return result;
}

// [[Rcpp::export]]
SEXP MergeRowsC(const SEXP& matrix, const std::vector<std::string>& groupings) {
  if(Rf_inherits(matrix, "dgCMatrix")){    
    arma::sp_mat mat = MergeRowsForSparse(as<arma::sp_mat>(S4(matrix)), groupings);  
    return(wrap(mat));
  } else if (Rf_isMatrix(matrix)){
    arma::mat mat =  MergeRowsForDense(as<arma::mat>(matrix), groupings);    
    return(wrap(mat));  
  } else {    
    stop("Unsupported matrix type");    
  }  
}  
