#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#define pi 3.141592653589793238462643383280

// Note: -----------------------------------------------------------------------
// *

// -----------------------------------------------------------------------------

// User-defined function: ------------------------------------------------------
// [[Rcpp::export]]
double log_g_ijk(int j, arma::vec zi, arma::vec gi, arma::vec w, arma::vec xi_k,
                 double b0g, double b1g){
  
  /*
   * Description: This is a function for calculating probability of the at-risk 
   *              indicator for the particular observations, p(gamma_ijk|.), in 
   *              a log scale.
   * Input: The index of the current variable (j), the observation i (zi), 
   *        the at-risk indicator for this observation (gi), the important 
   *        variable indicator (w), the xi that corresponding to the cluster 
   *        (xi_k), the hyperparameter of the at-risk indicator (b0g, b1g).
   *        
   */
  
  double result = 0.0;
  int g_ijk = gi[j];
  arma::vec xi_gw = gi % w % xi_k;
  arma::vec zi_xi = zi + xi_gw;
  
  result += std::log(R::beta(b0g + g_ijk, b1g + (1 - g_ijk)));
  result -= std::log(R::beta(b0g, b1g));
  result += std::lgamma(arma::accu(xi_gw));
  result -= arma::accu(arma::lgamma(xi_gw.rows(arma::find(xi_gw != 0))));
  result += arma::accu(arma::lgamma(zi_xi.rows(arma::find(zi_xi != 0))));
  result -= std::lgamma(arma::accu(zi + xi_gw));

  return result;
  
}

// [[Rcpp::export]]
arma::mat update_gamma(arma::mat z, arma::vec clus_assign, arma::vec w, 
                       arma::mat old_gamma, arma::mat xi, double b0g, double b1g){
  
  /*
   * Description: This is a function for calculating probability of the at-risk 
   *              indicator for the particular observations, p(gamma_ijk|.), in 
   *              a log scale.
   * Input: The index of the current variable (j), the observation i (zi), 
   *        the at-risk indicator for this observation (gi), the important 
   *        variable indicator (w), the xi that corresponding to the cluster 
   *        (xi_k), the hyperparameter of the at-risk indicator (b0g, b1g).
   *        
   */
  
  // Matrix for storing the final result
  arma::mat new_gamma(old_gamma);
  
  // Update the gamma for the only active variables
  arma::uvec active_var = arma::find(w == 1);
  
  // Loop through the observation to update the at-risk indicator
  for(int i = 0; i < z.n_rows; ++i){
    
    arma::vec zi = arma::conv_to<arma::vec>::from(z.row(i));
    arma::vec xi_i = arma::conv_to<arma::vec>::from(xi.row(clus_assign[i] - 1));
    
    // Consider only the active variables only
    for(int j = 0; j < active_var.size(); ++j){
      
      int wj = active_var[j];
      
      // We will update gamma for z(i, wj) = 0 only as the at-risk indicator 
      // for z(i, wj) = 0 can be both 0 and 1, while the at-risk indicator for 
      // z(i, wj) > 0 can be one only.
      if(z(i, wj) == 0){ 
        arma::vec old_gi = arma::conv_to<arma::vec>::from(new_gamma.row(i));
        arma::vec proposed_gi(old_gi);
        
        int g_ijk = old_gamma(i, wj);
        if(g_ijk == 0){
          proposed_gi.row(wj).fill(1);
        } else {
          proposed_gi.row(wj).fill(0);
        }
        
        // Calculate the log acceptance probability
        double logA = log_g_ijk(wj, zi, proposed_gi, w, xi_i, b0g, b1g) - 
          log_g_ijk(wj, zi, old_gi, w, xi_i, b0g, b1g);
        
        // Update g_ijk
        double logU = std::log(R::runif(0.0, 1.0));
        if(logU < logA){
          new_gamma.row(i) = proposed_gi.t();
        }
        
      }
      
    }
    
  }
  
  return new_gamma;
}