#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#define pi 3.141592653589793238462643383280

// Note: -----------------------------------------------------------------------
// * Cluster index starts with 0.
// * Update at-risk indicator (done)
// * Update important variable indicator (done)
// * Update beta (wip)
// * Update the cluster assignment
// -----------------------------------------------------------------------------

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
double log_wj(int j, arma::mat z, arma::mat gamma_mat, arma::vec w, 
              arma::mat beta, arma::uvec clus_assign, double b0w, double b1w){
  
  // Convert beta to xi by taking the exponential function.
  arma::mat xi = arma::exp(beta);
  double result = 0.0;
  
  // Get the value of wj
  int wj = w[j];
  
  // Calculate the gamma_w_xi factor
  arma::mat w_mat = arma::repelem(w, 1, z.n_rows); 
  w_mat = w_mat.t();
  arma::mat xi_obs = xi.rows(clus_assign);
  arma::mat gwx = gamma_mat % w_mat % xi_obs;
  arma::mat gwx_j = gwx.col(j);
  arma::vec z_gwx_j = gwx_j + z.col(j);
  
  result += R::lbeta(b0w + wj, b1w + (1 - wj));
  result += arma::accu(arma::lgamma(z_gwx_j.elem(arma::find(z_gwx_j != 0))));
  result -= arma::accu(arma::lgamma(gwx_j.elem(arma::find(gwx_j != 0))));
  result += arma::accu(arma::lgamma(arma::sum(gwx, 1)));
  result -= arma::accu(arma::lgamma(arma::sum(z, 1) + arma::sum(gwx, 1)));

  return result;
  
} 

// [[Rcpp::export]]
arma::vec log_beta_k(arma::vec beta_k, int ci, arma::mat z, arma::mat gamma_mat, 
                     arma::mat beta, arma::vec w, arma::uvec clus_assign, double s2){
  
  // Filter only the active variables
  arma::uvec active_var = arma::find(w == 1);
  arma::mat z_active = z.cols(active_var);
  arma::mat gamma_active = gamma_mat.cols(active_var);
  arma::mat xi_active = arma::exp(beta.cols(active_var));
  arma::vec beta_k_active = beta_k.rows(active_var);
  
  // Filter only the observations in the same cluster (ci)
  arma::uvec clus_index = arma::find(clus_assign == ci);
  z_active = z_active.rows(clus_index);
  gamma_active = gamma_active.rows(clus_index);
  xi_active = xi_active.row(ci);
  
  // Calculate the gamma_w_xi factor
  xi_active = arma::repelem(xi_active.t(), 1, z_active.n_rows);
  xi_active = xi_active.t();
  arma::mat gwx = gamma_active % xi_active;
  arma::mat z_gwx = z_active + gwx;
  
  arma::vec result(active_var.size(), arma::fill::value(0.0));
  Rcpp::NumericVector beta_rcpp = Rcpp::NumericVector(beta_k_active.begin(), beta_k_active.end());
  result += Rcpp::as<arma::vec>(Rcpp::wrap(Rcpp::dnorm4(beta_rcpp, 0, std::sqrt(s2), 1)));
  result += arma::vec(active_var.size(), arma::fill::value(arma::accu(arma::lgamma(arma::sum(gwx, 1)))));
  result -= arma::vec(active_var.size(), arma::fill::value(arma::accu(arma::lgamma(arma::sum(z_gwx, 1)))));
  
  // Replace 0 in the gmx with 1 as lgamma(1) = 0
  result += arma::sum(arma::lgamma(z_gwx.replace(0.0, 1.0)), 0).t();
  result -= arma::sum(arma::lgamma(gwx.replace(0.0, 1.0)), 0).t();
  
  return result;
  
}

// [[Rcpp::export]]
double log_w(arma::mat z, arma::mat g, arma::vec w, arma::mat xi,
               arma::uvec clus_assign, double b0w, double b1w){
  
  /*
   * Description: This is a function for calculating probability of the important
   *              variable indicator for the variable j, p(wj|.), in a log scale.
   * Input: The index of the current variable (j), the dataset (z), 
   *        the at-risk indicator for all observation (g), the important 
   *        variable indicator (w), the xi for all cluster (xi), the 
   *        hyperparameter of the important variable indicator (b0w, b1w).
   */
  
  double result = 0.0;
  
  arma::vec b0 = b0w + w;
  arma::vec b1 = b1w + (1 - w);
  
  result += arma::accu(arma::lgamma(b0) + arma::lgamma(b1) - arma::lgamma(b0 + b1));
  
  // Select xi for each observation based on its cluster
  arma::mat xi_obs = xi.rows(clus_assign);
  
  // Duplicate vector w
  arma::mat w_mat = arma::repelem(w, 1, z.n_rows); 
  w_mat = w_mat.t();
  
  // Calculate gamma_w_xi and z_gamma_w_xi for all observations
  arma::mat gwx = g % w_mat % xi_obs;
  arma::mat z_gwx = z + gwx;
  
  result += arma::accu(arma::lgamma(arma::sum(gwx, 1)));
  result -= arma::accu(arma::lgamma(gwx.elem(arma::find(gwx != 0))));
  result += arma::accu(arma::lgamma(z_gwx.elem(arma::find(z_gwx != 0))));
  result -= arma::accu(arma::lgamma(arma::sum(z_gwx, 1)));
  
  return result;
  
}

// [[Rcpp::export]]
arma::mat update_gamma(arma::mat z, arma::uvec clus_assign, arma::vec w, 
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
    arma::vec xi_i = arma::conv_to<arma::vec>::from(xi.row(clus_assign[i]));
    
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

// [[Rcpp::export]]
arma::vec update_w(arma::mat z, arma::uvec clus_assign, arma::vec old_w,
                   arma::mat gamma_mat, arma::mat beta, double b0w, double b1w){
  
  arma::vec proposed_w(old_w);
  
  // Get the index of the important and not important variables
  arma::uvec imp_index = arma::find(proposed_w == 1);
  arma::uvec n_imp_index = arma::find(proposed_w == 0);
  int adjusted_index = -1;
  int j_update = -1;
  
  // Propose a new vector of the important variable indicator
  if(imp_index.size() <= 1){ // If only one variable is active, it cannot delete that one out.
    adjusted_index = arma::randperm(n_imp_index.size(), 1)[0];
    j_update = n_imp_index[adjusted_index];
    proposed_w.row(j_update).fill(1);
  } else {
    adjusted_index = arma::randperm(old_w.size(), 1)[0];
    j_update = adjusted_index;
    proposed_w.row(j_update).fill(!(old_w[j_update]));
  }
  
  // Calculate the acceptance probability
  double logA = log_wj(j_update, z, gamma_mat, proposed_w, beta, clus_assign, b0w, b1w) - 
    log_wj(j_update, z, gamma_mat, old_w, beta, clus_assign, b0w, b1w);
  double logU = std::log(R::runif(0.0, 1.0));
  arma::vec new_w(old_w);
  if(logU <= logA){
    new_w = proposed_w;
  }
  
  return new_w;
  
}

// [[Rcpp::export]]
arma::mat update_w_test(arma::mat z, arma::uvec clus_assign, arma::vec old_w, 
                   arma::mat gamma_mat, arma::mat xi, double b0w, double b1w,
                   int iter_w, double trh){
  
  arma::mat result_w(iter_w, old_w.size(), arma::fill::value(-1)); 
  arma::vec final_w(5, arma::fill::value(-1));
  
  arma::vec previous_w(old_w);
  for(int t = 0; t < iter_w; ++t){
    
    arma::vec proposed_w(old_w); 
    
    // Get the index of the important and not important variables
    arma::uvec imp_index = arma::find(previous_w == 1);
    arma::uvec n_imp_index = arma::find(previous_w == 0);
    
    // Decide to perform add/delete or swap
    if(n_imp_index.size() == 0 or imp_index.size() == 0){
      arma::uvec adj_var_index = arma::randperm(previous_w.size(), 1);
      proposed_w.row(adj_var_index[0]).fill(!previous_w[adj_var_index[0]]);
    } else {
      double u = R::runif(0, 1);
      if(u <= 0.5){ // Add/Delete
        arma::uvec adj_var_index = arma::randperm(previous_w.size(), 1);
        proposed_w.row(adj_var_index[0]).fill(!previous_w[adj_var_index[0]]);
      } else { // Swap
        arma::uvec imp_i = arma::randperm(imp_index.size(), 1);
        arma::uvec n_imp_i = arma::randperm(n_imp_index.size(), 1);
        proposed_w.row(imp_index[imp_i[0]]).fill(0);
        proposed_w.row(n_imp_index[n_imp_i[0]]).fill(1);
      }
    }
    
    // Calculate the acceptance probability
    double logA = log_w(z, gamma_mat, proposed_w, xi, clus_assign, b0w, b1w) - 
      log_w(z, gamma_mat, previous_w, xi, clus_assign, b0w, b1w);
    double logU = std::log(R::runif(0, 1));
    if(logU <= logA){
      previous_w = proposed_w;
    }
    
    result_w.row(t) = proposed_w.t();
    
  }
  
  
  
  return result_w;
}

