#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#define pi 3.141592653589793238462643383280

// [[Rcpp::export]]
arma::vec rmultinom_1(arma::vec &probs_arma){
  
  /* Description: sample from the multinomial(N, probs).
   * Credit: https://gallery.rcpp.org/articles/recreating-rmultinom-and-rpois-with-rcpp/
   */
  
  Rcpp::NumericVector probs = Rcpp::NumericVector(probs_arma.begin(), probs_arma.end());
  unsigned int N = probs_arma.size();
  Rcpp::IntegerVector outcome(N);
  rmultinom(1, probs.begin(), N, outcome.begin());
  arma::vec mul_result = Rcpp::as<arma::vec>(Rcpp::wrap(outcome));

  return mul_result;
  
}

// [[Rcpp::export]]
arma::vec log_sum_exp(arma::vec log_unnorm_prob){
  
  /* Description: This function will calculate the normalized probability 
   *              by applying log-sum-exp trick on the log-scale probability.
   * Credit: https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
   */
  
  double max_elem = log_unnorm_prob.max();
  double t = log(0.00000000000000000001) - log(log_unnorm_prob.size());          
  
  for(int k = 0; k < log_unnorm_prob.size(); ++k){
    double prob_k = log_unnorm_prob.at(k) - max_elem;
    if(prob_k > t){
      log_unnorm_prob.row(k).fill(std::exp(prob_k));
    } else {
      log_unnorm_prob.row(k).fill(0.00000000000000000001);
    }
  } 
  
  // Normalize the vector
  return log_unnorm_prob/arma::accu(log_unnorm_prob);
}

// [[Rcpp::export]]
arma::vec logmar_k(arma::mat z, arma::mat atrisk, arma::rowvec beta_k){
  
  /* Calculate the log marginal for each observation in a log scale, given 
     that they are all in cluster k */
  
  arma::vec lmar(z.n_rows, arma::fill::zeros);
  arma::rowvec e_beta = arma::exp(beta_k);
  
  // At-risk indicator always equal 1 for non-zero z_ijk. It can be 0 if z_ijk = 0 only.
  // Data (z)
  lmar += arma::lgamma(arma::sum(z, 1) + 1);
  lmar -= arma::sum(arma::lgamma(z + 1), 1);
  
  // Hyperparameter part (beta)
  for(int i = 0; i < z.n_rows; ++i){
    arma::uvec atrisk_zero = arma::find(atrisk.row(i) == 1);
    arma::rowvec zi_ebeta = z.row(i) + e_beta;
    lmar[i] += arma::accu(arma::lgamma(zi_ebeta.cols(atrisk_zero)));
    lmar[i] -= std::lgamma(arma::accu(zi_ebeta.cols(atrisk_zero)));
  }
  
  return lmar;
  
}

// [[Rcpp::export]]
arma::mat logmar(arma::mat z, arma::mat atrisk, arma::mat beta_mat){
  
  /* Calculate the log marginal for all observations in all possible clusters */
  unsigned int Kmax = beta_mat.n_rows;
  arma::mat logmar_mat(z.n_rows, Kmax, arma::fill::zeros);
  
  // At-risk indicator always equal 1 for non-zero z_ijk. It can be 0 if z_ijk = 0 only.
  // Data (z)
  arma::vec sum_zi = arma::sum(z, 1); 
  logmar_mat += arma::repelem(arma::lgamma(sum_zi + 1), 1, Kmax);
  logmar_mat -= arma::repelem(arma::sum(arma::lgamma(z + 1), 1), 1, Kmax);
  
  // Cluster (Beta) -- Loop through every cluster for changing the beta
  for(int k = 0; k < Kmax; ++k){
    arma::rowvec e_beta_k = arma::exp(beta_mat.row(k));
    arma::vec marginal_part(z.n_rows, arma::fill::zeros);
    
    for(int i = 0; i < z.n_rows; ++i){
      arma::uvec atrisk_loc = arma::find(atrisk.row(i) == 1);
      arma::rowvec e_beta_k_zi = z.row(i) + e_beta_k;
      marginal_part[i] += arma::accu(arma::lgamma(e_beta_k_zi.cols(atrisk_loc)));
      marginal_part[i] -= std::lgamma(arma::accu(e_beta_k_zi.cols(atrisk_loc))); 
    }
    
    logmar_mat.col(k) += marginal_part;
    
  }
  
  return logmar_mat;
    
}

// Algorithms for updating parameters: -----------------------------------------

// At-risk

// Beta
// [[Rcpp::export]]
arma::mat update_beta(arma::mat z, arma::mat atrisk, arma::mat beta_old, 
                      arma::uvec ci, double mu, double s2, double s2_MH){
  
  /* Update the beta vector parameters (Parameter of information sharing 
     within the cluster) - We update only for the active clusters */
  
  // Active Clusters
  arma::uvec active_clus = arma::unique(ci);
  unsigned int Kpos = active_clus.size();
  
  // Update the beta vector for the active cluster.
  arma::mat beta_result(beta_old.n_rows, beta_old.n_cols, arma::fill::zeros);
  for(int kk = 0; kk < Kpos; ++kk){
    
    int k = active_clus[kk];
    arma::uvec ck = arma::find(ci == k);
    arma::rowvec bk_old = beta_old.row(k);
    
    // Proposed a new beta_k
    arma::rowvec bk_pro = arma::conv_to<arma::rowvec>::from(arma::mvnrnd(bk_old.t(), 
                                                                         s2_MH * arma::eye(z.n_cols, z.n_cols)));
    // Calculate logA
    double logA = 0.0;
    logA += arma::accu(arma::log_normpdf(bk_pro, mu, std::sqrt(s2)));
    logA += arma::accu(logmar_k(z, atrisk, bk_pro).rows(ck));
    logA -= arma::accu(arma::log_normpdf(bk_old, mu, std::sqrt(s2)));
    logA -= arma::accu(logmar_k(z, atrisk, bk_old).rows(ck));
    
    double logU = std::log(R::runif(0.01, 1.0));
    if(logU < logA){
      beta_result.row(k) = bk_pro;
    } else {
      beta_result.row(k) = bk_old;
    }
    
  }
  
  return beta_result;
  
} 

// Cluster Assignment
// [[Rcpp::export]]
Rcpp::List update_ci(unsigned int Kmax, arma::mat z, arma::mat atrisk, 
                     arma::mat beta_mat, arma::uvec ci_old, double theta, 
                     double mu, double s2){
  
  // Reallocate
  arma::uvec active_clus = arma::unique(ci_old); 
  //// Result from unique is already sorted from lowest to largest.
  unsigned int Kpos = active_clus.size();
  arma::mat lmar_realloc = logmar(z, atrisk, beta_mat);
  
  //// Find the number of the observations in each active cluster
  arma::vec nk(Kmax, arma::fill::zeros);
  for(int k = 0; k < Kmax; ++k){
    nk[k] += arma::accu(ci_old == k);
  }
  
  //// Start: Reallocation
  arma::uvec ci_realloc(ci_old);
  for(int i = 0; i < z.n_rows; ++i){
    int ci_old = ci_realloc[i];
    nk[ci_old] -= 1;
    arma::vec log_realloc_prob = arma::log(nk + theta) + lmar_realloc.row(i).t();
    
    // Change of index: From Kmax to the active only.
    arma::vec realloc_prob = log_sum_exp(log_realloc_prob.rows(active_clus));
    arma::vec ck_new_kk_vec = rmultinom_1(realloc_prob);
    arma::uvec kk_new_vec = arma::find(ck_new_kk_vec == 1);
    int kk_new = kk_new_vec[0];
    
    // Need to check the index back: From active to Kmax.
    int new_ci = active_clus[kk_new];
    nk[new_ci] += 1;
    ci_realloc[i] = new_ci;
  }
  
  
  // Split-Merge
  int split_index = -1;
  //// Start with determine to split (expand) or merge (collapse)
  //// samp_index: the index of two observations used for considering to split or merge.
  arma::uvec samp_ind = arma::randperm(z.n_rows, 2);
  while((Kpos == Kmax) and 
          (ci_realloc[samp_ind[0]] == ci_realloc[samp_ind[1]])){
    samp_ind = arma::randperm(z.n_rows, 2);
  }
  
  if(ci_realloc[samp_ind[0]] == ci_realloc[samp_ind[1]]){
    //// If split, proposed a new cluster with a new beta vector simulated from the prior
    split_index = 1;
    arma::uvec inactive_clus = arma::find(nk == 0);
    /// int new_active_clus = inactive_clus.row(arma::randperm(inactive_clus.size, 1));
    //// Dont't forget to adjust nk vec
    std::cout << inactive_clus << std::endl;
  } else {
    split_index = 0;
  }
  
  Rcpp::List result;
  result["ci_realloc"] = ci_realloc;
  result["split_index"] = split_index;
  result["samp_ind"] = samp_ind;
  return result;
  
}


// *****************************************************************************