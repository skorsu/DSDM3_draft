#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#define pi 3.141592653589793238462643383280

// Note: -----------------------------------------------------------------------
// * Cluster index starts with 0.
// (/) Debugging -- check both log probability function and sampler function.
// * Add the description for the function.
// * Update the overleaf for the corresponding parameters.
// (/) Instead of using xi as a input, using beta as a input then take the exp().
//
// * Check the function for the at-risk indicator
// ** log_prob_gamma_ijk: for the individual g_ijk (/)
// ** update_gamma: for all g (/)
// * Check the function for the important variable indicator
// ** log_wj (/)
// ** update_w (/)
// * Check the function for the beta
// ** log_beta_k: for every important variable in cluster k (/)
// ** update_beta: for all beta (/)
// * Check the function for the assignment update
// ** alloc (/)
// ** sm (/)
// ** update_tau_theta (/)
// -----------------------------------------------------------------------------

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
Rcpp::IntegerVector rmultinom_1(Rcpp::NumericVector &probs, unsigned int &N){
  
  /* Description: sample from the multinomial(N, probs).
   * Credit: https://gallery.rcpp.org/articles/recreating-rmultinom-and-rpois-with-rcpp/
   */
  
  Rcpp::IntegerVector outcome(N);
  rmultinom(1, probs.begin(), N, outcome.begin());
  return outcome;
}

// [[Rcpp::export]]
arma::vec adjust_tau(unsigned int K_max, arma::uvec clus_assign, 
                     arma::vec tau_vec){
  
  arma::vec a_tau = arma::zeros(K_max);
  
  arma::uvec new_active_clus = arma::unique(clus_assign);
  a_tau.elem(new_active_clus) = tau_vec.elem(new_active_clus);
  
  return a_tau;
}
  
// [[Rcpp::export]]
double log_g_ijk(int j, arma::vec zi, arma::vec gi, arma::vec w, arma::vec beta_k,
                 double b0g, double b1g){
  
  double result = 0.0;
  
  int g_ijk = gi[j];
  arma::vec gwx_i = gi % w % arma::exp(beta_k);
  arma::vec z_gwx_i = zi + gwx_i;
  
  result += std::log(R::beta(b0g + g_ijk, b1g + (1 - g_ijk)));
  result += std::lgamma(arma::accu(gwx_i));
  result -= arma::accu(arma::lgamma(gwx_i.replace(0.0, 1.0)));
  result -= std::lgamma(arma::accu(z_gwx_i));
  result += arma::accu(arma::lgamma(z_gwx_i.replace(0.0, 1.0)));
  
  return result;
  
}

// [[Rcpp::export]]
double log_w_j(int j, arma::mat z, arma::uvec clus_assign, arma::mat gamma_mat, 
               arma::vec w, arma::mat beta_mat, double b0w, double b1w){
  
  double result = 0.0;
  int wj = w[j];
  
  // Convert beta to xi by taking the exponential function.
  arma::mat xi = arma::exp(beta_mat);
  
  // Calculate the gamma_w_xi factor
  arma::mat w_mat = arma::repelem(w, 1, z.n_rows).t(); 
  arma::mat xi_obs = xi.rows(clus_assign);
  arma::mat gwx = gamma_mat % w_mat % xi_obs;
  arma::mat gwx_j = gwx.col(j);
  arma::vec z_gwx_j = gwx_j + z.col(j);
  
  result += R::lbeta(b0w + wj, b1w + (1 - wj));
  result += arma::accu(arma::lgamma(arma::sum(gwx, 1)));
  result -= arma::accu(arma::lgamma(arma::sum(z + gwx, 1)));
  
  result += arma::accu(arma::lgamma(z_gwx_j.replace(0.0, 1.0)));
  result -= arma::accu(arma::lgamma(gwx_j.replace(0.0, 1.0)));
  
  return result;
  
} 

// [[Rcpp::export]]
arma::vec log_beta_k(int k, arma::mat z, arma::uvec clus_assign, 
                     arma::mat gamma_mat, arma::vec w, arma::vec beta_k, 
                     double s2){
  
  // Filter only the observation which in the cluster k
  arma::uvec clus_index = arma::find(clus_assign == k);
  arma::mat z_active = z.rows(clus_index);
  arma::mat gamma_active = gamma_mat.rows(clus_index);
  
  // Filter only the important variables
  arma::uvec imp_var = arma::find(w == 1);
  z_active = z_active.cols(imp_var);
  gamma_active = gamma_active.cols(imp_var);
  arma::vec beta_active = beta_k.rows(imp_var);
  
  // Calculate the gamma_w_xi factor
  arma::mat xi_mat = arma::repelem(arma::exp(beta_active), 1, z_active.n_rows).t();
  arma::mat gwx = gamma_active % xi_mat;
  arma::mat z_gwx = z_active + gwx;
  
  Rcpp::NumericVector bb = Rcpp::NumericVector(beta_active.begin(), beta_active.end());
  Rcpp::NumericVector bb_d = Rcpp::dnorm4(bb, 0.0, std::sqrt(s2), true);
  
  arma::vec result = Rcpp::as<arma::vec>(Rcpp::wrap(bb_d));
  result += arma::vec(imp_var.size(), 
                      arma::fill::value(arma::accu(arma::lgamma(arma::sum(gwx, 1)))));
  result -= arma::vec(imp_var.size(), 
                      arma::fill::value(arma::accu(arma::lgamma(arma::sum(z_gwx, 1)))));
  
  result += arma::sum(arma::lgamma(z_gwx.replace(0.0, 1.0)), 0).t();
  result -= arma::sum(arma::lgamma(gwx.replace(0.0, 1.0)), 0).t();
  
  return result;
  
}

// [[Rcpp::export]]
Rcpp::List realloc(unsigned int K_max, arma::mat z, arma::uvec clus_assign, 
                   arma::mat gamma_mat, arma::vec w, arma::mat beta_mat,
                   arma::vec tau, arma::vec theta){
  
  arma::uvec active_clus = arma::unique(clus_assign);
  unsigned int K_pos = active_clus.size();
  arma::vec theta_clus = theta.rows(active_clus);
  
  // Filter only the important variable related components.
  arma::uvec imp_var = arma::find(w == 1);
  arma::mat z_active = z.cols(imp_var);
  arma::mat gamma_active = gamma_mat.cols(imp_var);
  arma::mat xi_active = arma::exp(beta_mat.cols(imp_var));
  
  arma::mat log_prob(clus_assign.size(), active_clus.size(), arma::fill::zeros);
  arma::vec n_clus(active_clus.size(), arma::fill::zeros);
  
  for(int kk = 0; kk < K_pos; ++kk){
    
    // Filter the component for the cluster k
    int k = active_clus[kk];
    arma::mat xi_mat = arma::repelem(xi_active.row(k).t(), 1, z_active.n_rows).t();
    arma::mat gwx_k = gamma_active % xi_mat;
    arma::mat z_gwx_k = z_active + gwx_k;
    
    // Calculate P(zi|ci = k) for all possible clusters in cluster k.
    arma::vec log_z(clus_assign.size(), arma::fill::zeros);
    log_z += arma::lgamma(arma::sum(gwx_k, 1));
    log_z -= arma::lgamma(arma::sum(z_gwx_k, 1));
    log_z += arma::sum(arma::lgamma(z_gwx_k.replace(0.0, 1.0)), 1);
    log_z -= arma::sum(arma::lgamma(gwx_k.replace(0.0, 1.0)), 1);
    
    log_prob.col(kk) = log_z;
    
    // Count the number of the observation in the cluster k
    arma::uvec nk = arma::find(clus_assign == k);
    n_clus.row(kk).fill(nk.size());
    
  }
  
  // Reallocation
  for(int i = 0; i < clus_assign.size(); ++i){
    
    // Exclude the ci from the count
    int ci = clus_assign[i];
    arma::uvec k = arma::find(active_clus == ci);
    n_clus[k[0]] -= 1;
    
    // Normalize log probability and reallocate the observation
    arma::vec norm_prob = log_sum_exp(arma::log(n_clus + theta_clus) + log_prob.row(i).t());
    log_prob.row(i) = norm_prob.t();
    Rcpp::NumericVector rcpp_prob = Rcpp::NumericVector(norm_prob.begin(), norm_prob.end());
    arma::vec assign = Rcpp::as<arma::vec>(Rcpp::wrap(rmultinom_1(rcpp_prob, K_pos)));
    
    arma::uvec new_k = arma::find(assign == 1);
    clus_assign[i] = active_clus[new_k[0]];
    arma::uvec index = arma::find(active_clus == clus_assign[i]);
    n_clus[index[0]] += 1;
    
  }
  
  // Update tau vector based on the new cluster assignment
  tau = adjust_tau(K_max, clus_assign, tau);
  
  Rcpp::List result;
  result["alloc_prob"] = log_prob;
  result["assign"] = clus_assign;
  result["tau"] = tau;
  return result;
  
}

// [[Rcpp::export]]
Rcpp::List realloc_sm(arma::mat z, arma::uvec clus_assign, 
                      arma::mat gamma_mat, arma::vec w, arma::mat beta_mat,
                      arma::uvec S, arma::uvec clus_sm){
  
  unsigned int K_sam = 2;
  
  // Filter only the important variable related components.
  arma::uvec imp_var = arma::find(w == 1);
  arma::mat z_active = z.cols(imp_var);
  arma::mat gamma_active = gamma_mat.cols(imp_var);
  arma::mat xi_active = arma::exp(beta_mat.cols(imp_var));
  
  arma::mat log_prob(S.size(), 2, arma::fill::zeros);
  arma::vec n_clus(2, arma::fill::zeros);
  
  for(int kk = 0; kk < 2; ++kk){
    
    // Filter the component for the cluster k
    int k = clus_sm[kk];
    arma::mat xi_mat = arma::repelem(xi_active.row(k).t(), 1, z_active.n_rows).t();
    arma::mat gwx_k = gamma_active % xi_mat;
    arma::mat z_gwx_k = z_active + gwx_k;
    
    // Calculate P(zi|ci = k) for all possible clusters in cluster k.
    arma::vec log_z(clus_assign.size(), arma::fill::zeros);
    log_z += arma::lgamma(arma::sum(gwx_k, 1));
    log_z -= arma::lgamma(arma::sum(z_gwx_k, 1));
    log_z += arma::sum(arma::lgamma(z_gwx_k.replace(0.0, 1.0)), 1);
    log_z -= arma::sum(arma::lgamma(gwx_k.replace(0.0, 1.0)), 1);
    
    log_prob.col(kk) = log_z.rows(S);
    
    // Count the number of the observation in the cluster k
    arma::uvec nk = arma::find(clus_assign == k);
    n_clus.row(kk).fill(nk.size());
    
  }
  
  // Reallocation
  for(int s = 0; s < S.size(); ++s){
    
    // Exclude the ci from the count
    int i = S[s];
    int ci = clus_assign[i];
    arma::uvec k = arma::find(clus_sm == ci);
    n_clus[k[0]] -= 1;
    
    // Normalize log probability and reallocate the observation
    arma::vec norm_prob = log_sum_exp(arma::log(n_clus) + log_prob.row(s).t());
    log_prob.row(s) = norm_prob.t();
    Rcpp::NumericVector rcpp_prob = Rcpp::NumericVector(norm_prob.begin(), norm_prob.end());
    arma::vec assign = Rcpp::as<arma::vec>(Rcpp::wrap(rmultinom_1(rcpp_prob, K_sam)));
    
    arma::uvec new_k = arma::find(assign == 1);
    clus_assign[i] = clus_sm[new_k[0]];
    arma::uvec index = arma::find(clus_sm == clus_assign[i]);
    n_clus[index[0]] += 1;
    
  }

  Rcpp::List result;
  result["alloc_prob"] = log_prob;
  result["clus_assign"] = clus_assign;
  return result;
  
}

// [[Rcpp::export]]
double log_proposal(arma::mat z, arma::uvec clus_target, arma::uvec clus_init, 
                    arma::mat gamma_mat, arma::vec w, arma::mat beta_mat,
                    arma::uvec S, arma::uvec clus_sm){
  
  double result = 0.0;
  
  // Filter only the important variable related components.
  arma::uvec imp_var = arma::find(w == 1);
  arma::mat z_active = z.cols(imp_var);
  arma::mat gamma_active = gamma_mat.cols(imp_var);
  arma::mat xi_active = arma::exp(beta_mat.cols(imp_var));
  
  arma::mat log_prob(S.size(), 2, arma::fill::zeros);
  arma::vec n_clus(2, arma::fill::zeros);
  
  for(int kk = 0; kk < 2; ++kk){
    
    // Filter the component for the cluster k
    int k = clus_sm[kk];
    arma::mat xi_mat = arma::repelem(xi_active.row(k).t(), 1, z_active.n_rows).t();
    arma::mat gwx_k = gamma_active % xi_mat;
    arma::mat z_gwx_k = z_active + gwx_k;
    
    // Calculate P(zi|ci = k) for all possible clusters in cluster k.
    arma::vec log_z(clus_init.size(), arma::fill::zeros);
    log_z += arma::lgamma(arma::sum(gwx_k, 1));
    log_z -= arma::lgamma(arma::sum(z_gwx_k, 1));
    log_z += arma::sum(arma::lgamma(z_gwx_k.replace(0.0, 1.0)), 1);
    log_z -= arma::sum(arma::lgamma(gwx_k.replace(0.0, 1.0)), 1);
    
    log_prob.col(kk) = log_z.rows(S);
    
    // Count the number of the observation in the cluster k
    arma::uvec nk = arma::find(clus_init == k);
    n_clus.row(kk).fill(nk.size());
    
  }
  
  // Calculate the log proposal
  for(int s = 0; s < S.size(); ++s){
    
    // Exclude the ci from the count
    int i = S[s];
    int ci = clus_init[i];
    arma::uvec k = arma::find(clus_sm == ci);
    n_clus[k[0]] -= 1;
    
    // Normalize log probability and reallocate the observation
    arma::vec norm_prob = log_sum_exp(arma::log(n_clus) + log_prob.row(s).t());
    log_prob.row(s) = norm_prob.t();
    
    // Log Proposal
    int ci_target = clus_target[i];
    arma::uvec index = arma::find(clus_sm == clus_init[i]);
    result += std::log(norm_prob[index[0]]);

    // Reassign the observation based on the clus_target
    clus_init[i] = ci_target;
    n_clus[index[0]] += 1;
    
  }
  
  return result;
  
}

// [[Rcpp::export]]
Rcpp::List sm(unsigned int K_max, arma::mat z, arma::uvec clus_assign, 
              arma::mat gamma_mat, arma::vec w, arma::mat beta_mat,
              arma::vec tau, arma::vec theta, unsigned int launch_iter,
              double b0c, double b1c){
  
  arma::uvec active_clus = arma::unique(clus_assign);
  unsigned int K_pos = active_clus.size();
  unsigned int n = clus_assign.size();
  
  int split_ind = -1; // Indicate that we will perform split or merge
  
  arma::uvec launch_assign(clus_assign);
  arma::vec launch_tau(tau);
  
  // Select two observations to determine whether we will split or merge.
  arma::uvec samp_obs = arma::randperm(n, 2);
  
  // If all cluster is already active, we can perform only merge.
  while((K_pos == K_max) and (launch_assign[samp_obs[0]] == launch_assign[samp_obs[1]])){
    samp_obs = arma::randperm(n, 2); // sample two observations again until we get merge.
  }
  arma::uvec samp_clus = launch_assign.rows(samp_obs);
  
  int new_clus = -1;
  // Create the split indicator. If we split, split_ind = 1. Otherwise, split_ind = 0.
  if(samp_clus[0] == samp_clus[1]){
    split_ind = 1; // Split
    
    arma::uvec inactive_index = arma::find(tau == 0);
    arma::uvec new_clus_index = arma::randperm(inactive_index.size(), 1);
    new_clus = inactive_index[new_clus_index[0]]; // new active cluster
    samp_clus.row(0).fill(new_clus);
    
    launch_assign.row(samp_obs[0]).fill(new_clus); // set ci_launch to be a new cluster
    launch_tau.row(new_clus).fill(R::rgamma(theta[new_clus], 1.0));
  } else { 
    split_ind = 0;
  }
  
  // Create a set S := {same cluster as samp_obs, but not index_samp}
  arma::uvec S = arma::find(clus_assign == clus_assign[samp_obs[0]] or clus_assign == clus_assign[samp_obs[1]]);
  arma::uvec samp_obs_index = arma::find(S == samp_obs[0] or S == samp_obs[1]);
  S.shed_rows(samp_obs_index);
  
  // Randomly assign observation in S to be ci_launch or cj_launch.
  arma::vec init_ind = arma::randu(S.size(), arma::distr_param(0, 1));
  launch_assign.rows(S) = samp_clus.rows(arma::conv_to<arma::uvec>::from(init_ind > 0.5));
  
  // Launch Step
  for(int t = 0; t < launch_iter; ++t){
    Rcpp::List l_product = realloc_sm(z, launch_assign, gamma_mat, w, beta_mat, S, samp_clus);
    arma::uvec ll = l_product["clus_assign"];
    launch_assign = ll;
  }
  
  // Split-Merge
  arma::uvec proposed_assign(launch_assign);
  
  if(split_ind == 1){ // Perform another launch step
    Rcpp::List p_product = realloc_sm(z, launch_assign, gamma_mat, w, beta_mat, S, samp_clus);
    arma::uvec pp = p_product["clus_assign"];
    proposed_assign = pp;
  } else { // Perform Merge
    proposed_assign.rows(S).fill(samp_clus[1]);
    proposed_assign.rows(samp_obs).fill(samp_clus[1]);
  }
  
  // LogA Calculation
  double logA = 0.0; // log of acceptance probability
  
  arma::uvec imp_var = arma::find(w == 1);
  arma::mat z_active = z.cols(imp_var);
  arma::mat gamma_active = gamma_mat.cols(imp_var);
  arma::mat xi_active = arma::exp(beta_mat.cols(imp_var));
  
  // Marginal Probability (Data): Proposed
  arma::mat gwx_proposed = gamma_active % xi_active.rows(proposed_assign);
  arma::mat z_gwx_proposed = z_active + gwx_proposed;
  arma::vec sum_gwx_proposed = arma::sum(gwx_proposed, 1);
  arma::vec sum_z_gwx_proposed = arma::sum(z_gwx_proposed, 1);
  
  logA += arma::accu(arma::lgamma(sum_gwx_proposed.replace(0, 1)));
  logA -= arma::accu(arma::lgamma(gwx_proposed.replace(0, 1)));
  logA += arma::accu(arma::lgamma(z_gwx_proposed.replace(0, 1)));
  logA -= arma::accu(arma::lgamma(sum_z_gwx_proposed.replace(0, 1)));
  
  // Marginal Probability (Data): Initial
  arma::mat gwx_init = gamma_active % xi_active.rows(clus_assign);
  arma::mat z_gwx_init = z_active + gwx_init;
  arma::vec sum_gwx_init = arma::sum(gwx_init, 1);
  arma::vec sum_z_gwx_init = arma::sum(z_gwx_init, 1);
  
  logA -= arma::accu(arma::lgamma(sum_gwx_init.replace(0, 1)));
  logA += arma::accu(arma::lgamma(gwx_init.replace(0, 1)));
  logA -= arma::accu(arma::lgamma(z_gwx_init.replace(0, 1)));
  logA += arma::accu(arma::lgamma(sum_z_gwx_init.replace(0, 1)));
  
  // Cluster Assignment part
  double sum_theta_init = 0.0;
  double sum_theta_proposed = 0.0;
  
  for(int k = 0; k < K_max; ++k){
    
    arma::uvec n_init = arma::find(clus_assign == k);
    arma::uvec n_proposed = arma::find(proposed_assign == k);
    
    if(n_proposed.size() > 0){
      logA += std::lgamma(n_proposed.size() + theta[k]);
      logA -= std::lgamma(theta[k]);
      sum_theta_proposed += theta[k];
    }
    
    if(n_init.size() > 0){
      logA -= std::lgamma(n_init.size() + theta[k]);
      logA += std::lgamma(theta[k]);
      sum_theta_init += theta[k];
    }
    
  }
  
  // Proposed
  logA += std::lgamma(sum_theta_proposed);
  logA -= std::lgamma(clus_assign.size() + sum_theta_proposed);
  
  // Initial
  logA -= std::lgamma(sum_theta_init);
  logA += std::lgamma(clus_assign.size() + sum_theta_init);
  
  // Proposal
  logA += log_proposal(z, launch_assign, proposed_assign, gamma_mat, w, 
                       beta_mat, S, samp_clus);
  
  if(split_ind == 1){ // Split
    logA -= log_proposal(z, proposed_assign, launch_assign, gamma_mat, w, 
                         beta_mat, S, samp_clus);
  }
  
  logA += (((2 * split_ind) - 1) * (std::log(b0c) - std::log(b1c)));
  
  // MH
  int accept_proposed = 0; // indicator for accepting the proposed assignment.
  arma::uvec new_assign(clus_assign);
  arma::vec new_tau(tau);
  double logU = std::log(R::runif(0.0, 1.0));
  if(logU <= logA){
    accept_proposed += 1;
    new_assign = proposed_assign;
    new_tau = adjust_tau(K_max, new_assign, launch_tau);
  }
  
  Rcpp::List result;
  result["split_ind"] = split_ind;
  result["logA"] = logA;
  result["accept_proposed"] = accept_proposed;
  result["assign"] = new_assign;
  result["tau"] = new_tau;
  return result;
  
}

// [[Rcpp::export]]
arma::mat update_gamma(arma::mat z, arma::uvec clus_assign, arma::mat gamma_mat, 
                       arma::vec w, arma::mat beta_mat, double b0g, double b1g){
  
  // Loop through observation
  for(int i = 0; i < clus_assign.size(); ++i){ 
    
    // Get the information for the observation i
    arma::vec zi_now = z.row(i).t();
    arma::vec beta_now = beta_mat.row(clus_assign[i]).t();
    arma::uvec zi_zero_w_one = arma::find(((zi_now == 0) % w) == 1);
    
    // Loop through z_ijk = 0 for only wj = 1
    for(int jj = 0; jj < zi_zero_w_one.size(); ++jj){ 
      
      // Propose a new gamma_ijk for jth location
      int j = zi_zero_w_one[jj];
      arma::vec old_g = gamma_mat.row(i).t();
      arma::vec proposed_g(old_g);
      proposed_g.row(j).fill(1 - old_g[j]);
      
      // MH
      double logA = 0.0;
      logA += log_g_ijk(j, zi_now, proposed_g, w, beta_now, b0g, b1g);
      logA -= log_g_ijk(j, zi_now, old_g, w, beta_now, b0g, b1g);
      double logU = std::log(R::runif(0.0, 1.0));
      
      if(logU <= logA){
        gamma_mat.row(i) = proposed_g.t();
      }
      
    }
    
  }
  
  return gamma_mat;
  
}

// [[Rcpp::export]]
arma::vec update_w(arma::mat z, arma::uvec clus_assign, arma::mat gamma_mat,
                   arma::vec w, arma::mat beta_mat, double b0w, double b1w){
  
  arma::vec proposed_w(w);
  
  // Get the index of the important and not important variables
  arma::uvec imp_index = arma::find(proposed_w == 1);
  arma::uvec unimp_index = arma::find(proposed_w == 0);
  int adjusted_index = -1;
  int j_update = -1;
  
  // Propose a new vector of the important variable indicator
  if(imp_index.size() <= 1){ // If only one variable is active, it cannot delete that one out.
    adjusted_index = arma::randperm(unimp_index.size(), 1)[0];
    j_update = unimp_index[adjusted_index];
    proposed_w.row(j_update).fill(1);
  } else {
    adjusted_index = arma::randperm(w.size(), 1)[0];
    j_update = adjusted_index;
    proposed_w.row(j_update).fill(1 - w[j_update]);
  }
  
  // MH
  double logA = 0.0;
  logA += log_w_j(j_update, z, clus_assign, gamma_mat, proposed_w, beta_mat, 
                  b0w, b1w);
  logA -= log_w_j(j_update, z, clus_assign, gamma_mat, w, beta_mat, b0w, b1w);
  double logU = std::log(R::runif(0.0, 1.0));
  
  if(logU <= logA){
    w = proposed_w;
  }

  return w;
  
}

// [[Rcpp::export]]
arma::mat update_beta(arma::mat z, arma::uvec clus_assign, arma::mat gamma_mat, 
                      arma::vec w, arma::mat beta_mat, double MH_var, double s2){
  
  arma::uvec active_clus = arma::unique(clus_assign);
  arma::uvec imp_var = arma::find(w == 1);
  
  // Loop through the active cluster
  for(int kk = 0; kk < active_clus.size(); ++kk){
    int k = active_clus[kk];
    arma::vec proposed_beta_k = beta_mat.row(k).t();
    
    // Propose a new value of beta
    for(int jj = 0; jj < imp_var.size(); ++jj){
      int j = imp_var[jj];
      double new_beta = R::rnorm(proposed_beta_k[j], std::sqrt(MH_var));
      proposed_beta_k.row(j).fill(new_beta);
    }
    
    // Calculate the log probability for all beta given that it is for the 
    // important variables.
    
    arma::vec logA(imp_var.size(), arma::fill::zeros);
    logA += log_beta_k(k, z, clus_assign, gamma_mat, w,  proposed_beta_k, s2);
    logA -= log_beta_k(k, z, clus_assign, gamma_mat, w,  beta_mat.row(k).t(), s2);
    
    arma::vec logU = arma::log(arma::randu(imp_var.size()));
    
    // MH
    for(int jj = 0; jj < imp_var.size(); ++jj){
      
      if(logU[jj] <= logA[jj]){
        int j = imp_var[jj];
        beta_mat(k, j) = proposed_beta_k[jj];
      }
      
    }
    
  }
  
  return beta_mat;

}

// [[Rcpp::export]]
Rcpp::List update_theta_u(arma::uvec clus_assign, arma::vec tau, 
                          arma::vec theta, double old_U){
  
  arma::vec new_tau(tau); 
  
  // update alpha vector
  arma::uvec active_clus = arma::unique(clus_assign);
  for(int k = 0; k < active_clus.size(); ++k){
    int current_c = active_clus[k];
    arma::uvec nk = arma::find(clus_assign == current_c);
    double scale_gamma = 1/(1 + old_U); // change the rate to scale parameter
    new_tau.row(current_c).fill(R::rgamma(nk.size() + theta[current_c], scale_gamma));
  }
  
  // update U
  int n = clus_assign.size();
  double scale_u = 1/arma::accu(new_tau);
  double new_U = R::rgamma(n, scale_u);
  
  Rcpp::List result;
  result["new_tau"] = new_tau;
  result["new_U"] = new_U;
  return result;
  
}

