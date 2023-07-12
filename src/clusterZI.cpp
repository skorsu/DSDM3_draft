#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#define pi 3.141592653589793238462643383280

// Note: -----------------------------------------------------------------------
// * Cluster index starts with 0.
// * Update at-risk indicator (done)
// * Update important variable indicator (done)
// * Update beta (done)
// * Update the cluster assignment
// ** reallocation (done)
// ** sm (wip)
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
arma::vec adjust_tau(int K_max, arma::uvec clus_assign, arma::vec tau_vec){
  arma::vec a_tau = arma::zeros(K_max);
  
  arma::uvec new_active_clus = arma::unique(clus_assign);
  a_tau.elem(new_active_clus) = tau_vec.elem(new_active_clus);
  
  return a_tau;
}
  
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
                     arma::vec w, arma::uvec clus_assign, double s2){
  
  // Filter only the active variables
  arma::uvec active_var = arma::find(w == 1);
  arma::mat z_active = z.cols(active_var);
  arma::mat gamma_active = gamma_mat.cols(active_var);
  arma::vec xi_active = arma::exp(beta_k.rows(active_var));
  arma::vec beta_k_active = beta_k.rows(active_var);

  // Filter only the observations in the same cluster (ci)
  arma::uvec clus_index = arma::find(clus_assign == ci);
  z_active = z_active.rows(clus_index);
  gamma_active = gamma_active.rows(clus_index);

  // Calculate the gamma_w_xi factor
  arma::mat xi_mat = arma::repelem(xi_active, 1, z_active.n_rows);
  xi_mat = xi_mat.t();
  
  arma::mat gwx = gamma_active % xi_mat;
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
Rcpp::List realloc(arma::mat z, arma::uvec clus_assign, arma::vec w, 
                   arma::mat gamma_mat, arma::mat beta, arma::vec tau, 
                   arma::vec theta){
  
  arma::uvec new_assign(clus_assign);
  arma::uvec active_clus = arma::unique(clus_assign);
  unsigned int K_pos = active_clus.size();
  
  // Number of the observation in each cluster
  arma::vec n_clus(K_pos, arma::fill::value(0.0));
  for(int i = 0; i < clus_assign.size(); ++i){
    arma::uvec index = arma::find(active_clus == clus_assign[i]);
    n_clus[index[0]] += 1;
  }

  // Filter only the active variables
  arma::uvec active_var = arma::find(w == 1);
  arma::mat z_active = z.cols(active_var);
  arma::mat gamma_active = gamma_mat.cols(active_var);
  arma::mat xi_active = arma::exp(beta.cols(active_var));
  
  // Select only the component for the active clusters
  xi_active = xi_active.rows(active_clus);
  arma::vec theta_active = theta.rows(active_clus);
  
  for(int i = 0; i < clus_assign.size(); ++i){
    
    // Adjust the n_clus vector
    arma::uvec index = arma::find(active_clus == new_assign[i]);
    n_clus[index[0]] -= 1;
    
    // Calculate gwx for the observation i in each active cluster
    arma::vec g_i = arma::conv_to<arma::vec>::from(gamma_active.row(i));
    arma::mat gwx_i = xi_active % arma::repelem(g_i, 1, K_pos).t();
    arma::mat z_gwx_i = gwx_i + arma::repelem(z_active.row(i), K_pos, 1);
    arma::vec sum_gwx_i = arma::sum(gwx_i, 1);
    arma::vec sum_z_gwx_i = arma::sum(z_gwx_i, 1);
    
    // Proportional of the log probability
    arma::vec log_prob = arma::lgamma(sum_gwx_i.replace(0.0, 1.0)) -
      arma::sum(arma::lgamma(gwx_i.replace(0.0, 1.0)), 1) +
      arma::sum(arma::lgamma(z_gwx_i.replace(0.0, 1.0)), 1) -
      arma::lgamma(sum_z_gwx_i.replace(0.0, 1.0)) +
      arma::log(n_clus + theta_active);
    arma::vec prob = log_sum_exp(log_prob); 
    
    // Reassign
    Rcpp::NumericVector rcpp_prob = Rcpp::NumericVector(prob.begin(), prob.end());
    arma::vec assign = Rcpp::as<arma::vec>(Rcpp::wrap(rmultinom_1(rcpp_prob, K_pos)));
    
    arma::uvec new_index = arma::find(assign == 1);
    new_assign[i] = active_clus[new_index[0]];
    index = arma::find(active_clus == new_assign[i]);
    n_clus[index[0]] += 1;
    
  }
  
  // Adjust tau
  unsigned int K_max = beta.n_rows;
  arma::vec new_tau = adjust_tau(K_max, new_assign, tau);
  
  Rcpp::List result;
  result["new_assign"] = new_assign;
  result["new_tau"] = new_tau;
  return result;
  
}

// [[Rcpp::export]]
Rcpp::List sm(arma::mat z, arma::uvec clus_assign, arma::vec w, 
              arma::mat gamma_mat, arma::mat beta, arma::vec tau, 
              arma::vec theta, int launch_iter){
  
  arma::uvec new_assign(clus_assign);
  arma::uvec active_clus = arma::unique(clus_assign);
  unsigned int K_pos = active_clus.size();
  unsigned int K_max = beta.n_rows;
  unsigned int n = clus_assign.size();
  
  int split_ind = -1; // Indicate that we will perform split or merge
  int accept_proposed = 0; // indicator for accepting the proposed assignment.
  double logA = 0.0; // log of acceptance probability
  
  arma::uvec launch_assign(clus_assign);
  arma::vec launch_tau(tau);
  
  // Select two observations to determine whether we will split or merge.
  arma::uvec samp_obs = arma::randperm(n, 2);
  
  // If all cluster is already active, we can perform only merge.
  while((K_pos == K_max) and (clus_assign[samp_obs[0]] == clus_assign[samp_obs[1]])){
    samp_obs = arma::randperm(n, 2); // sample two observations again until we get merge.
  }
  arma::uvec samp_clus = clus_assign.rows(samp_obs);
  
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
  
  // Perform a launch step
  // Number of the observation in each cluster
  arma::vec n_clus(2, arma::fill::value(1.0));
  for(int s = 0; s < S.size(); ++s){
    int cs = launch_assign[S[s]];
    arma::uvec index = arma::find(samp_clus == cs);
    n_clus[index[0]] += 1;
  }
  
  // Filter only the active variables
  arma::uvec active_var = arma::find(w == 1);
  arma::mat z_active = z.cols(active_var);
  arma::mat gamma_active = gamma_mat.cols(active_var);
  arma::mat xi_active = arma::exp(beta.cols(active_var));
  
  // Select only the component for the launch clusters
  xi_active = xi_active.rows(samp_clus);
  
  arma::vec g_i;
  arma::mat gwx_i;
  arma::mat z_gwx_i;
  arma::vec sum_gwx_i;
  arma::vec sum_z_gwx_i;
  
  for(int t = 0; t < launch_iter; ++t){
    for(int s = 0; s < S.size(); ++s){
      
      // Adjust the n_clus vector
      arma::uvec index = arma::find(samp_clus == launch_assign[S[s]]);
      n_clus[index[0]] -= 1;
      
      // Calculate gwx for the observation i in each active cluster
      g_i = arma::conv_to<arma::vec>::from(gamma_active.row(S[s]));
      gwx_i = xi_active % arma::repelem(g_i, 1, 2).t();
      z_gwx_i = gwx_i + arma::repelem(z_active.row(S[s]), 2, 1);
      sum_gwx_i = arma::sum(gwx_i, 1);
      sum_z_gwx_i = arma::sum(z_gwx_i, 1);
      
      // Proportional of the log probability
      arma::vec log_prob = arma::lgamma(sum_gwx_i.replace(0.0, 1.0)) -
        arma::sum(arma::lgamma(gwx_i.replace(0.0, 1.0)), 1) +
        arma::sum(arma::lgamma(z_gwx_i.replace(0.0, 1.0)), 1) -
        arma::lgamma(sum_z_gwx_i.replace(0.0, 1.0)) +
        arma::log(n_clus);
      arma::vec prob = log_sum_exp(log_prob);

      // Reassign
      Rcpp::NumericVector rcpp_prob = Rcpp::NumericVector(prob.begin(), prob.end());
      arma::vec assign = Rcpp::as<arma::vec>(Rcpp::wrap(rmultinom_1(rcpp_prob, K_pos)));

      arma::uvec new_index = arma::find(assign == 1);
      launch_assign[S[s]] = samp_clus[new_index[0]];
      index = arma::find(samp_clus == launch_assign[S[s]]);
      n_clus[index[0]] += 1;
      }
    }
  
  // Split-Merge
  arma::uvec sm_assign(launch_assign);
  
  if(split_ind == 1){
    
    for(int s = 0; s < S.size(); ++s){
      
      // Adjust the n_clus vector
      arma::uvec index = arma::find(samp_clus == sm_assign[S[s]]);
      n_clus[index[0]] -= 1;
      
      // Calculate gwx for the observation i in each active cluster
      g_i = arma::conv_to<arma::vec>::from(gamma_active.row(S[s]));
      gwx_i = xi_active % arma::repelem(g_i, 1, 2).t();
      z_gwx_i = gwx_i + arma::repelem(z_active.row(S[s]), 2, 1);
      sum_gwx_i = arma::sum(gwx_i, 1);
      sum_z_gwx_i = arma::sum(z_gwx_i, 1);
      
      // Proportional of the log probability
      arma::vec log_prob = arma::lgamma(sum_gwx_i.replace(0.0, 1.0)) -
        arma::sum(arma::lgamma(gwx_i.replace(0.0, 1.0)), 1) +
        arma::sum(arma::lgamma(z_gwx_i.replace(0.0, 1.0)), 1) -
        arma::lgamma(sum_z_gwx_i.replace(0.0, 1.0)) +
        arma::log(n_clus);
      arma::vec prob = log_sum_exp(log_prob);
      
      // Reassign
      Rcpp::NumericVector rcpp_prob = Rcpp::NumericVector(prob.begin(), prob.end());
      arma::vec assign = Rcpp::as<arma::vec>(Rcpp::wrap(rmultinom_1(rcpp_prob, K_pos)));
      
      arma::uvec new_index = arma::find(assign == 1);
      sm_assign[S[s]] = samp_clus[new_index[0]];
      index = arma::find(samp_clus == sm_assign[S[s]]);
      n_clus[index[0]] += 1;
    }
  } else {
    sm_assign.rows(S).fill(samp_clus[1]);
    sm_assign.rows(samp_obs).fill(samp_clus[1]);
    n_clus[0] -= (S.size() + 2);
    n_clus[1] += (S.size() + 2);
  }
  
  // Calculate the acceptance probability
  arma::vec n_old(K_max, arma::fill::value(0.0));
  arma::vec n_proposed(K_max, arma::fill::value(0.0)); 
  for(int i = 0; i < clus_assign.size(); ++i){
    n_old[clus_assign[i]] += 1;
    n_proposed[sm_assign[i]] += 1;
  }
  arma::uvec unique_old = arma::unique(clus_assign);
  arma::uvec unique_proposed = arma::unique(sm_assign);
  
  // Proposed: ci part
  logA += std::lgamma(arma::accu(theta.rows(unique_proposed)));
  logA -= arma::accu(arma::lgamma(theta.rows(unique_proposed)));
  
  arma::vec tn_proposed = n_proposed.rows(unique_proposed) + theta.rows(unique_proposed);
  
  logA += arma::accu(arma::lgamma(tn_proposed));
  logA -= std::lgamma(arma::accu(tn_proposed));
  
  // Old: ci part
  logA -= std::lgamma(arma::accu(theta.rows(unique_old)));
  logA += arma::accu(arma::lgamma(theta.rows(unique_old)));
  
  arma::vec tn_old = n_old.rows(unique_old) + theta.rows(unique_old);
  
  logA -= arma::accu(arma::lgamma(tn_old));
  logA += std::lgamma(arma::accu(tn_old));
  
  // Proposed: zi part
  
  arma::mat gwx_proposed = gamma_active % xi_active.rows(sm_assign);
  arma::mat z_gwx_proposed = z_active + gwx_proposed;
  arma::vec sum_gwx_proposed = arma::sum(gwx_proposed, 1);
  arma::vec sum_z_gwx_proposed = arma::sum(z_gwx_proposed, 1);
  
  logA += arma::accu(arma::lgamma(sum_gwx_proposed.replace(0, 1)));
  logA -= arma::accu(arma::lgamma(gwx_proposed.replace(0, 1)));
  logA += arma::accu(arma::lgamma(z_gwx_proposed.replace(0, 1)));
  logA -= arma::accu(arma::lgamma(sum_z_gwx_proposed.replace(0, 1)));
  
  // Old: zi part
  
  arma::mat gwx_old = gamma_active % xi_active.rows(clus_assign);
  arma::mat z_gwx_old = z_active + gwx_old;
  arma::vec sum_gwx_old = arma::sum(gwx_old, 1);
  arma::vec sum_z_gwx_old = arma::sum(z_gwx_old, 1);
  
  logA -= arma::accu(arma::lgamma(sum_gwx_old.replace(0, 1)));
  logA += arma::accu(arma::lgamma(gwx_old.replace(0, 1)));
  logA -= arma::accu(arma::lgamma(z_gwx_old.replace(0, 1)));
  logA += arma::accu(arma::lgamma(sum_z_gwx_old.replace(0, 1)));

  Rcpp::List result;
  result["split_ind"] = split_ind;
  result["launch_assign"] = launch_assign;
  result["sm_assign"] = sm_assign;
  return result;
  
}

Rcpp::List update_theta_u(arma::uvec clus_assign, arma::vec tau_vec, 
                          double old_U, arma::vec theta){
  
  arma::vec new_tau(tau_vec); 
  
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
arma::mat update_beta(arma::mat z, arma::uvec clus_assign, arma::vec w,
                      arma::mat gamma_mat, arma::mat old_beta, double mh_var,
                      double s2){
  
  arma::mat proposed_beta(old_beta);
  arma::mat new_beta(old_beta);
  
  // Filter only the active variables and active clusters
  arma::uvec active_var = arma::find(w == 1);
  arma::uvec active_clus = arma::conv_to<arma::uvec>::from(arma::unique(clus_assign));
  
  // (1) Sample a new beta
  for(int k = 0; k < active_clus.size(); ++k){
    
    int cc = active_clus[k];
    
    for(int j = 0; j < active_var.size(); ++j){
      
      int cv = active_var[j];
      double pb = R::rnorm(old_beta(cc, cv), std::sqrt(mh_var)); 
      proposed_beta.row(cc).col(cv).fill(pb);
      
    }
    
    // (2) Calculate the logA for all active betas in the cluster k
    arma::vec beta_k_p = arma::conv_to<arma::vec>::from(proposed_beta.row(cc));
    arma::vec beta_k_o = arma::conv_to<arma::vec>::from(old_beta.row(cc));
    arma::vec logA = log_beta_k(beta_k_p, cc, z, gamma_mat, w, clus_assign, s2) -
      log_beta_k(beta_k_o, cc, z, gamma_mat, w, clus_assign, s2);
    arma::vec logU = arma::log(arma::randu(logA.size()));
    
    // (3) Determine
    for(int j = 0; j < active_var.size(); ++j){
      
      int cv = active_var[j];
      
      if(logU[j] <= logA[j]){
        new_beta.row(cc).col(cv).fill(beta_k_p[j]);
      }
      
    }
    
  }
  
  return new_beta;
}
