#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#define pi 3.141592653589793238462643383280

// Note: -----------------------------------------------------------------------
// * Cluster index starts with 0.
// * Debugging -- check both log probability function and sampler function.
// * Add the description for the function.
// * Update the overleaf for the corresponding parameters.
// * Instead of using xi as a input, using beta as a input then take the exp().
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

