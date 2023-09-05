#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#define pi 3.141592653589793238462643383280

// Note: -----------------------------------------------------------------------
// * Debugging (Dimension of the matrix/vector)
// -----------------------------------------------------------------------------

Rcpp::IntegerVector rmultinom_1(Rcpp::NumericVector &probs, unsigned int &N){
  
  /* Description: sample from the multinomial(N, probs).
   * Credit: https://gallery.rcpp.org/articles/recreating-rmultinom-and-rpois-with-rcpp/
   */
  
  Rcpp::IntegerVector outcome(N);
  rmultinom(1, probs.begin(), N, outcome.begin());
  return outcome;
}

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
Rcpp::List adjust_tau_beta(arma::mat beta_mat, arma::vec tau_vec,
                           arma::uvec clus_assign){
  
  /* Adjust tau and beta: let it be 0 for inactive cluster */
  
  arma::uvec active_clus = arma::unique(clus_assign);
  
  arma::vec new_tau_vec(tau_vec.size(), arma::fill::zeros);
  arma::mat new_beta_mat(beta_mat.n_rows, beta_mat.n_cols, arma::fill::zeros);
  
  new_tau_vec.elem(active_clus) = tau_vec.elem(active_clus);
  new_beta_mat.rows(active_clus) = beta_mat.rows(active_clus);
  
  Rcpp::List result;
  result["tau"] = new_tau_vec;
  result["beta"] = new_beta_mat;
  return result;
  
}

// [[Rcpp::export]]
double log_marginal(arma::vec zi, arma::vec gmi, arma::vec beta_k){
  
  /* Calculate the proportional of the marginal probability in a log scale. */
  
  double result = 0.0;
  
  arma::uvec gm1 = arma::find(gmi == 1);
  arma::vec xi_k = arma::exp(beta_k);
  
  result += std::lgamma(arma::accu(xi_k.rows(gm1)));
  result -= arma::accu(arma::lgamma(xi_k.rows(gm1)));
  result += arma::accu(arma::lgamma(zi.rows(gm1) + xi_k.rows(gm1)));
  result -= std::lgamma(arma::accu(zi.rows(gm1) + xi_k.rows(gm1)));
    
  return result;
  
}

// [[Rcpp::export]]
arma::uvec realloc_sm(arma::mat z, arma::uvec clus_assign, arma::mat gamma_mat, 
                      arma::mat beta_mat, arma::uvec S, arma::uvec clus_sm){
  
  arma::uvec new_assign(clus_assign); 
  arma::vec nk(2, arma::fill::zeros);
  unsigned int K_sm = clus_sm.size();
  
  for(int ss = 0; ss < S.size(); ++ss){
    int s = S[ss];
    arma::uvec n_index = arma::find(clus_sm == new_assign[s]);
    nk[n_index[0]] += 1;
  }
  
  for(int ss = 0; ss < S.size(); ++ss){

    int s = S[ss];
    arma::uvec n_index = arma::find(clus_sm == new_assign[s]);
    nk[n_index[0]] -= 1;

    arma::vec zs = z.row(s).t();
    arma::vec gms = gamma_mat.row(s).t();
    arma::vec log_prob(2, arma::fill::zeros);

    for(int kk = 0; kk <= 1; ++kk){
      int k = clus_sm[kk];
      log_prob[kk] += log_marginal(zs, gms, beta_mat.row(k).t());
      log_prob[kk] += std::log(nk[kk]);
    }

    arma::vec realloc_prob = log_sum_exp(log_prob);
    Rcpp::NumericVector realloc_rcpp = Rcpp::NumericVector(realloc_prob.begin(),
                                                           realloc_prob.end());
    Rcpp::IntegerVector realloc_rcpp_index = rmultinom_1(realloc_rcpp, K_sm);
    arma::vec realloc_index = Rcpp::as<arma::vec>(Rcpp::wrap(realloc_rcpp_index));

    // New assign
    arma::uvec new_ck = arma::find(realloc_index == 1);
    new_assign.row(s).fill(clus_sm[new_ck[0]]);

    nk[new_ck[0]] += 1;

  }
  
  return new_assign;
  
}

// *****************************************************************************
// [[Rcpp::export]]
arma::mat update_at_risk(arma::mat z, arma::uvec clus_assign, arma::mat gamma_mat, 
                         arma::mat beta_mat, double r0g, double r1g){
  
  /* Update the at-risk matrix. */
  
  arma::mat gamma_new(z.n_rows, z.n_cols, arma::fill::value(-1));
  
  for(int i = 0; i < z.n_rows; ++i){
    
    arma::vec zi = z.row(i).t();
    arma::vec gm_i = gamma_mat.row(i).t();
    arma::vec beta_k = beta_mat.row(clus_assign[i]).t(); 
    
    arma::uvec zi0 = arma::find(zi == 0);
    
    for(int k = 0; k < zi0.size(); ++k){
      int pp_gmk = 1 - gm_i[zi0[k]];
      arma::vec proposed_k(gm_i);
      proposed_k[zi0[k]] = pp_gmk;
      
      // Calculate logA
      double logA = 0.0;
      logA += R::lbeta(r0g + pp_gmk, r1g + (1 - pp_gmk));
      logA += log_marginal(zi, proposed_k, beta_k);
      logA -= R::lbeta(r0g + gm_i[zi0[k]], r1g + (1 - gm_i[zi0[k]]));
      logA -= log_marginal(zi, gm_i, beta_k);
      
      // MH
      double logU = std::log(R::runif(0.0, 1.0));
      if(logU <= logA){
        gm_i = proposed_k;
      }

    }
    
    gamma_new.row(i) = gm_i.t();
    
  }
  
  return gamma_new;
  
}

// [[Rcpp::export]]
arma::mat update_beta(arma::mat z, arma::uvec clus_assign, arma::mat gamma_mat, 
                      arma::mat beta_mat, double mu, double s2, double s2_MH){
  
  /* Update the beta matrix. */
  
  arma::mat beta_new(beta_mat);
  arma::uvec active_clus = arma::unique(clus_assign);
  
  arma::mat s2_MH_mat(beta_mat.n_cols, beta_mat.n_cols, arma::fill::eye);
  s2_MH_mat = s2_MH_mat * std::sqrt(s2_MH);
  
  for(int kk = 0; kk < active_clus.size(); ++kk){
    int k = active_clus[kk];
    double logA = 0.0;
    
    // Propose a new beta_k
    arma::vec proposed_beta = arma::mvnrnd(beta_mat.row(k).t(), s2_MH_mat);
    
    // Calculate logA
    logA += arma::accu(arma::log_normpdf(proposed_beta, mu, std::sqrt(s2)));
    logA -= arma::accu(arma::log_normpdf(beta_mat.row(k).t(), mu, std::sqrt(s2)));
    
    arma::uvec index_k = arma::find(clus_assign == k);
    
    for(int ii = 0; ii < index_k.size(); ++ii){
      int i = index_k[ii];
      logA += log_marginal(z.row(i).t(), gamma_mat.row(i).t(), proposed_beta);
      logA -= log_marginal(z.row(i).t(), gamma_mat.row(i).t(), beta_mat.row(k).t());
    }
    
    // MH
    double logU = std::log(R::runif(0.0, 1.0));
    if(logU <= logA){
      beta_new.row(k) = proposed_beta.t();
    }
    
  }
  
  return beta_new;
  
}

// [[Rcpp::export]]
Rcpp::List realloc(arma::mat z, arma::uvec clus_assign,
                   arma::mat gamma_mat, arma::mat beta_mat,
                   arma::vec tau_vec, arma::vec theta_vec){
  
  /* Reallocate */
  
  arma::uvec new_assign(clus_assign);
  arma::uvec active_clus = arma::unique(clus_assign);
  unsigned int K_max = active_clus.size();
  arma::vec nk(active_clus.size(), arma::fill::zeros);
  
  // Count the number of element in each cluster
  for(int i = 0; i < z.n_rows; ++i){
    arma::uvec n_index = arma::find(active_clus == new_assign[i]);
    nk[n_index[0]] += 1;
  }
  
  // Reallocate
  for(int i = 0; i < z.n_rows; ++i){
    
    arma::uvec n_index = arma::find(active_clus == new_assign[i]);
    nk[n_index[0]] -= 1;
    
    arma::vec zi = z.row(i).t();
    arma::vec gmi = gamma_mat.row(i).t();
    arma::vec log_prob(K_max, arma::fill::zeros);
    
    for(int kk = 0; kk < K_max; ++kk){
      int k = active_clus[kk];
      log_prob[kk] += log_marginal(zi, gmi, beta_mat.row(k).t());
      log_prob[kk] += std::log(theta_vec[k] + nk[kk]);
    }
    
    arma::vec realloc_prob = log_sum_exp(log_prob);
    Rcpp::NumericVector realloc_rcpp = Rcpp::NumericVector(realloc_prob.begin(), 
                                                           realloc_prob.end());
    Rcpp::IntegerVector realloc_rcpp_index = rmultinom_1(realloc_rcpp, K_max);
    arma::vec realloc_index = Rcpp::as<arma::vec>(Rcpp::wrap(realloc_rcpp_index));
    
    // New assign
    arma::uvec new_ck = arma::find(realloc_index == 1);
    new_assign.row(i).fill(active_clus[new_ck[0]]);
    
    nk[new_ck[0]] += 1;
    
  }
  
  // Adjust tau and beta
  Rcpp::List new_tb = adjust_tau_beta(beta_mat, tau_vec, new_assign);
  arma::mat new_beta = new_tb["beta"];
  arma::vec new_tau = new_tb["tau"];
  
  Rcpp::List result;
  result["assign"] = new_assign;
  result["tau"] = new_tau;
  result["beta"] = new_beta;
  return result;
  
}

// [[Rcpp::export]]
Rcpp::List sm(unsigned int K_max, arma::mat z, arma::uvec clus_assign,
              arma::mat gamma_mat, arma::mat beta_mat, arma::vec tau_vec, 
              arma::vec theta_vec, unsigned int launch_iter,
              double mu, double s2, double r0c, double r1c){
  
  /* Expand/Collapse the cluster space via Split-Merge */
  unsigned int n = z.n_rows;
  arma::uvec active_clus = arma::unique(clus_assign);
  unsigned int K_pos = active_clus.size();
  int expand_ind = -1;
  
  // Decide to expand (split) or collapse (merge)
  arma::uvec samp_ind = arma::randperm(n, 2);
  while((K_pos == K_max) and 
          (clus_assign[samp_ind[0]] == clus_assign[samp_ind[1]])){
    samp_ind = arma::randperm(n, 2);
  }
  
  // Create a set S
  arma::uvec samp_clus = clus_assign.rows(samp_ind);
  arma::uvec S = arma::find(clus_assign == samp_clus[0] or
                              clus_assign == samp_clus[1]);
  arma::uvec notS = arma::find(S != samp_ind[0] and S != samp_ind[1]);
  S = S.rows(notS);
  
  arma::uvec launch_assign(clus_assign);
  arma::vec launch_tau(tau_vec);
  arma::mat launch_beta(beta_mat);
  
  if(samp_clus[0] == samp_clus[1]){ // Split
    expand_ind = 1;
    arma::uvec inactive_clus = arma::find(tau_vec == 0);
    int new_ck = inactive_clus[arma::randperm(inactive_clus.size(), 1)[0]];
    launch_assign.row(samp_ind[0]).fill(new_ck);
    samp_clus.row(0).fill(new_ck);
    launch_tau.row(new_ck).fill(R::rgamma(tau_vec[new_ck], 1.0));
    launch_beta.row(new_ck) = arma::randn(z.n_cols, arma::distr_param(mu, std::sqrt(s2))).t(); 
  } else { // Merge
    expand_ind = 0;
  }
  
  // Perform a launch step
  arma::vec rand_index = arma::randu(S.size());
  launch_assign.rows(S) = samp_clus.rows((rand_index >= 0.5));
  for(int t = 0; t <= launch_iter; ++t){
    launch_assign = realloc_sm(z, launch_assign, gamma_mat, launch_beta, S, samp_clus);
  }
  
  // Perform last SM
  arma::vec x;
  if(expand_ind == 1){
    
  }
  
  Rcpp::List result;
  result["samp_ind"] = samp_ind;
  result["samp_clus"] = samp_clus;
  result["launch_assign"] = launch_assign;
  result["launch_tau"] = launch_tau;
  result["launch_beta"] = launch_beta;
  result["expand_ind"] = expand_ind;
  result["S"] = S;
  return result;
  
}


// *****************************************************************************
// [[Rcpp::export]]
arma::cube beta_mat_update(unsigned int K, unsigned int iter, arma::mat z, 
                           arma::uvec clus_assign, double mu, double s2, 
                           double s2_MH){
  
  /* Try: only beta */
  
  arma::cube result(K, z.n_cols, iter);
  arma::mat gm(z.n_rows, z.n_cols, arma::fill::ones);
  
  // Initialize the beta matrix
  arma::mat b_init(K, z.n_cols, arma::fill::zeros);
  arma::mat b_mcmc(b_init);
  
  for(int t = 0; t < iter; ++t){
    b_mcmc = update_beta(z, clus_assign, gm, b_init, mu, s2, s2_MH);
    result.slice(t) = b_mcmc;
    b_init = b_mcmc;
  }
  
  return result;
  
}

// [[Rcpp::export]]
Rcpp::List beta_ar_update(unsigned int K, unsigned int iter, arma::mat z, 
                          arma::uvec clus_assign, double r0g, double r1g, 
                          double mu, double s2, double s2_MH){
  
  /* Try: both beta and at-risk */
  
  arma::cube at_risk_mat(z.n_rows, z.n_cols, iter);
  arma::cube beta_mat(K, z.n_cols, iter);
  
  // Initialize
  arma::mat gm_init(z.n_rows, z.n_cols, arma::fill::ones);
  arma::mat gm_mcmc(gm_init);
  arma::mat b_init(K, z.n_cols, arma::fill::zeros);
  arma::mat b_mcmc(b_init);
  
  for(int t = 0; t < iter; ++t){
    gm_mcmc = update_at_risk(z, clus_assign, gm_init, b_init, r0g, r1g);
    b_mcmc = update_beta(z, clus_assign, gm_mcmc, b_init, mu, s2, s2_MH);
    
    at_risk_mat.slice(t) = gm_mcmc;
    beta_mat.slice(t) = b_mcmc;
    
    gm_init = gm_mcmc;
    b_init = b_mcmc;
  }
  
  Rcpp::List result;
  result["gamma"] = at_risk_mat;
  result["beta"] = beta_mat;
  return result;
  
}

// *****************************************************************************



