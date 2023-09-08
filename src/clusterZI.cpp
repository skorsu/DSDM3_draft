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
  
  /* Reallocation algorithm for the split merge */
  
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

// [[Rcpp::export]]
double log_proposal(arma::uvec clus_after, arma::uvec clus_before, 
                    arma::mat z, arma::mat gamma_mat, arma::mat beta_mat, 
                    arma::uvec S, arma::uvec clus_sm){
  
  /* Calculate the proposal probability, p(after|before), in a log scale */
  
  double log_val = 0.0;
  
  arma::vec nk(2, arma::fill::zeros);
  
  for(int ss = 0; ss < S.size(); ++ss){
    int s = S[ss];
    arma::uvec n_index = arma::find(clus_sm == clus_before[s]);
    nk[n_index[0]] += 1;
  }
  
  
  for(int ss = 0; ss < S.size(); ++ss){
    int s = S[ss];
    arma::uvec old_index = arma::find(clus_sm == clus_before[s]);
    nk[old_index[0]] -= 1;
    
    // Calculate the reallocation probability
    arma::vec zs = z.row(s).t();
    arma::vec gms = gamma_mat.row(s).t();
    arma::vec log_prob(2, arma::fill::zeros);
    
    for(int kk = 0; kk <= 1; ++kk){
      int k = clus_sm[kk];
      log_prob[kk] += log_marginal(zs, gms, beta_mat.row(k).t());
      log_prob[kk] += std::log(nk[kk]);
    }
    
    arma::vec prob = log_sum_exp(log_prob);
    arma::uvec new_index = arma::find(clus_sm == clus_after[s]);
    log_val += std::log(prob[new_index[0]]);
    nk[new_index[0]] += 1;
  }
  
  return log_val;
  
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
    launch_tau.row(new_ck) = R::rgamma(theta_vec[new_ck], 1.0);
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
  arma::uvec proposed_assign(launch_assign);
  if(expand_ind == 1){
    proposed_assign = realloc_sm(z, launch_assign, gamma_mat, launch_beta, S, samp_clus);
  } else {
    proposed_assign.rows(S).fill(samp_clus[1]);
    proposed_assign.rows(samp_ind).fill(samp_clus[1]);
  }
  
  Rcpp::List proposed_tb =  adjust_tau_beta(launch_beta, launch_tau, proposed_assign);
  arma::mat proposed_beta = proposed_tb["beta"];
  arma::vec proposed_tau = proposed_tb["tau"];
  
  // MH
  double logA = 0.0;
  arma::vec nk_old(K_max, arma::fill::zeros);
  arma::vec nk_proposed(K_max, arma::fill::zeros);
  
  for(int i = 0; i < z.n_rows; ++i){
    arma::vec zi = z.row(i).t();
    arma::vec gmi = gamma_mat.row(i).t();
    logA += log_marginal(zi, gmi, proposed_beta.row(proposed_assign[i]).t());
    logA -= log_marginal(zi, gmi, beta_mat.row(clus_assign[i]).t());
    
    nk_old[clus_assign[i]] += 1;
    nk_proposed[proposed_assign[i]] += 1;
    
  }
  
  arma::uvec proposed_active = arma::find(nk_proposed > 0); 
  logA += std::lgamma(arma::accu(theta_vec.rows(proposed_active)));
  logA -= arma::accu(arma::lgamma(theta_vec.rows(proposed_active)));
  logA += arma::accu(arma::lgamma(nk_proposed.rows(proposed_active) + theta_vec.rows(proposed_active)));
  logA -= std::lgamma(arma::accu(nk_proposed.rows(proposed_active) + theta_vec.rows(proposed_active)));
  
  logA -= std::lgamma(arma::accu(theta_vec.rows(active_clus)));
  logA += arma::accu(arma::lgamma(theta_vec.rows(active_clus)));
  logA -= arma::accu(arma::lgamma(nk_old.rows(active_clus) + theta_vec.rows(active_clus)));
  logA += std::lgamma(arma::accu(nk_old.rows(active_clus) + theta_vec.rows(active_clus)));
  
  logA += arma::accu(arma::log_normpdf(proposed_beta, mu, std::sqrt(s2)));
  logA -= arma::accu(arma::log_normpdf(beta_mat, mu, std::sqrt(s2)));
  
  logA += log_proposal(launch_assign, proposed_assign, z, gamma_mat, 
                       launch_beta, S, samp_clus);
  if(expand_ind == 1){
    logA -= log_proposal(proposed_assign, launch_assign, z, gamma_mat, 
                         launch_beta, S, samp_clus);
  }
  
  // MH
  double logU = std::log(R::runif(0.0, 1.0));
  int sm_accept = 0;
  arma::uvec new_assign(clus_assign);
  arma::vec new_tau(tau_vec);
  arma::mat new_beta(beta_mat);
  if(logU <= logA){
    sm_accept += 1;
    new_assign = proposed_assign;
    new_beta = proposed_beta;
    new_tau = proposed_tau;
  }
  
  Rcpp::List result;
  result["S"] = S;
  result["logA"] = logA;
  result["expand_ind"] = expand_ind;
  result["sm_accept"] = sm_accept;
  result["assign"] = new_assign;
  result["tau"] = new_tau;
  result["beta"] = new_beta;
  return result;
  
}

// [[Rcpp::export]]
Rcpp::List update_tau(arma::uvec clus_assign, arma::vec tau_vec, 
                      arma::vec theta_vec, double U){
  
  /* Update tau and U */
  
  arma::vec new_tau(tau_vec);
  arma::uvec active_clus = arma::unique(clus_assign);
  double scale_U = 1/(1 + U);
  
  for(int kk = 0; kk < active_clus.size(); ++kk){
    int k = active_clus[kk];
    arma::uvec nk = arma::find(clus_assign == k);
    new_tau.row(k).fill(R::rgamma(nk.size() + theta_vec[k], scale_U)); 
  }
  
  double scale_u = 1/arma::accu(new_tau);
  double new_U = R::rgamma(clus_assign.size(), scale_u);
  
  Rcpp::List result;
  result["tau"] = new_tau;
  result["U"] = new_U;
  return result;
}

// *****************************************************************************
// [[Rcpp::export]]
arma::mat DM_DM(unsigned int iter, unsigned int K_max, arma::mat z,
                arma::vec theta_vec, double MH_var, double mu, double s2,
                int print_iter){
  
  /* This is one of our competitive model. We have to specify the number of 
  clusters, and we did not update the at-risk indicator. */
  
  arma::mat clus_iter(iter, z.n_rows, arma::fill::value(-1));
  
  // Initial the cluster assignment and beta matrix
  arma::uvec ci_init(z.n_rows, arma::fill::zeros);
  arma::mat beta_init(K_max, z.n_cols, arma::fill::ones);
  
  arma::mat gamma_mat(z.n_rows, z.n_cols, arma::fill::ones);
  
  // MCMC object
  arma::mat beta_mcmc(beta_init);
  
  for(int t = 0; t < iter; ++t){
    
    // Update beta
    beta_mcmc = update_beta(z, ci_init, gamma_mat, beta_init, mu, s2, MH_var);
    
    // Reallocate
    arma::uvec ci_mcmc(ci_init);
    arma::vec nk(K_max, arma::fill::zeros);
    
    // Count the number of element in each cluster
    for(int i = 0; i < z.n_rows; ++i){
      nk[ci_init[i]] += 1;
    }
    
    for(int i = 0; i < z.n_rows; ++i){
      
      nk[ci_init[i]] -= 1;
      
      arma::vec zi = z.row(i).t();
      arma::vec gmi = gamma_mat.row(i).t();
      arma::vec log_prob(K_max, arma::fill::zeros);
      
      for(int k = 0; k < K_max; ++k){
        log_prob[k] += log_marginal(zi, gmi, beta_mcmc.row(k).t());
        log_prob[k] += std::log(theta_vec[k] + nk[k]);
      }
      
      arma::vec realloc_prob = log_sum_exp(log_prob);
      Rcpp::NumericVector realloc_rcpp = Rcpp::NumericVector(realloc_prob.begin(), 
                                                             realloc_prob.end());
      Rcpp::IntegerVector realloc_rcpp_index = rmultinom_1(realloc_rcpp, K_max);
      arma::vec realloc_index = Rcpp::as<arma::vec>(Rcpp::wrap(realloc_rcpp_index));
      
      // New assign
      arma::uvec new_ck = arma::find(realloc_index == 1);
      ci_mcmc.row(i).fill(new_ck[0]);
      
      nk[ci_mcmc[i]] += 1;
      
    }
    
    for(int i = 0; i < z.n_rows; ++i){
      clus_iter.row(t).col(i).fill(ci_mcmc[i]);
    }
    
    ci_init = ci_mcmc;
    beta_init = beta_mcmc;
    
    // Print the result
    if(((t + 1) - (floor((t + 1)/print_iter) * print_iter)) == 0){
      std::cout << "Iter: " << (t+1) << " - Done!" << std::endl;
    }
    
  }
  
  return clus_iter;
  
}


// [[Rcpp::export]]
Rcpp::List DM_ZIDM(unsigned int iter, unsigned int K_max, arma::mat z,
                   arma::vec theta_vec, unsigned int launch_iter,
                   double MH_var, double mu, double s2, 
                   double r0c, double r1c, int print_iter){
  
  /* This is one of our competitive model. We include the SM for the cluster
     space, but we did not update the at-risk indicator. */
  
  // Store the result
  arma::mat clus_iter(iter, z.n_rows, arma::fill::value(-1));
  arma::cube beta_iter(K_max, z.n_cols, iter);
  arma::vec sm_iter(iter, arma::fill::zeros); 
  arma::vec accept_iter(iter, arma::fill::zeros);
  arma::vec logA_sm_iter(iter, arma::fill::zeros);
  
  arma::mat gamma_mat(z.n_rows, z.n_cols, arma::fill::ones);
  
  // Initialize
  arma::uvec ci_init(z.n_rows, arma::fill::zeros);
  arma::mat beta_init(K_max, z.n_cols, arma::fill::ones);
  arma::vec tau_init(K_max, arma::fill::zeros);
  tau_init.row(0).fill(R::rgamma(theta_vec[0], 1.0));
  double U_init = R::rgamma(z.n_rows, 1/(arma::accu(tau_init)));
  
  // MCMC object
  arma::mat beta_mcmc(beta_init);
  Rcpp::List realloc_List;
  Rcpp::List sm_List;
  Rcpp::List tau_List;
  
  // Begin
  for(int t = 0; t < iter; ++t){
    
    // Update beta
    beta_mcmc = update_beta(z, ci_init, gamma_mat, beta_init, mu, s2, MH_var);
    
    // Reallocate
    realloc_List = realloc(z, ci_init, gamma_mat, beta_mcmc, tau_init, theta_vec);
    arma::uvec ci_realloc = realloc_List["assign"];
    arma::vec tau_realloc = realloc_List["tau"];
    arma::mat beta_realloc = realloc_List["beta"];
    
    // Split-Merge
    sm_List = sm(K_max, z, ci_realloc, gamma_mat, beta_mcmc, tau_realloc, 
                 theta_vec, launch_iter, mu, s2, r0c, r1c);
    
    logA_sm_iter.row(t).fill(sm_List["logA"]);
    sm_iter.row(t).fill(sm_List["expand_ind"]);
    accept_iter.row(t).fill(sm_List["sm_accept"]);
    arma::uvec ci_sm = sm_List["assign"];
    arma::vec tau_sm = sm_List["tau"];
    arma::mat beta_sm = sm_List["beta"];
    
    // Update tau and U
    tau_List = update_tau(ci_sm, tau_sm, theta_vec, U_init);
    arma::vec tau_U = tau_List["tau"];
    double U_U = tau_List["U"];
    
    // Record the result
    beta_iter.slice(t) = beta_sm;
    
    for(int i = 0; i < z.n_rows; ++i){
      clus_iter.row(t).col(i).fill(ci_sm[i]);
    }
    
    // Update the initial value for the next iteration
    ci_init = ci_sm;
    beta_init = beta_sm;
    tau_init = tau_U;
    U_init = U_U;
    
    // Print the result
    if(((t + 1) - (floor((t + 1)/print_iter) * print_iter)) == 0){
      std::cout << "Iter: " << (t+1) << " - Done!" << std::endl;
    }
    
  }
  
  // Result
  Rcpp::List result;
  result["assign"] = clus_iter;
  result["sm"] = sm_iter;
  result["accept_iter"] = accept_iter;
  return result;
  
}

// [[Rcpp::export]]
Rcpp::List ZIDM_ZIDM(unsigned int iter, unsigned int K_max, arma::mat z,
                     arma::vec theta_vec, unsigned int launch_iter,
                     double MH_var, double mu, double s2, double r0g, double r1g, 
                     double r0c, double r1c, int print_iter){
  
  /* This is our model. Update at-risk indicator and include the SM for 
     the cluster space. */
  
  // Store the result
  arma::mat clus_iter(iter, z.n_rows, arma::fill::value(-1));
  arma::cube gamma_iter(z.n_rows, z.n_cols, iter);
  arma::cube beta_iter(K_max, z.n_cols, iter);
  arma::vec sm_iter(iter, arma::fill::zeros); 
  arma::vec accept_iter(iter, arma::fill::zeros);
  arma::vec logA_sm_iter(iter, arma::fill::zeros);
  
  // Initialize
  arma::uvec ci_init(z.n_rows, arma::fill::zeros);
  arma::mat gamma_init(z.n_rows, z.n_cols, arma::fill::ones);
  arma::mat beta_init(K_max, z.n_cols, arma::fill::ones);
  arma::vec tau_init(K_max, arma::fill::zeros);
  tau_init.row(0).fill(R::rgamma(theta_vec[0], 1.0));
  double U_init = R::rgamma(z.n_rows, 1/(arma::accu(tau_init)));
  
  // MCMC object
  arma::mat gamma_mcmc(gamma_init);
  arma::mat beta_mcmc(beta_init);
  Rcpp::List realloc_List;
  Rcpp::List sm_List;
  Rcpp::List tau_List;
  
  // Begin
  for(int t = 0; t < iter; ++t){
    
    // Update at-risk
    gamma_mcmc = update_at_risk(z, ci_init, gamma_init, beta_init, r0g, r1g);
    
    // Update beta
    beta_mcmc = update_beta(z, ci_init, gamma_mcmc, beta_init, mu, s2, MH_var);
    
    // Reallocate
    realloc_List = realloc(z, ci_init, gamma_mcmc, beta_mcmc, tau_init, theta_vec);
    arma::uvec ci_realloc = realloc_List["assign"];
    arma::vec tau_realloc = realloc_List["tau"];
    arma::mat beta_realloc = realloc_List["beta"];
    
    // Split-Merge
    sm_List = sm(K_max, z, ci_realloc, gamma_mcmc, beta_mcmc, tau_realloc, 
                 theta_vec, launch_iter, mu, s2, r0c, r1c);
    
    logA_sm_iter.row(t).fill(sm_List["logA"]);
    sm_iter.row(t).fill(sm_List["expand_ind"]);
    accept_iter.row(t).fill(sm_List["sm_accept"]);
    arma::uvec ci_sm = sm_List["assign"];
    arma::vec tau_sm = sm_List["tau"];
    arma::mat beta_sm = sm_List["beta"];
    
    // Update tau and U
    tau_List = update_tau(ci_sm, tau_sm, theta_vec, U_init);
    arma::vec tau_U = tau_List["tau"];
    double U_U = tau_List["U"];
    
    // Record the result
    gamma_iter.slice(t) = gamma_mcmc;
    beta_iter.slice(t) = beta_sm;
    
    for(int i = 0; i < z.n_rows; ++i){
      clus_iter.row(t).col(i).fill(ci_sm[i]);
    }
    
    // Update the initial value for the next iteration
    ci_init = ci_sm;
    gamma_init = gamma_mcmc;
    beta_init = beta_sm;
    tau_init = tau_U;
    U_init = U_U;
    
    // Print the result
    if(((t + 1) - (floor((t + 1)/print_iter) * print_iter)) == 0){
      std::cout << "Iter: " << (t+1) << " - Done!" << std::endl;
    }
    
  }
  
  // Result
  Rcpp::List result;
  // result["gamma"] = gamma_iter;
  // result["beta"] = beta_iter;
  result["assign"] = clus_iter;
  result["sm"] = sm_iter;
  // result["logA"] = logA_sm_iter;
  result["accept_iter"] = accept_iter;
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
  arma::mat b_init(K, z.n_cols, arma::fill::ones);
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
  arma::mat b_init(K, z.n_cols, arma::fill::ones);
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