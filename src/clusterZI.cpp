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
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
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
double log_g_ijk(int j, arma::vec zi, arma::vec gi, arma::vec w, arma::vec beta_k,
                 double r0g, double r1g){
  
  int g_ijk = gi[j];
  
  // Filter only the important variables
  arma::uvec imp_var = arma::find(w == 1);
  arma::vec zi_active = zi.rows(imp_var);
  arma::vec gi_active = gi.rows(imp_var);
  arma::vec xi_active = arma::exp(beta_k.rows(imp_var));
  
  // Calculate gwx
  arma::vec gwx_i = gi_active % xi_active;
  arma::vec z_gwx_i = zi_active + gwx_i;
  
  // Calculate the probability
  double result = 0.0;
  result += std::log(R::beta(r0g + g_ijk, r1g + (1 - g_ijk)));
  result += std::lgamma(arma::accu(gwx_i));
  result -= std::lgamma(arma::accu(z_gwx_i));
  
  return result;
  
}

// [[Rcpp::export]]
double log_w(arma::mat z, arma::uvec clus_assign, arma::mat gamma_mat, 
             arma::vec w, arma::mat beta_mat, double r0w, double r1w){
  
  double result = 0.0;
  
  result += (arma::accu(arma::lgamma(r0w + w) + arma::lgamma(r0w + (1 - w))));
  result -= (w.size() * std::lgamma(r0w + r1w + 1));
  
  // Filter only the important variables
  arma::uvec imp_var = arma::find(w == 1);
  arma::mat z_active = z.cols(imp_var);
  arma::mat gamma_active = gamma_mat.cols(imp_var);
  arma::mat xi_active = arma::exp(beta_mat.cols(imp_var));
  xi_active = xi_active.rows(clus_assign);
  
  // Calculate gwx
  arma::mat gwx = gamma_active % xi_active;
  arma::mat z_gwx = z_active + gwx;
  arma::vec sum_gwx = arma::sum(gwx, 1);
  arma::vec sum_z_gwx = arma::sum(z_gwx, 1);

  result += arma::accu(arma::lgamma(sum_gwx.replace(0.0, 1.0)));
  result -= arma::accu(arma::lgamma(sum_z_gwx.replace(0.0, 1.0)));
  result += arma::accu(arma::lgamma(z_gwx.replace(0.0, 1.0)));
  result -= arma::accu(arma::lgamma(gwx.replace(0.0, 1.0)));
  
  return result;
  
} 

// [[Rcpp::export]]
double log_beta_k(int k, arma::mat z, arma::uvec clus_assign, 
                  arma::mat gamma_mat, arma::vec w, arma::vec beta_k, 
                  double s2){
  
  // Filter only the important variables
  arma::uvec imp_var = arma::find(w == 1); 
  arma::mat z_active = z.cols(imp_var);
  arma::mat gamma_active = gamma_mat.cols(imp_var);
  arma::vec xi_active = arma::exp(beta_k.rows(imp_var));
  
  // Filter only the observation which in the cluster k
  arma::uvec clus_index = arma::find(clus_assign == k);
  z_active = z_active.rows(clus_index);
  gamma_active = gamma_active.rows(clus_index);
  
  // Calculate the probability of beta_jk
  Rcpp::NumericVector bb = Rcpp::NumericVector(beta_k.begin(), beta_k.end());
  Rcpp::NumericVector bb_d = Rcpp::dnorm4(bb, 0.0, std::sqrt(s2), true);
  
  // Calculate gwx
  arma::mat xi_mat = arma::repelem(xi_active, 1, z_active.n_rows).t();
  arma::mat gwx = gamma_active % xi_mat;
  arma::mat z_gwx = z_active + gwx;
  arma::vec sum_gwx = arma::sum(gwx, 1);
  arma::vec sum_z_gwx = arma::sum(z_gwx, 1);
  
  // Calculate the conditional probability
  double result = Rcpp::sum(bb_d);
  
  result += arma::accu(arma::lgamma(sum_gwx.replace(0.0, 1.0)));
  result -= arma::accu(arma::lgamma(sum_z_gwx.replace(0.0, 1.0)));
  result += arma::accu(arma::lgamma(z_gwx.replace(0.0, 1.0)));
  result -= arma::accu(arma::lgamma(gwx.replace(0.0, 1.0)));
  
  return result;
  
}

// [[Rcpp::export]]
arma::uvec realloc_sm(arma::mat z, arma::uvec clus_assign, 
                      arma::mat gamma_mat, arma::vec w, arma::mat beta_mat,
                      arma::uvec S, arma::uvec clus_sm){
  
  // Note: Similar to realloc function. (Loop through S and active_clus = clus_sm)
  
  unsigned int K = 2;
  
  // Filter only the important variables
  arma::uvec imp_var = arma::find(w == 1); 
  arma::mat z_active = z.cols(imp_var);
  arma::mat gamma_active = gamma_mat.cols(imp_var);
  arma::mat xi_active = arma::exp(beta_mat.cols(imp_var));
  
  // Count the number of the observation for each cluster
  arma::vec n_clus(K, arma::fill::zeros); 
  for(int ss = 0; ss < S.size(); ++ss){
    int s = S[ss];
    arma::uvec index_ci = arma::find(clus_sm == clus_assign[s]);
    n_clus[index_ci[0]] += 1;
  } 
  
  // Calculate the posterior predictive (= marginal distribution)
  arma::mat gwx(z_active.n_rows, imp_var.size(), arma::fill::zeros); 
  arma::mat z_gwx(z_active.n_rows, imp_var.size(), arma::fill::zeros);
  arma::vec sum_gwx(z_active.n_rows, arma::fill::zeros);
  arma::vec sum_z_gwx(z_active.n_rows, arma::fill::zeros);
  
  arma::mat log_postpred(z_active.n_rows, K, arma::fill::zeros);
  
  for(int kk = 0; kk < K; ++kk){
    
    int k = clus_sm[kk];
    
    gwx = gamma_active % (arma::repelem(xi_active.row(k).t(), 1, z_active.n_rows).t());
    z_gwx = z_active + gwx;
    sum_gwx = arma::sum(gwx, 1);
    sum_z_gwx = arma::sum(z_gwx, 1);
    
    arma::vec log_k(z_active.n_rows, arma::fill::zeros);
    log_k += arma::lgamma(sum_gwx.replace(0.0, 1.0));
    log_k -= arma::lgamma(sum_z_gwx.replace(0.0, 1.0));
    log_k += arma::sum(arma::lgamma(z_gwx.replace(0.0, 1.0)), 1);
    log_k -= arma::sum(arma::lgamma(gwx.replace(0.0, 1.0)), 1);
    
    log_postpred.col(kk) = log_k;
    
  }
  
  // Reallocate
  for(int ss = 0; ss < S.size(); ++ss){
    
    int s = S[ss];
    
    // Remove the i from the count
    int ci = clus_assign[s];
    arma::uvec index_ci_old = arma::find(clus_sm == ci);
    n_clus[index_ci_old[0]] -= 1;
    
    // Calculate the reallocation probability
    arma::vec log_realloc_prob = log_postpred.row(s).t();
    log_realloc_prob += arma::log(n_clus);
    arma::vec realloc_prob = log_sum_exp(log_realloc_prob);
    
    // Sampling the new cluster assignment
    Rcpp::NumericVector rprob = Rcpp::NumericVector(realloc_prob.begin(), realloc_prob.end());
    Rcpp::IntegerVector newci = rmultinom_1(rprob, K);
    arma::uvec new_ci = Rcpp::as<arma::uvec>(Rcpp::wrap(newci));
    arma::uvec index_ci_new = arma::find(new_ci == 1);
    
    // Reassign
    int ci_new = clus_sm[index_ci_new[0]];
    clus_assign[s] = ci_new;
    n_clus[index_ci_new[0]] += 1;
  } 

  return clus_assign;

}

// [[Rcpp::export]]
double log_proposal(arma::mat z, arma::uvec clus_target, arma::uvec clus_init, 
                    arma::mat gamma_mat, arma::vec w, arma::mat beta_mat,
                    arma::uvec S, arma::uvec clus_sm){
  
  double result = 0.0;
  
  unsigned int K = 2;
  
  // Filter only the important variables
  arma::uvec imp_var = arma::find(w == 1); 
  arma::mat z_active = z.cols(imp_var);
  arma::mat gamma_active = gamma_mat.cols(imp_var);
  arma::mat xi_active = arma::exp(beta_mat.cols(imp_var));
  
  // Count the number of the observation for each cluster from the initial
  // cluster assignment
  arma::vec n_clus(K, arma::fill::zeros); 
  for(int ss = 0; ss < S.size(); ++ss){
    int s = S[ss];
    arma::uvec index_ci = arma::find(clus_sm == clus_init[s]);
    n_clus[index_ci[0]] += 1;
  }  
  
  // Calculate the posterior predictive (= marginal distribution)
  arma::mat gwx(z_active.n_rows, imp_var.size(), arma::fill::zeros); 
  arma::mat z_gwx(z_active.n_rows, imp_var.size(), arma::fill::zeros);
  arma::vec sum_gwx(z_active.n_rows, arma::fill::zeros);
  arma::vec sum_z_gwx(z_active.n_rows, arma::fill::zeros);
  
  arma::mat log_postpred(z_active.n_rows, K, arma::fill::zeros);
  
  for(int kk = 0; kk < K; ++kk){
    
    int k = clus_sm[kk];
    
    gwx = gamma_active % (arma::repelem(xi_active.row(k).t(), 1, z_active.n_rows).t());
    z_gwx = z_active + gwx;
    sum_gwx = arma::sum(gwx, 1);
    sum_z_gwx = arma::sum(z_gwx, 1);
    
    arma::vec log_k(z_active.n_rows, arma::fill::zeros);
    log_k += arma::lgamma(sum_gwx.replace(0.0, 1.0));
    log_k -= arma::lgamma(sum_z_gwx.replace(0.0, 1.0));
    log_k += arma::sum(arma::lgamma(z_gwx.replace(0.0, 1.0)), 1);
    log_k -= arma::sum(arma::lgamma(gwx.replace(0.0, 1.0)), 1);
    
    log_postpred.col(kk) = log_k;
    
  }
  
  // Calculate the log proposal by looping only the observation in S
  for(int ss = 0; ss < S.size(); ++ss){
    
    int s = S[ss];
    int ci_old = clus_init[s];
    int ci_new = clus_target[s];
    
    arma::uvec index_ci_old = arma::find(clus_sm == ci_old);
    arma::uvec index_ci_new = arma::find(clus_sm == ci_new);
    
    // Remove the ci_old
    n_clus[index_ci_old[0]] -= 1;
    
    // Calculate the reallocation probability
    arma::vec log_realloc_prob = log_postpred.row(s).t();
    log_realloc_prob += arma::log(n_clus);
    arma::vec realloc_prob = log_sum_exp(log_realloc_prob);
    
    // Calculate the log proposal for the observation u
    result += std::log(realloc_prob[index_ci_new[0]]);
    
    // Add the ci_new back
    n_clus[index_ci_new[0]] += 1;
    
  }
  
  return result;
  
}

// [[Rcpp::export]]
arma::mat update_gamma(arma::mat z, arma::uvec clus_assign, arma::mat gamma_mat,
                       arma::vec w, arma::mat beta_mat, double r0g, double r1g){
  
  // Loop through the observation
  for(int i = 0; i < clus_assign.size(); ++i){
    
    arma::vec zi = z.row(i).t();
    arma::vec beta_k = beta_mat.row(clus_assign[i]).t();
    
    // Find the location of the zero of this particular observations
    
    arma::uvec zijk_zero = arma::find(w == 1 and z.row(i).t() == 0);
    
    if(zijk_zero.size() != 0){
      for(int jj = 0; jj < zijk_zero.size(); ++jj){
        int j = zijk_zero[jj];
        
        // Proposed a new at-risk indicator for (i, j)
        arma::vec gi = gamma_mat.row(i).t();
        int proposed_g_ijk = 1 - gamma_mat(i, j);
        
        // Calculate logA
        double logA = ((-1.0) * log_g_ijk(j, zi, gi, w, beta_k, r0g, r1g)); // Old
        gi.row(j).fill(proposed_g_ijk);
        logA += log_g_ijk(j, zi, gi, w, beta_k, r0g, r1g); // Proposed
        
        // MH
        double logU = std::log(R::runif(0.0, 1.0));
        if(logU <= logA){
          gamma_mat.row(i).col(j).fill(proposed_g_ijk);
        }
        
      }
      
    }
    
  }
  
  return gamma_mat;
  
}

// [[Rcpp::export]]
arma::vec update_w(arma::mat z, arma::uvec clus_assign, arma::mat gamma_mat,
                   arma::vec w, arma::mat beta_mat, double r0w, double r1w){
  
  arma::vec proposed_w(w); 
  
  // Propose a new vector of the important variable indicator
  arma::uvec imp_var = arma::find(w == 1);
  arma::uvec unimp_var = arma::find(w == 0);
  int proposed_index = -1;
  
  if(imp_var.size() <= 1){
    proposed_index = unimp_var[arma::randperm(unimp_var.size())[0]];
    proposed_w[proposed_index] = 1 - proposed_w[proposed_index];
  } else {
    proposed_index = arma::randperm(w.size())[0];
    proposed_w[proposed_index] = 1 - proposed_w[proposed_index];
  }
  
  // Calculate logA
  double logA = log_w(z, clus_assign, gamma_mat, proposed_w, beta_mat, r0w, r1w);
  logA -= log_w(z, clus_assign, gamma_mat, w, beta_mat, r0w, r1w);
  
  // MH
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
  arma::mat MH_var_mat = MH_var * arma::eye(w.size(), w.size());
  
  // Loop through the active cluster
  for(int kk = 0; kk < active_clus.size(); ++kk){
    int k = active_clus[kk];
   
    // Propose a new beta
    arma::vec proposed_beta = arma::mvnrnd(beta_mat.row(k).t(), MH_var_mat);
    
    // Calculate logA
    double logA = 0.0;
    logA += log_beta_k(k, z, clus_assign, gamma_mat, w, proposed_beta, s2);
    logA -= log_beta_k(k, z, clus_assign, gamma_mat, w, beta_mat.row(k).t(), s2);
    
    // MH
    double logU = std::log(R::runif(0.0, 1.0));
    if(logU <= logA){
      beta_mat.row(k) = proposed_beta.t();
    }
    
  }
  
  return beta_mat;
  
}

// Realloc without SM
// [[Rcpp::export]]
Rcpp::List realloc_no_sm(unsigned int K_max, arma::mat z, arma::uvec clus_assign, 
                         arma::mat gamma_mat, arma::vec w, arma::mat beta_mat,
                         arma::vec theta_vec){
  
  // Filter only the important variables
  arma::uvec imp_var = arma::find(w == 1); 
  
  arma::mat z_active = z.cols(imp_var);
  arma::mat gamma_active = gamma_mat.cols(imp_var);
  arma::mat xi_active = arma::exp(beta_mat.cols(imp_var));

  // Count the number of the observation for each cluster
  arma::vec n_clus(K_max, arma::fill::zeros); 
  for(int i = 0; i < clus_assign.size(); ++i){
    n_clus[clus_assign[i]] += 1;
  }
  
  // Calculate the posterior predictive (= marginal distribution)
  // I will create a matrix that store the marginal distribution for all 
  // observation in every clusters.
  
  arma::mat gwx(z_active.n_rows, imp_var.size(), arma::fill::zeros); 
  arma::mat z_gwx(z_active.n_rows, imp_var.size(), arma::fill::zeros);
  arma::vec sum_gwx(z_active.n_rows, arma::fill::zeros);
  arma::vec sum_z_gwx(z_active.n_rows, arma::fill::zeros);
  
  arma::mat log_postpred(z_active.n_rows, K_max, arma::fill::zeros);
  
  for(int k = 0; k < K_max; ++k){
    
    gwx = gamma_active % (arma::repelem(xi_active.row(k).t(), 1, z_active.n_rows).t());
    
    z_gwx = z_active + gwx;
    sum_gwx = arma::sum(gwx, 1);
    sum_z_gwx = arma::sum(z_gwx, 1);
    
    arma::vec log_k(z_active.n_rows, arma::fill::zeros);
    log_k += arma::lgamma(sum_gwx.replace(0.0, 1.0));
    
    log_k -= arma::lgamma(sum_z_gwx.replace(0.0, 1.0));
    log_k += arma::sum(arma::lgamma(z_gwx.replace(0.0, 1.0)), 1);
    log_k -= arma::sum(arma::lgamma(gwx.replace(0.0, 1.0)), 1);
    
    log_postpred.col(k) = log_k;
    
  }

  // Reallocate
  for(int i = 0; i < clus_assign.size(); ++i){
    // Remove the i from the count
    int ci = clus_assign[i];
    n_clus[ci] -= 1;
    
    // Calculate the reallocation probability
    arma::vec log_realloc_prob = log_postpred.row(i).t();
    log_realloc_prob += arma::log(n_clus + theta_vec);
    arma::vec realloc_prob = log_sum_exp(log_realloc_prob);
      
    // Sampling the new cluster assignment
    Rcpp::NumericVector rprob = Rcpp::NumericVector(realloc_prob.begin(), realloc_prob.end());
    Rcpp::IntegerVector newci = rmultinom_1(rprob, K_max);
    arma::uvec new_ci = Rcpp::as<arma::uvec>(Rcpp::wrap(newci));
    arma::uvec index_new_c = arma::find(new_ci == 1);
    
    // Reassign
    clus_assign[i] = index_new_c[0];
    n_clus[index_new_c[0]] += 1;
  }
  
  Rcpp::List result;
  result["clus_assign"] = clus_assign;
  result["z_active"] = z_active;
  result["n_clus"] = n_clus;
  result["log_postpred"] = log_postpred;
  return result;
  
}

// [[Rcpp::export]]
Rcpp::List clus_no_SM(unsigned int iter, unsigned int K_max, arma::mat z, 
                      arma::mat gamma_mat, arma::vec w,
                      double MH_var, double s2, arma::vec theta_vec){
  
  Rcpp::List result;
  
  // Create objects for storing the result.
  arma::cube beta_result(K_max, w.size(), iter);
  arma::mat ci_result(z.n_rows, iter, arma::fill::value(-1));
  
  // Initialize the cluster assignment and beta matrix
  arma::mat beta_init(K_max, w.size(), arma::fill::zeros);
  arma::uvec ci_init(z.n_rows, arma::fill::zeros);
  
  arma::mat beta_mcmc(beta_init);
  
  // Begin
  for(int t = 0; t < iter; ++t){
    
    // Update
    beta_mcmc = update_beta(z, ci_init, gamma_mat, w, beta_init, MH_var, s2);
    Rcpp::List ci_List = realloc_no_sm(K_max, z, ci_init, gamma_mat, w, 
                                       beta_mcmc, theta_vec);
    arma::uvec ci_mcmc = ci_List["clus_assign"];
    
    beta_init = beta_mcmc;
    ci_init = ci_mcmc;
    beta_result.slice(t) = beta_init;
    ci_result.col(t) = arma::conv_to<arma::vec>::from(ci_init);
  }
  
  result["beta"] = beta_result;
  result["ci"] = ci_result.t();
  return result;
}

// [[Rcpp::export]]
Rcpp::List realloc(arma::mat z, arma::uvec clus_assign,
                   arma::mat gamma_mat, arma::vec w, arma::mat beta_mat,
                   arma::vec tau_vec, arma::vec theta_vec){

  arma::uvec active_clus = arma::unique(clus_assign);
  unsigned int K_pos = active_clus.size();

  // Filter only the important variables
  arma::uvec imp_var = arma::find(w == 1);
  arma::mat z_active = z.cols(imp_var);
  arma::mat gamma_active = gamma_mat.cols(imp_var);
  arma::mat xi_active = arma::exp(beta_mat.cols(imp_var));

  // Count the number of the observation for each cluster
  arma::uvec n_clus(K_pos, arma::fill::zeros);
  for(int i = 0; i < clus_assign.size(); ++i){
    arma::uvec index_ci = arma::find(active_clus == clus_assign[i]);
    n_clus[index_ci[0]] += 1;
  }

  // Calculate the posterior predictive (= marginal distribution)
  arma::mat gwx(z_active.n_rows, imp_var.size(), arma::fill::zeros);
  arma::mat z_gwx(z_active.n_rows, imp_var.size(), arma::fill::zeros);
  arma::vec sum_gwx(z_active.n_rows, arma::fill::zeros);
  arma::vec sum_z_gwx(z_active.n_rows, arma::fill::zeros);

  arma::mat log_postpred(z_active.n_rows, K_pos, arma::fill::zeros);

  for(int kk = 0; kk < K_pos; ++kk){

    int k = active_clus[kk];

    gwx = gamma_active % (arma::repelem(xi_active.row(k).t(), 1, z_active.n_rows).t());
    z_gwx = z_active + gwx;
    sum_gwx = arma::sum(gwx, 1);
    sum_z_gwx = arma::sum(z_gwx, 1);

    arma::vec log_k(z_active.n_rows, arma::fill::zeros);
    log_k += arma::lgamma(sum_gwx.replace(0.0, 1.0));
    log_k -= arma::lgamma(sum_z_gwx.replace(0.0, 1.0));
    log_k += arma::sum(arma::lgamma(z_gwx.replace(0.0, 1.0)), 1);
    log_k -= arma::sum(arma::lgamma(gwx.replace(0.0, 1.0)), 1);

    log_postpred.col(kk) = log_k;

  }

  // Reallocate
  for(int i = 0; i < clus_assign.size(); ++i){
    // Remove the i from the count
    int ci = clus_assign[i];
    arma::uvec index_ci_old = arma::find(active_clus == ci);
    n_clus[index_ci_old[0]] -= 1;

    // Calculate the reallocation probability
    arma::vec log_realloc_prob = log_postpred.row(i).t();
    log_realloc_prob += arma::log(n_clus + theta_vec.rows(active_clus));
    arma::vec realloc_prob = log_sum_exp(log_realloc_prob);

    // Sampling the new cluster assignment
    Rcpp::NumericVector rprob = Rcpp::NumericVector(realloc_prob.begin(), realloc_prob.end());
    Rcpp::IntegerVector newci = rmultinom_1(rprob, K_pos);
    arma::uvec new_ci = Rcpp::as<arma::uvec>(Rcpp::wrap(newci));
    arma::uvec index_ci_new = arma::find(new_ci == 1);

    // Reassign
    int ci_new = active_clus[index_ci_new[0]];
    clus_assign[i] = ci_new;
    n_clus[index_ci_new[0]] += 1;
  }

  // Adjust tau and beta
  Rcpp::List updated_tb = adjust_tau_beta(beta_mat, tau_vec, clus_assign);
  arma::vec new_tau = updated_tb["tau"];
  arma::mat new_beta = updated_tb["beta"];
  
  Rcpp::List result;
  result["log_postpred"] = log_postpred;
  result["clus_assign"] = clus_assign;
  result["tau"] = new_tau;
  result["beta"] = new_beta;
  return result;

}

// [[Rcpp::export]]
Rcpp::List sm(unsigned int K_max, arma::mat z, arma::uvec clus_assign,
              arma::mat gamma_mat, arma::vec w, arma::mat beta_mat,
              arma::vec tau_vec, arma::vec theta_vec, unsigned int launch_iter,
              double s2, double r0c, double r1c){

  arma::uvec active_clus = arma::unique(clus_assign);
  unsigned int K_pos = active_clus.size();
  unsigned int n = clus_assign.size();
  int split_ind = -1; // Indicate that we will perform split or merge

  arma::uvec launch_assign(clus_assign);
  arma::vec launch_tau(tau_vec);
  arma::mat launch_beta(beta_mat);

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

    arma::uvec inactive_index = arma::find(tau_vec == 0);
    arma::uvec new_clus_index = arma::randperm(inactive_index.size(), 1);
    new_clus = inactive_index[new_clus_index[0]]; // new active cluster
    samp_clus.row(0).fill(new_clus);

    launch_assign.row(samp_obs[0]).fill(new_clus); // set ci_launch to be a new cluster
    launch_tau.row(new_clus).fill(R::rgamma(theta_vec[new_clus], 1.0));
    
    launch_beta.row(new_clus) = arma::randn(w.size(), arma::distr_param(0.0, std::sqrt(s2))).t();
    
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
    arma::uvec ll = realloc_sm(z, launch_assign, gamma_mat, w, launch_beta, S, samp_clus);
    launch_assign = ll;
  }

  // Split-Merge
  arma::uvec proposed_assign(launch_assign);

  if(split_ind == 1){ // Perform another launch step
    arma::uvec pp = realloc_sm(z, launch_assign, gamma_mat, w, launch_beta, S, samp_clus);
    proposed_assign = pp;
    } else { // Perform Merge
      proposed_assign.rows(S).fill(samp_clus[1]);
      proposed_assign.rows(samp_obs).fill(samp_clus[1]);
   }
    
  Rcpp::List proposed_tb = adjust_tau_beta(launch_beta, launch_tau, proposed_assign);
  arma::vec proposed_tau = proposed_tb["tau"];
  arma::mat proposed_beta = proposed_tb["beta"];
  
  // LogA
  double logA = 0.0;
  
  // Filter only the important variables
  arma::uvec imp_var = arma::find(w == 1);
  arma::mat z_active = z.cols(imp_var);
  arma::mat gamma_active = gamma_mat.cols(imp_var);
  
  arma::mat xi_old_active = arma::exp(beta_mat.cols(imp_var));
  xi_old_active = xi_old_active.rows(clus_assign);
  arma::mat xi_pp_active = arma::exp(proposed_beta.cols(imp_var));
  xi_pp_active = xi_pp_active.rows(proposed_assign);
  
  // Marginal Data: Proposed Cluster Assignment
  arma::mat gwx_pp = gamma_active % xi_pp_active;
  arma::mat z_gwx_pp = z_active + gwx_pp;
  arma::vec sum_gwx_pp = arma::sum(gwx_pp, 1);
  arma::vec sum_z_gwx_pp = arma::sum(z_gwx_pp, 1);
  
  logA += arma::accu(lgamma(sum_gwx_pp.replace(0.0, 1.0)) -
    lgamma(sum_z_gwx_pp.replace(0.0, 1.0)) +
    arma::sum(lgamma(z_gwx_pp.replace(0.0, 1.0)), 1) -
    arma::sum(lgamma(gwx_pp.replace(0.0, 1.0)), 1));
  
  // Marginal Data: Old Cluster Assignment
  arma::mat gwx_old = gamma_active % xi_old_active;
  arma::mat z_gwx_old = z_active + gwx_old;
  arma::vec sum_gwx_old = arma::sum(gwx_old, 1);
  arma::vec sum_z_gwx_old = arma::sum(z_gwx_old, 1);
  
  logA -= arma::accu(lgamma(sum_gwx_old.replace(0.0, 1.0)) -
    lgamma(sum_z_gwx_old.replace(0.0, 1.0)) +
    arma::sum(lgamma(z_gwx_old.replace(0.0, 1.0)), 1) -
    arma::sum(lgamma(gwx_old.replace(0.0, 1.0)), 1));
  
  arma::vec n_old(K_max, arma::fill::zeros); 
  arma::vec n_pp(K_max, arma::fill::zeros);
  
  for(int i = 0; i < clus_assign.size(); ++i){
    n_old[clus_assign[i]] += 1;
    n_pp[proposed_assign[i]] += 1;
  }
  
  // Marginal Cluster: Proposed Cluster Assignment
  arma::vec theta_pp = theta_vec.rows(arma::find(n_pp > 0));
  arma::vec n_theta_pp = n_pp.rows(arma::find(n_pp > 0)) + theta_pp;
  
  logA += (std::lgamma(arma::accu(theta_pp)) - std::lgamma(arma::accu(n_theta_pp)) +
    arma::accu(arma::lgamma(n_theta_pp)) - arma::accu(arma::lgamma(theta_pp)));
  
  // Marginal Cluster: Old Cluster Assignment
  arma::vec theta_old = theta_vec.rows(arma::find(n_old > 0));
  arma::vec n_theta_old = n_old.rows(arma::find(n_old > 0)) + theta_old;
  
  logA -= (std::lgamma(arma::accu(theta_old)) - std::lgamma(arma::accu(n_theta_old)) +
    arma::accu(arma::lgamma(n_theta_old)) - arma::accu(arma::lgamma(theta_old)));
  
  // Beta: Proposed Cluster Assignment
  logA += arma::accu(arma::log_normpdf(proposed_beta, 0.0, std::sqrt(s2)));
  
  // Beta: Proposed Cluster Assignment
  logA -= arma::accu(arma::log_normpdf(beta_mat, 0.0, std::sqrt(s2)));
  
  // Theta
  logA += (((2 * split_ind) - 1) * (std::log(r0c) - std::log(r1c)));
  
  // Proposal
  logA += log_proposal(z, launch_assign, proposed_assign, gamma_mat, w, 
                       launch_beta, S, samp_clus);
  if(split_ind == 1){
    logA -= log_proposal(z, proposed_assign, launch_assign, gamma_mat, w, 
                         launch_beta, S, samp_clus);
  }
  
  // MH
  double logU = std::log(R::runif(0.0, 1.0));
  int accept_MH = 0;
  arma::uvec new_assign(clus_assign);
  arma::vec new_tau(tau_vec);
  arma::mat new_beta(beta_mat);
  
  if(logU <= logA){
    accept_MH += 1;
    new_assign = proposed_assign;
    new_tau = proposed_tau;
    new_beta = proposed_beta;
  }
  
  Rcpp::List result;
  result["split_ind"] = split_ind;
  result["accept_MH"] = accept_MH;
  result["logA"] = logA;
  result["clus_assign"] = new_assign;
  result["tau"] = new_tau;
  result["beta"] = new_beta;
  return result;

}

// [[Rcpp::export]]
Rcpp::List update_tau_u(arma::uvec clus_assign, arma::vec tau_vec, 
                        arma::vec theta_vec, double old_U){

  arma::vec new_tau(tau_vec);

  // update alpha vector
  arma::uvec active_clus = arma::unique(clus_assign);
  for(int k = 0; k < active_clus.size(); ++k){
    int current_c = active_clus[k];
    arma::uvec nk = arma::find(clus_assign == current_c);
    double scale_gamma = 1/(1 + old_U); // change the rate to scale parameter
    new_tau.row(current_c).fill(R::rgamma(nk.size() + theta_vec[current_c], scale_gamma));
  }

  // update U
  int n = clus_assign.size();
  double scale_u = 1/arma::accu(new_tau);
  double new_U = R::rgamma(n, scale_u);

  Rcpp::List result;
  result["tau"] = new_tau;
  result["U"] = new_U;
  return result;

}

// [[Rcpp::export]]
Rcpp::List full_func(unsigned int iter, unsigned int K_max, arma::mat z,
                     arma::mat gamma_mat, arma::vec w,
                     arma::vec theta, unsigned int launch_iter,
                     double MH_var, double s2, double r0c, double r1c){

  // Object for storing the result
  arma::cube beta_result(K_max, w.size(), iter);
  arma::mat ci_result(z.n_rows, iter);
  arma::vec sm_result(iter, arma::fill::zeros);
  arma::vec accept_result(iter, arma::fill::zeros);
  arma::vec logA_sm_result(iter, arma::fill::zeros);

  // Initialize beta and cluster assignment
  arma::mat beta_init(K_max, w.size(), arma::fill::zeros);
  arma::uvec ci_init(z.n_rows, arma::fill::zeros);

  // Initialize tau
  arma::vec tau_init(K_max, arma::fill::zeros);
  tau_init[0] = R::rgamma(theta[0], 1.0);
  double U_init = R::rgamma(z.n_rows, 1/(arma::accu(tau_init)));

  arma::mat beta_mcmc(beta_init);
  Rcpp::List realloc_list;
  Rcpp::List sm_list;
  Rcpp::List param_list;

  // Begin
  for(int t = 0; t < iter; ++t){
    // Update beta
    beta_mcmc = update_beta(z, ci_init, gamma_mat, w, beta_init, MH_var, s2);

    // Update Cluster assignment: (1) Reallocate
    realloc_list = realloc(z, ci_init, gamma_mat, w, beta_mcmc, tau_init, theta);
    arma::uvec ci_realloc = realloc_list["clus_assign"];
    arma::vec tau_realloc = realloc_list["tau"];
    arma::mat beta_realloc = realloc_list["beta"];

    // Update Cluster assignment: (2) SM
    sm_list = sm(K_max, z, ci_realloc, gamma_mat, w, beta_realloc, tau_realloc,
                 theta, launch_iter, s2, r0c, r1c);
    arma::uvec ci_sm = sm_list["clus_assign"];
    arma::vec tau_sm = sm_list["tau"];
    arma::mat beta_sm = sm_list["beta"];
    sm_result.row(t).fill(sm_list["split_ind"]);
    accept_result.row(t).fill(sm_list["accept_MH"]);
    logA_sm_result.row(t).fill(sm_list["logA"]);
    
    // Update Cluster assignment: (3) Update tau and U
    param_list = update_tau_u(ci_sm, tau_sm, theta, U_init);
    double new_U = param_list["U"];
    arma::vec new_tau = param_list["tau"];

    // Record the result
    beta_init = beta_sm;
    ci_init = ci_sm;
    tau_init = new_tau;
    U_init = new_U;
    
    beta_result.slice(t) = beta_init;
    
    for(int i = 0; i < ci_init.size(); ++i){
      ci_result.row(i).col(t).fill(ci_init[i]);
    }
    
  }

  Rcpp::List result;
  result["beta"] = beta_result;
  result["ci"] = ci_result;
  result["accept_result"] = accept_result;
  result["sm_result"] = sm_result;
  result["logA"] = logA_sm_result;
  return result;
  
}

// [[Rcpp::export]]
Rcpp::List full_func_at_risk(unsigned int iter, unsigned int K_max, arma::mat z,
                             arma::vec w, arma::vec theta, unsigned int launch_iter,
                             double MH_var, double s2, double r0g, double r1g, 
                             double r0c, double r1c){
  
  // Object for storing the result
  arma::cube beta_result(K_max, w.size(), iter);
  arma::mat ci_result(z.n_rows, iter);
  arma::vec sm_result(iter, arma::fill::zeros);
  arma::vec accept_result(iter, arma::fill::zeros);
  arma::vec logA_sm_result(iter, arma::fill::zeros);
  
  // Initialize the at-risk indicator
  arma::mat gamma_init(z.n_rows, z.n_cols, arma::fill::ones);
   
  // Initialize beta and cluster assignment
  arma::mat beta_init(K_max, w.size(), arma::fill::zeros);
  arma::uvec ci_init(z.n_rows, arma::fill::zeros);
   
  // Initialize tau
  arma::vec tau_init(K_max, arma::fill::zeros);
  tau_init[0] = R::rgamma(theta[0], 1.0);
  double U_init = R::rgamma(z.n_rows, 1/(arma::accu(tau_init)));
  
  arma::mat gamma_mcmc(gamma_init);
  arma::mat beta_mcmc(beta_init);
  Rcpp::List realloc_list;
  Rcpp::List sm_list;
  Rcpp::List param_list;
   
  // Begin
  for(int t = 0; t < iter; ++t){
    
    // Update at-risk indicator
    gamma_mcmc = update_gamma(z, ci_init, gamma_init, w, beta_init, r0g, r1g);
    
    // Update beta
    beta_mcmc = update_beta(z, ci_init, gamma_mcmc, w, beta_init, MH_var, s2);
     
    // Update Cluster assignment: (1) Reallocate
    realloc_list = realloc(z, ci_init, gamma_mcmc, w, beta_mcmc, tau_init, theta);
    arma::uvec ci_realloc = realloc_list["clus_assign"];
    arma::vec tau_realloc = realloc_list["tau"];
    arma::mat beta_realloc = realloc_list["beta"];
     
    // Update Cluster assignment: (2) SM
    sm_list = sm(K_max, z, ci_realloc, gamma_mcmc, w, beta_realloc, tau_realloc,
                 theta, launch_iter, s2, r0c, r1c);
    arma::uvec ci_sm = sm_list["clus_assign"];
    arma::vec tau_sm = sm_list["tau"];
    arma::mat beta_sm = sm_list["beta"];
    sm_result.row(t).fill(sm_list["split_ind"]);
    accept_result.row(t).fill(sm_list["accept_MH"]);
    logA_sm_result.row(t).fill(sm_list["logA"]);
    
    // Update Cluster assignment: (3) Update tau and U
    param_list = update_tau_u(ci_sm, tau_sm, theta, U_init);
    double new_U = param_list["U"];
    arma::vec new_tau = param_list["tau"];
     
    // Record the result
    gamma_init = gamma_mcmc;
    beta_init = beta_sm;
    ci_init = ci_sm;
    tau_init = new_tau;
    U_init = new_U;
    
    beta_result.slice(t) = beta_init;
    
    for(int i = 0; i < ci_init.size(); ++i){
      ci_result.row(i).col(t).fill(ci_init[i]);
    } 
    
  } 
  
  Rcpp::List result;
  result["beta"] = beta_result;
  result["ci"] = ci_result;
  result["accept_result"] = accept_result;
  result["sm_result"] = sm_result;
  result["logA"] = logA_sm_result;
  return result;
  
}

Rcpp::List full_func_1(unsigned int iter, unsigned int K_max, arma::mat z,
                       arma::mat gamma_mat, arma::vec w,
                       arma::vec theta, unsigned int launch_iter,
                       double MH_var, double s2, double r0c, double r1c){
  
  // Object for storing the result
  arma::cube beta_result(K_max, w.size(), iter);
  arma::mat ci_result(z.n_rows, iter);
  arma::vec sm_result(iter, arma::fill::zeros);
  arma::vec accept_result(iter, arma::fill::zeros);
  arma::vec logA_sm_result(iter, arma::fill::zeros);
  
  // Initial the cluster assignment
  arma::vec n_init(K_max, arma::fill::zeros);
  arma::uvec ci_init(z.n_rows, arma::fill::zeros); 
  for(int i = 0; i < z.n_rows; ++i){
    ci_init[i] = arma::randperm(K_max)[0];
    n_init[ci_init[i]] += 1;
  }
  
  // Initialize tau and beta
  arma::uvec active_init = arma::unique(ci_init);
  arma::vec tau_init(K_max, arma::fill::zeros);
  arma::mat beta_init(K_max, w.size(), arma::fill::zeros);
  
  for(int kk = 0; kk < active_init.size(); ++kk){
    int k = active_init[kk];
    tau_init[k] = R::rgamma(theta[k], 1.0);
    beta_init.row(kk) = arma::randn(w.size(), arma::distr_param(0.0, std::sqrt(s2))).t();
  }
  
  // Initialize U
  double U_init = R::rgamma(z.n_rows, 1/(arma::accu(tau_init)));
  
  arma::mat beta_mcmc(beta_init);
  Rcpp::List realloc_list;
  Rcpp::List sm_list;
  Rcpp::List param_list;
   
  // Begin
  for(int t = 0; t < iter; ++t){
    // Update beta
    beta_mcmc = update_beta(z, ci_init, gamma_mat, w, beta_init, MH_var, s2);
     
    // Update Cluster assignment: (1) Reallocate
    realloc_list = realloc(z, ci_init, gamma_mat, w, beta_mcmc, tau_init, theta);
    arma::uvec ci_realloc = realloc_list["clus_assign"];
    arma::vec tau_realloc = realloc_list["tau"];
    arma::mat beta_realloc = realloc_list["beta"];
     
    // Update Cluster assignment: (2) SM
    sm_list = sm(K_max, z, ci_realloc, gamma_mat, w, beta_realloc, tau_realloc,
                 theta, launch_iter, s2, r0c, r1c);
    arma::uvec ci_sm = sm_list["clus_assign"];
    arma::vec tau_sm = sm_list["tau"];
    arma::mat beta_sm = sm_list["beta"];
    sm_result.row(t).fill(sm_list["split_ind"]);
    accept_result.row(t).fill(sm_list["accept_MH"]);
    logA_sm_result.row(t).fill(sm_list["logA"]);
    
    // Update Cluster assignment: (3) Update tau and U
    param_list = update_tau_u(ci_sm, tau_sm, theta, U_init);
    double new_U = param_list["U"];
    arma::vec new_tau = param_list["tau"];
     
    // Record the result
    beta_init = beta_sm;
    ci_init = ci_sm;
    tau_init = new_tau;
    U_init = new_U;
    
    beta_result.slice(t) = beta_init;
    
    for(int i = 0; i < ci_init.size(); ++i){
      ci_result.row(i).col(t).fill(ci_init[i]);
    } 
    
  } 
  
  Rcpp::List result;
  result["beta"] = beta_result;
  result["ci"] = ci_result;
  result["accept_result"] = accept_result;
  result["sm_result"] = sm_result;
  result["logA"] = logA_sm_result;
  result["init_clus"] = active_init.size();
  return result;
  
} 
