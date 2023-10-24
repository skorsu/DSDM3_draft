#include "RcppArmadillo.h"
#define pi 3.141592653589793238462643383280

// [[Rcpp::depends(RcppArmadillo)]]

// Log Probability: ------------------------------------------------------------

// [[Rcpp::export]]
double log_marginal(arma::vec zi, arma::vec gamma_ik, arma::vec beta_k){
  /* Calculate the marginal distribution of the data in a log scale */
  
  double log_prob = 0.0;

  arma::vec zijk_active = zi % gamma_ik;
  arma::vec xi_active = arma::exp(beta_k) % gamma_ik;
  
  arma::vec zx_active = zijk_active + xi_active;
  zx_active = zx_active.rows(arma::find(zx_active != 0));
  
  log_prob += std::lgamma(arma::accu(zijk_active) + 1);
  log_prob -= arma::accu(arma::lgamma(zijk_active + 1));
  log_prob += arma::accu(arma::lgamma(zx_active));
  log_prob -= std::lgamma(arma::accu(zx_active));

  return log_prob;
  
}

// [[Rcpp::export]]
double log_atrisk(arma::vec zi, arma::vec gamma_ik, arma::vec beta_k,
                  double r0g, double r1g){
  /* Calculate the probability of the at-risk vector */
  
  double log_prob = log_marginal(zi, gamma_ik, beta_k);
  
  arma::vec g0 = gamma_ik.rows(arma::find(zi == 0));
  log_prob += arma::accu(arma::lgamma(g0 + r0g) + arma::lgamma((1 - g0) + r1g));
  
  return log_prob;
  
}

// [[Rcpp::export]]
double log_betak(int k, arma::mat z, arma::mat gamma_mat, arma::vec beta_k,
                 arma::uvec ci, double mu_beta, double s2_beta){
  
  double log_prob = 0.0;
  
  arma::mat zk = z.rows(arma::find(ci == k));
  arma::mat gk = gamma_mat.rows(arma::find(ci == k));
  
  for(int i = 0; i < zk.n_rows; ++i){
    log_prob += log_marginal(zk.row(i).t(), gk.row(i).t(), beta_k);
  }
  
  log_prob += arma::accu(arma::log_normpdf(beta_k, mu_beta, std::sqrt(s2_beta)));
  
  return log_prob;
  
}

// Updating Parameters: --------------------------------------------------------
// [[Rcpp::export]]
arma::cube update_atrisk(arma::uvec ci, arma::mat z, arma::cube gamma_old, 
                         arma::mat beta_mat, double r0g, double r1g){
  /* Update the at-risk indicator for all zero across all active clusters */
  arma::cube gamma_new(gamma_old);
  arma::uvec active_clus = arma::unique(ci);
  arma::umat z0 = z == 0;
  
  for(int kk = 0; kk < active_clus.size(); ++kk){ // Loop through active cluster
    arma::mat gmat = gamma_old.slice(active_clus[kk]);
    arma::vec bk = beta_mat.row(active_clus[kk]).t(); 
    
    for(int i = 0; i < z.n_rows; ++i){ // Loop through observation
      arma::vec zi = z.row(i).t(); 
      arma::uvec z0i = arma::find(z0.row(i) == 1);
      arma::vec gik_old = gmat.row(i).t();

      for(int jj = 0; jj < z0i.size(); ++jj){ // Loop through the variable
        // Proposed at-risk
        arma::vec gik_proposed(gik_old);
        gik_proposed[z0i[jj]] = 1 - gik_old[z0i[jj]];
        // Calculate log A
        double logA = log_atrisk(zi, gik_proposed, bk, r0g, r1g) -
          log_atrisk(zi, gik_old, bk, r0g, r1g);
        // Decide to accept?
        double logU = std::log(R::runif(0.0, 1.0));
        if(logU <= logA){
          gik_old = gik_proposed;
        }
      }

      gmat.row(i) = gik_old.t();
      
    }
    
    gamma_new.slice(active_clus[kk]) = gmat;
    
  }
  
  return gamma_new;
  
}

// [[Rcpp::export]]
arma::mat update_beta(arma::uvec ci, arma::mat z, arma::cube gamma_cube, 
                      arma::mat beta_old, double mu_beta, double s2_beta,
                      double s2_MH){
  
  arma::mat beta_new(beta_old);
  arma::uvec active_clus = arma::unique(ci);
  arma::mat s2_MH_mat = std::sqrt(s2_MH) * arma::eye(beta_old.n_cols, beta_old.n_cols); 
  
  for(int kk = 0; kk < active_clus.size(); ++kk){
    // Proposed a new beta
    arma::vec betak_old = beta_old.row(active_clus[kk]).t();
    arma::vec betak_proposed = arma::mvnrnd(betak_old, s2_MH_mat);
    // Calculate the log A
    double logA = log_betak(active_clus[kk], z, gamma_cube.slice(active_clus[kk]), 
                            betak_proposed, ci, mu_beta, s2_beta) - 
                  log_betak(active_clus[kk], z, gamma_cube.slice(active_clus[kk]), 
                            betak_old, ci, mu_beta, s2_beta);
    // Decide to accept?
    double logU = std::log(R::runif(0.0, 1.0));
    if(logU <= logA){
      beta_new.row(active_clus[kk]) = betak_proposed.t();
    }
  }
  
  return beta_new;
  
}

