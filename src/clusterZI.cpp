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

