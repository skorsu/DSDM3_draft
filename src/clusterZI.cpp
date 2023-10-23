#include "RcppArmadillo.h"
#define pi 3.141592653589793238462643383280

// [[Rcpp::depends(RcppArmadillo)]]

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
  
  double log_prob = log_marginal(zi, gamma_ik, beta_k);
  
  arma::vec g0 = gamma_ik.rows(arma::find(zi == 0));
  log_prob += arma::accu(arma::lgamma(g0 + r0g) + arma::lgamma((1 - g0) + r1g));
  
  return log_prob;
  
}