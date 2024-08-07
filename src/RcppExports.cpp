// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rmultinom_1
arma::vec rmultinom_1(arma::vec& probs_arma);
RcppExport SEXP _ClusterZI_rmultinom_1(SEXP probs_armaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type probs_arma(probs_armaSEXP);
    rcpp_result_gen = Rcpp::wrap(rmultinom_1(probs_arma));
    return rcpp_result_gen;
END_RCPP
}
// log_sum_exp
arma::vec log_sum_exp(arma::vec log_unnorm_prob);
RcppExport SEXP _ClusterZI_log_sum_exp(SEXP log_unnorm_probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type log_unnorm_prob(log_unnorm_probSEXP);
    rcpp_result_gen = Rcpp::wrap(log_sum_exp(log_unnorm_prob));
    return rcpp_result_gen;
END_RCPP
}
// log_marginal
double log_marginal(arma::vec zi, arma::vec gmi, arma::vec beta_k);
RcppExport SEXP _ClusterZI_log_marginal(SEXP ziSEXP, SEXP gmiSEXP, SEXP beta_kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type zi(ziSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gmi(gmiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_k(beta_kSEXP);
    rcpp_result_gen = Rcpp::wrap(log_marginal(zi, gmi, beta_k));
    return rcpp_result_gen;
END_RCPP
}
// logmar_ik
double logmar_ik(arma::rowvec zi, arma::rowvec atrisk_i, arma::rowvec beta_k);
RcppExport SEXP _ClusterZI_logmar_ik(SEXP ziSEXP, SEXP atrisk_iSEXP, SEXP beta_kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type zi(ziSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type atrisk_i(atrisk_iSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type beta_k(beta_kSEXP);
    rcpp_result_gen = Rcpp::wrap(logmar_ik(zi, atrisk_i, beta_k));
    return rcpp_result_gen;
END_RCPP
}
// logmar_k
arma::vec logmar_k(arma::mat z, arma::mat atrisk, arma::rowvec beta_k);
RcppExport SEXP _ClusterZI_logmar_k(SEXP zSEXP, SEXP atriskSEXP, SEXP beta_kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type atrisk(atriskSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type beta_k(beta_kSEXP);
    rcpp_result_gen = Rcpp::wrap(logmar_k(z, atrisk, beta_k));
    return rcpp_result_gen;
END_RCPP
}
// logmar
arma::mat logmar(arma::mat z, arma::mat atrisk, arma::mat beta_mat);
RcppExport SEXP _ClusterZI_logmar(SEXP zSEXP, SEXP atriskSEXP, SEXP beta_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type atrisk(atriskSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_mat(beta_matSEXP);
    rcpp_result_gen = Rcpp::wrap(logmar(z, atrisk, beta_mat));
    return rcpp_result_gen;
END_RCPP
}
// log_atrisk
double log_atrisk(arma::rowvec zi, arma::rowvec atrisk_i, arma::rowvec beta_k, double r0g, double r1g);
RcppExport SEXP _ClusterZI_log_atrisk(SEXP ziSEXP, SEXP atrisk_iSEXP, SEXP beta_kSEXP, SEXP r0gSEXP, SEXP r1gSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type zi(ziSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type atrisk_i(atrisk_iSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type beta_k(beta_kSEXP);
    Rcpp::traits::input_parameter< double >::type r0g(r0gSEXP);
    Rcpp::traits::input_parameter< double >::type r1g(r1gSEXP);
    rcpp_result_gen = Rcpp::wrap(log_atrisk(zi, atrisk_i, beta_k, r0g, r1g));
    return rcpp_result_gen;
END_RCPP
}
// launch_mcmc
Rcpp::List launch_mcmc(arma::mat z, arma::mat atrisk, arma::mat beta_mat, arma::uvec ci_old, unsigned int launch_iter, arma::uvec S, arma::uvec samp_clus, arma::vec nk);
RcppExport SEXP _ClusterZI_launch_mcmc(SEXP zSEXP, SEXP atriskSEXP, SEXP beta_matSEXP, SEXP ci_oldSEXP, SEXP launch_iterSEXP, SEXP SSEXP, SEXP samp_clusSEXP, SEXP nkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type atrisk(atriskSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_mat(beta_matSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ci_old(ci_oldSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type launch_iter(launch_iterSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type samp_clus(samp_clusSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nk(nkSEXP);
    rcpp_result_gen = Rcpp::wrap(launch_mcmc(z, atrisk, beta_mat, ci_old, launch_iter, S, samp_clus, nk));
    return rcpp_result_gen;
END_RCPP
}
// log_proposal
double log_proposal(arma::mat z, arma::mat atrisk, arma::mat beta_mat, arma::uvec ci_A, arma::uvec ci_B, arma::uvec S, arma::uvec samp_clus);
RcppExport SEXP _ClusterZI_log_proposal(SEXP zSEXP, SEXP atriskSEXP, SEXP beta_matSEXP, SEXP ci_ASEXP, SEXP ci_BSEXP, SEXP SSEXP, SEXP samp_clusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type atrisk(atriskSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_mat(beta_matSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ci_A(ci_ASEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ci_B(ci_BSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type samp_clus(samp_clusSEXP);
    rcpp_result_gen = Rcpp::wrap(log_proposal(z, atrisk, beta_mat, ci_A, ci_B, S, samp_clus));
    return rcpp_result_gen;
END_RCPP
}
// split_beta
arma::mat split_beta(unsigned int nbeta, unsigned int clus_split, unsigned int clus_new, arma::mat beta_old, double mu, double s2);
RcppExport SEXP _ClusterZI_split_beta(SEXP nbetaSEXP, SEXP clus_splitSEXP, SEXP clus_newSEXP, SEXP beta_oldSEXP, SEXP muSEXP, SEXP s2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type nbeta(nbetaSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type clus_split(clus_splitSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type clus_new(clus_newSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_old(beta_oldSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    rcpp_result_gen = Rcpp::wrap(split_beta(nbeta, clus_split, clus_new, beta_old, mu, s2));
    return rcpp_result_gen;
END_RCPP
}
// update_atrisk
arma::mat update_atrisk(arma::mat z, arma::mat atrisk_old, arma::mat beta_mat, arma::uvec ci, double r0g, double r1g);
RcppExport SEXP _ClusterZI_update_atrisk(SEXP zSEXP, SEXP atrisk_oldSEXP, SEXP beta_matSEXP, SEXP ciSEXP, SEXP r0gSEXP, SEXP r1gSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type atrisk_old(atrisk_oldSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_mat(beta_matSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ci(ciSEXP);
    Rcpp::traits::input_parameter< double >::type r0g(r0gSEXP);
    Rcpp::traits::input_parameter< double >::type r1g(r1gSEXP);
    rcpp_result_gen = Rcpp::wrap(update_atrisk(z, atrisk_old, beta_mat, ci, r0g, r1g));
    return rcpp_result_gen;
END_RCPP
}
// update_beta
Rcpp::List update_beta(arma::mat z, arma::mat atrisk, arma::mat beta_old, arma::uvec ci, double mu, double s2, double s2_MH);
RcppExport SEXP _ClusterZI_update_beta(SEXP zSEXP, SEXP atriskSEXP, SEXP beta_oldSEXP, SEXP ciSEXP, SEXP muSEXP, SEXP s2SEXP, SEXP s2_MHSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type atrisk(atriskSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_old(beta_oldSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ci(ciSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< double >::type s2_MH(s2_MHSEXP);
    rcpp_result_gen = Rcpp::wrap(update_beta(z, atrisk, beta_old, ci, mu, s2, s2_MH));
    return rcpp_result_gen;
END_RCPP
}
// update_beta_adaptive
Rcpp::List update_beta_adaptive(unsigned int t, unsigned int t_threshold, arma::mat z, arma::mat atrisk, arma::cube beta_record, arma::mat beta_init, arma::uvec ci, double mu, double s2, double s2_MH);
RcppExport SEXP _ClusterZI_update_beta_adaptive(SEXP tSEXP, SEXP t_thresholdSEXP, SEXP zSEXP, SEXP atriskSEXP, SEXP beta_recordSEXP, SEXP beta_initSEXP, SEXP ciSEXP, SEXP muSEXP, SEXP s2SEXP, SEXP s2_MHSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type t(tSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type t_threshold(t_thresholdSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type atrisk(atriskSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type beta_record(beta_recordSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ci(ciSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< double >::type s2_MH(s2_MHSEXP);
    rcpp_result_gen = Rcpp::wrap(update_beta_adaptive(t, t_threshold, z, atrisk, beta_record, beta_init, ci, mu, s2, s2_MH));
    return rcpp_result_gen;
END_RCPP
}
// realloc_full
arma::uvec realloc_full(unsigned int Kmax, arma::mat z, arma::mat atrisk, arma::mat beta_mat, arma::uvec ci_old, double theta);
RcppExport SEXP _ClusterZI_realloc_full(SEXP KmaxSEXP, SEXP zSEXP, SEXP atriskSEXP, SEXP beta_matSEXP, SEXP ci_oldSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type Kmax(KmaxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type atrisk(atriskSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_mat(beta_matSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ci_old(ci_oldSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(realloc_full(Kmax, z, atrisk, beta_mat, ci_old, theta));
    return rcpp_result_gen;
END_RCPP
}
// sm
Rcpp::List sm(unsigned int Kmax, unsigned int nbeta_split, arma::mat z, arma::mat atrisk, arma::mat beta_old, arma::uvec ci_old, double theta, double mu, double s2, unsigned int launch_iter, double r0c, double r1c);
RcppExport SEXP _ClusterZI_sm(SEXP KmaxSEXP, SEXP nbeta_splitSEXP, SEXP zSEXP, SEXP atriskSEXP, SEXP beta_oldSEXP, SEXP ci_oldSEXP, SEXP thetaSEXP, SEXP muSEXP, SEXP s2SEXP, SEXP launch_iterSEXP, SEXP r0cSEXP, SEXP r1cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type Kmax(KmaxSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nbeta_split(nbeta_splitSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type atrisk(atriskSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_old(beta_oldSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ci_old(ci_oldSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< unsigned int >::type launch_iter(launch_iterSEXP);
    Rcpp::traits::input_parameter< double >::type r0c(r0cSEXP);
    Rcpp::traits::input_parameter< double >::type r1c(r1cSEXP);
    rcpp_result_gen = Rcpp::wrap(sm(Kmax, nbeta_split, z, atrisk, beta_old, ci_old, theta, mu, s2, launch_iter, r0c, r1c));
    return rcpp_result_gen;
END_RCPP
}
// mod
Rcpp::List mod(unsigned int iter, unsigned int Kmax, unsigned int nbeta_split, arma::mat z, arma::mat atrisk_init, arma::mat beta_init, arma::uvec ci_init, double theta, double mu, double s2, double s2_MH, unsigned int launch_iter, double r0g, double r1g, double r0c, double r1c, unsigned int thin);
RcppExport SEXP _ClusterZI_mod(SEXP iterSEXP, SEXP KmaxSEXP, SEXP nbeta_splitSEXP, SEXP zSEXP, SEXP atrisk_initSEXP, SEXP beta_initSEXP, SEXP ci_initSEXP, SEXP thetaSEXP, SEXP muSEXP, SEXP s2SEXP, SEXP s2_MHSEXP, SEXP launch_iterSEXP, SEXP r0gSEXP, SEXP r1gSEXP, SEXP r0cSEXP, SEXP r1cSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type Kmax(KmaxSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nbeta_split(nbeta_splitSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type atrisk_init(atrisk_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ci_init(ci_initSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< double >::type s2_MH(s2_MHSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type launch_iter(launch_iterSEXP);
    Rcpp::traits::input_parameter< double >::type r0g(r0gSEXP);
    Rcpp::traits::input_parameter< double >::type r1g(r1gSEXP);
    Rcpp::traits::input_parameter< double >::type r0c(r0cSEXP);
    Rcpp::traits::input_parameter< double >::type r1c(r1cSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(mod(iter, Kmax, nbeta_split, z, atrisk_init, beta_init, ci_init, theta, mu, s2, s2_MH, launch_iter, r0g, r1g, r0c, r1c, thin));
    return rcpp_result_gen;
END_RCPP
}
// mod_adaptive
Rcpp::List mod_adaptive(unsigned int iter, unsigned int Kmax, unsigned int nbeta_split, arma::mat z, arma::mat atrisk_init, arma::mat beta_init, arma::uvec ci_init, double theta, double mu, double s2, double s2_MH, unsigned int t_thres, unsigned int launch_iter, double r0g, double r1g, double r0c, double r1c, unsigned int thin);
RcppExport SEXP _ClusterZI_mod_adaptive(SEXP iterSEXP, SEXP KmaxSEXP, SEXP nbeta_splitSEXP, SEXP zSEXP, SEXP atrisk_initSEXP, SEXP beta_initSEXP, SEXP ci_initSEXP, SEXP thetaSEXP, SEXP muSEXP, SEXP s2SEXP, SEXP s2_MHSEXP, SEXP t_thresSEXP, SEXP launch_iterSEXP, SEXP r0gSEXP, SEXP r1gSEXP, SEXP r0cSEXP, SEXP r1cSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type Kmax(KmaxSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nbeta_split(nbeta_splitSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type atrisk_init(atrisk_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ci_init(ci_initSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< double >::type s2_MH(s2_MHSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type t_thres(t_thresSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type launch_iter(launch_iterSEXP);
    Rcpp::traits::input_parameter< double >::type r0g(r0gSEXP);
    Rcpp::traits::input_parameter< double >::type r1g(r1gSEXP);
    Rcpp::traits::input_parameter< double >::type r0c(r0cSEXP);
    Rcpp::traits::input_parameter< double >::type r1c(r1cSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(mod_adaptive(iter, Kmax, nbeta_split, z, atrisk_init, beta_init, ci_init, theta, mu, s2, s2_MH, t_thres, launch_iter, r0g, r1g, r0c, r1c, thin));
    return rcpp_result_gen;
END_RCPP
}
// realloc_dm
arma::uvec realloc_dm(unsigned int Kmax, arma::mat z, arma::mat atrisk, arma::mat beta_mat, arma::uvec ci_old, double theta);
RcppExport SEXP _ClusterZI_realloc_dm(SEXP KmaxSEXP, SEXP zSEXP, SEXP atriskSEXP, SEXP beta_matSEXP, SEXP ci_oldSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type Kmax(KmaxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type atrisk(atriskSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_mat(beta_matSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ci_old(ci_oldSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(realloc_dm(Kmax, z, atrisk, beta_mat, ci_old, theta));
    return rcpp_result_gen;
END_RCPP
}
// DM_DM
Rcpp::List DM_DM(unsigned int iter, unsigned int Kmax, arma::mat z, arma::mat beta_init, arma::uvec ci_init, double theta, double mu, double s2, double s2_MH, unsigned int t_thres, unsigned int thin);
RcppExport SEXP _ClusterZI_DM_DM(SEXP iterSEXP, SEXP KmaxSEXP, SEXP zSEXP, SEXP beta_initSEXP, SEXP ci_initSEXP, SEXP thetaSEXP, SEXP muSEXP, SEXP s2SEXP, SEXP s2_MHSEXP, SEXP t_thresSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type Kmax(KmaxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ci_init(ci_initSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< double >::type s2_MH(s2_MHSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type t_thres(t_thresSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(DM_DM(iter, Kmax, z, beta_init, ci_init, theta, mu, s2, s2_MH, t_thres, thin));
    return rcpp_result_gen;
END_RCPP
}
// DM_ZIDM
Rcpp::List DM_ZIDM(unsigned int iter, unsigned int Kmax, unsigned int nbeta_split, arma::mat z, arma::mat beta_init, arma::uvec ci_init, double theta, double mu, double s2, double s2_MH, unsigned int t_thres, unsigned int launch_iter, double r0c, double r1c, unsigned int thin);
RcppExport SEXP _ClusterZI_DM_ZIDM(SEXP iterSEXP, SEXP KmaxSEXP, SEXP nbeta_splitSEXP, SEXP zSEXP, SEXP beta_initSEXP, SEXP ci_initSEXP, SEXP thetaSEXP, SEXP muSEXP, SEXP s2SEXP, SEXP s2_MHSEXP, SEXP t_thresSEXP, SEXP launch_iterSEXP, SEXP r0cSEXP, SEXP r1cSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type Kmax(KmaxSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nbeta_split(nbeta_splitSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ci_init(ci_initSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< double >::type s2_MH(s2_MHSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type t_thres(t_thresSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type launch_iter(launch_iterSEXP);
    Rcpp::traits::input_parameter< double >::type r0c(r0cSEXP);
    Rcpp::traits::input_parameter< double >::type r1c(r1cSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(DM_ZIDM(iter, Kmax, nbeta_split, z, beta_init, ci_init, theta, mu, s2, s2_MH, t_thres, launch_iter, r0c, r1c, thin));
    return rcpp_result_gen;
END_RCPP
}
// ZIDM_DM
Rcpp::List ZIDM_DM(unsigned int iter, unsigned int Kmax, arma::mat z, arma::mat atrisk_init, arma::mat beta_init, arma::uvec ci_init, double theta, double mu, double s2, double s2_MH, unsigned int t_thres, double r0g, double r1g, unsigned int thin);
RcppExport SEXP _ClusterZI_ZIDM_DM(SEXP iterSEXP, SEXP KmaxSEXP, SEXP zSEXP, SEXP atrisk_initSEXP, SEXP beta_initSEXP, SEXP ci_initSEXP, SEXP thetaSEXP, SEXP muSEXP, SEXP s2SEXP, SEXP s2_MHSEXP, SEXP t_thresSEXP, SEXP r0gSEXP, SEXP r1gSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type Kmax(KmaxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type atrisk_init(atrisk_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ci_init(ci_initSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< double >::type s2_MH(s2_MHSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type t_thres(t_thresSEXP);
    Rcpp::traits::input_parameter< double >::type r0g(r0gSEXP);
    Rcpp::traits::input_parameter< double >::type r1g(r1gSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(ZIDM_DM(iter, Kmax, z, atrisk_init, beta_init, ci_init, theta, mu, s2, s2_MH, t_thres, r0g, r1g, thin));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _ClusterZI_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _ClusterZI_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _ClusterZI_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _ClusterZI_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ClusterZI_rmultinom_1", (DL_FUNC) &_ClusterZI_rmultinom_1, 1},
    {"_ClusterZI_log_sum_exp", (DL_FUNC) &_ClusterZI_log_sum_exp, 1},
    {"_ClusterZI_log_marginal", (DL_FUNC) &_ClusterZI_log_marginal, 3},
    {"_ClusterZI_logmar_ik", (DL_FUNC) &_ClusterZI_logmar_ik, 3},
    {"_ClusterZI_logmar_k", (DL_FUNC) &_ClusterZI_logmar_k, 3},
    {"_ClusterZI_logmar", (DL_FUNC) &_ClusterZI_logmar, 3},
    {"_ClusterZI_log_atrisk", (DL_FUNC) &_ClusterZI_log_atrisk, 5},
    {"_ClusterZI_launch_mcmc", (DL_FUNC) &_ClusterZI_launch_mcmc, 8},
    {"_ClusterZI_log_proposal", (DL_FUNC) &_ClusterZI_log_proposal, 7},
    {"_ClusterZI_split_beta", (DL_FUNC) &_ClusterZI_split_beta, 6},
    {"_ClusterZI_update_atrisk", (DL_FUNC) &_ClusterZI_update_atrisk, 6},
    {"_ClusterZI_update_beta", (DL_FUNC) &_ClusterZI_update_beta, 7},
    {"_ClusterZI_update_beta_adaptive", (DL_FUNC) &_ClusterZI_update_beta_adaptive, 10},
    {"_ClusterZI_realloc_full", (DL_FUNC) &_ClusterZI_realloc_full, 6},
    {"_ClusterZI_sm", (DL_FUNC) &_ClusterZI_sm, 12},
    {"_ClusterZI_mod", (DL_FUNC) &_ClusterZI_mod, 17},
    {"_ClusterZI_mod_adaptive", (DL_FUNC) &_ClusterZI_mod_adaptive, 18},
    {"_ClusterZI_realloc_dm", (DL_FUNC) &_ClusterZI_realloc_dm, 6},
    {"_ClusterZI_DM_DM", (DL_FUNC) &_ClusterZI_DM_DM, 11},
    {"_ClusterZI_DM_ZIDM", (DL_FUNC) &_ClusterZI_DM_ZIDM, 15},
    {"_ClusterZI_ZIDM_DM", (DL_FUNC) &_ClusterZI_ZIDM_DM, 14},
    {"_ClusterZI_rcpparma_hello_world", (DL_FUNC) &_ClusterZI_rcpparma_hello_world, 0},
    {"_ClusterZI_rcpparma_outerproduct", (DL_FUNC) &_ClusterZI_rcpparma_outerproduct, 1},
    {"_ClusterZI_rcpparma_innerproduct", (DL_FUNC) &_ClusterZI_rcpparma_innerproduct, 1},
    {"_ClusterZI_rcpparma_bothproducts", (DL_FUNC) &_ClusterZI_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_ClusterZI(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
