#include <RcppArmadillo.h>
#include "BayesGgmHelper.h"
#include "progress.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace Rcpp;

//' Bayesian adaptive graphical ridge block Gibbs sampler for Gaussian graphical models.
//'
//' Implements a Bayesian adaptive graphical ridge block Gibbs sampler to simulate the
//' posterior distribution of the precision matrix for Gaussian graphical models.
//'
//' @param X A numeric matrix, assumed to be generated from a multivariate Gaussian distribution.
//' @param burnin An integer specifying the number of burn-in iterations.
//' @param iterations An integer specifying the length of the Markov chain after the burn-in iterations.
//' @param a A double specifying the value of the shape parameter for the inverse gamma prior.
//' @param b A double specifying the value of the scale parameter for the inverse gamma prior.
//' @param verbose A logical determining whether the progress of the MCMC sampler should be displayed.
//' @return A list containing precision `Omega` and covariance `Sigma` matrices
//' from the Markov chains.
//' @examples
//'# Generate true precision matrix:
//'p             <- 10
//'n             <- 500
//' OmegaTrue    <- pracma::Toeplitz(c(0.7^rep(1:p-1)))
//' SigTrue      <- pracma::inv(OmegaTrue)
//'# Generate expected value vector:
//'mu            <- rep(0,p)
//'# Generate multivariate normal distribution:
//'set.seed(123)
//'X             <- MASS::mvrnorm(n, mu = mu, Sigma = SigTrue)
//'# Generate posterior distribution:
//'posterior     <- blockBAGR(X, iterations = 1000, burnin = 500)
//'# Estimated precision matrix using the mean of the posterior:
//'OmegaEst      <- apply(simplify2array(posterior$Omega), 1:2, mean)
//' @export
// [[Rcpp::export]]
List blockBAGR(arma::mat X, int burnin, int iterations, bool verbose = true, double a = 1, double b = 1e-2){
  // variable declarations and initialisations
  int totIter, n, p;
  totIter = burnin + iterations;
  n = X.n_rows;
  p = X.n_cols;

  Progress prog(totIter, verbose);

  arma::mat S(p, p), Sig(p,p), C(p, p), ind_noi_all(p-1, p), tau(p, p, arma::fill::zeros);
  S = X.t()* X;
  arma::mat filler(p, p, arma::fill::value(1e-12));
  Sig = arma::eye<arma::mat>(p, p) + filler;
  C = arma::eye<arma::mat>(p, p) + filler;
  arma::mat indmx(p,p);
  int idx_counter = 0;
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < p; j++) {
      indmx(j, i) = idx_counter;
      idx_counter += 1;
    }
  }
  arma::vec upperind = extract_upperind(indmx);
  arma::vec lowerind = extract_upperind(indmx.t());
  ind_noi_all.fill(NA_REAL);
  arma::vec tau_temp(p-1);
  arma::mat Sig11(p-1, p-1);
  arma::vec Sig12(p-1);
  arma::vec S_ind(p-1);
  arma::mat invC11(p-1, p-1);
  arma::mat Ci(p-1, p-1);
  arma::mat Ci_chol(p-1, p-1);
  arma::vec mu_i(p-1);
  arma::vec beta(p-1);
  List SigmaMatList(iterations);
  List OmegaMatList(iterations);
  arma::vec permInt(p);
  double sig_ii = 0.556533162 -p*0.0042253;

  for(int i = 0; i < p; i++){
    permInt(i) = i;
  }

  for (int i = 0; i < p; i++) {
    if (i == 0){
      ind_noi_all.col(i) = permInt.subvec(i+1,p-1);//;
    }
    if (i == p-1){
      ind_noi_all.col(i) = permInt.subvec(0,p-2);
    }
    if ((i != p-1) & (i != 0)){
      arma::vec A, B;
      A = permInt.subvec(0,i-1);
      B = permInt.subvec(i+1,p-1);
      ind_noi_all.col(i) = join_vert(A, B);
    }
  }

  int idx = 0;
  for (int iter=0; iter<totIter; iter++){
    if (Progress::check_abort())
      return -1.0;

    if (verbose){
      prog.increment(); // update progress
    }

    arma::vec low_tri_tau_vals(upperind.size());

    for (int i=0; i<low_tri_tau_vals.size(); i++){
      double a_post = 0.5 + a;
      double b_post = std::pow(C(upperind(i)), 2)/2.0 + b;
      low_tri_tau_vals(i)= 1.0/rrgamma(1, a_post, 1.0/b_post)(0);
    }

    for(int i = 0; i < upperind.size(); i++)
    {
      tau(upperind(i)) = low_tri_tau_vals(i);
    }
    for(int i = 0; i < lowerind.size(); i++)
    {
      tau(lowerind(i)) = low_tri_tau_vals(i);
    }

    for (int i=0; i<p; i++){
      arma::uvec ind_noi = arma::conv_to<arma::uvec>::from(ind_noi_all.col(i));
      for(int idx=0; idx<ind_noi.size(); idx++){
        int j = ind_noi(idx);
        tau_temp(idx) = tau(j,i);
        Sig12(idx) = Sig(j,i);
        S_ind(idx) = S(j,i);
      }
      for(int idx1=0; idx1<ind_noi.size(); idx1++){
        int i = ind_noi(idx1);
        for(int idx2=0; idx2<ind_noi.size(); idx2++){
          int j = ind_noi(idx2);
          Sig11(idx1,idx2) = Sig(i,j);
        }
      }
      invC11 = Sig11 - (Sig12*Sig12.t())/Sig(i,i);
      Ci = ((S(i,i)*std::pow(sig_ii, 2)+1.0)/std::pow(sig_ii, 2))*invC11 + diagmat(1/tau_temp);
      Ci = (Ci + Ci.t())/2;
      Ci_chol = chol(Ci) + arma::eye<arma::mat>(p-1, p-1)*1e-8;
      mu_i = -inv(Ci)*S_ind;
      arma::vec rnorm = Rcpp::rnorm(p-1, 0, 1);
      beta = mu_i + inv(Ci_chol)*rnorm;
      for(int idx=0; idx<ind_noi.size(); idx++){
        int j = ind_noi(idx);
        C(j,i) = beta(idx);
        C(i,j) = beta(idx);
      }
      NumericVector gam = rrgamma(1, n/2 + 1.0, 1.0/((((S(i,i)*std::pow(sig_ii, 2)+1.0)/std::pow(sig_ii, 2)))/2.0));
      C(i,i) = gam(0) + arma::as_scalar(beta.t() * invC11 * beta);
      Sig = inv(C);
    }

    if (iter>=burnin){
      OmegaMatList[idx] = C;
      SigmaMatList[idx] = Sig;
      idx += 1;
    }
  }
  // Create the final result list
  return Rcpp::List::create(Rcpp::Named("Omega") = OmegaMatList,
                            Rcpp::Named("Sigma") = SigmaMatList);
}


