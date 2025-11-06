#ifndef UTILS_H
#define UTILS_H

arma::mat psinv(const arma::mat A, const int max_iter, const double tol);
double abs3(double x);
int WhichMax(arma::vec vX);
int rando_index(arma::vec p);
arma::vec zero_bound(arma::vec parvec, double zbound);
Rcpp::List logisticReg(arma::vec vY, arma::mat vZ, arma::vec SmoothProb,
                 arma::vec vBeta, int maxIter, double tol);
Rcpp::List NR_step_cov(arma::mat mX, arma::mat mbeta, arma::mat mU);
Rcpp::List NR_step_covIT(arma::mat mX, arma::mat mbeta, arma::mat mU, double tol, int maxIt);
Rcpp::List NR_step_covIT_LS(arma::mat mX, arma::mat mbeta, arma::mat mU, double dC, double tol, int maxIt);
Rcpp::List NR_step_covIT_wei(arma::mat mX, arma::mat mbeta, arma::mat mU, arma::vec vWei, double tol, int maxIt);
Rcpp::List NR_step_covIT_wei_LS(arma::mat mX, arma::mat mbeta, arma::mat mU, double dC, arma::vec vWei, double tol, int maxIt);
arma::mat vecTomatClass(arma::vec vClass);
Rcpp::List AvgMarginalEff(arma::mat beta, arma::mat P_ij, arma::vec weights);
arma::vec grad_MLTLCA(arma::vec parvec, arma::vec vY, arma::vec vD, arma::vec vPW_N, arma::mat mPX, arma::mat mPMX);
#endif
