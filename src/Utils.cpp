#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include "SafeFunctions.h"

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat psinv(const arma::mat A, int max_iter = 1000, double tol = 2.220446e-16) {
  int k;
  arma::mat AH = A.t();          // Hermitian transpose
  AH = arma::conj(AH);
  arma::mat absB = arma::abs(A * AH);
  
  // Compute step size alpha
  arma::vec row_sums = arma::sum(absB, 1);
  double m = row_sums.max();
  double alpha = 1.0 / m;           // any alpha < 2/m works
  
  // Initialize Y
  // arma::cx_mat Y = alpha * arma::conj(A.t());
  arma::mat Y = alpha * AH;
  arma::mat Ynext;
  
  // Precompute identity
  arma::mat I = arma::eye(A.n_rows, A.n_rows);
  
  for (k = 0; k < max_iter; ++k) {
    Ynext = Y * (2.0 * I - A * Y);
    double diff = arma::abs(Ynext - Y).max();
    if (diff < tol) {
      // Rcpp::Rcout << "Converged after " << k + 1 << " steps" << std::endl;
      break;
    }
    Y = Ynext;
  }
  
  return Ynext;
}







double abs3(double x){
  double abs_x = x;
  if (abs_x < 0) {
    abs_x = -abs_x;
  }
  return abs_x;
}

int WhichMax(arma::vec vX){
  
  int iK = vX.size();
  int k;
  
  int iMax = 0;
  double dMax = vX(0);
  
  for (k = 1; k < iK; k++) {
    if (vX(k) > dMax) {
      iMax = k;
      dMax = vX(k);
    }
  }
  
  return iMax;
}


//[[Rcpp::export]]
int rando_index(arma::vec p) {
  int N = p.size();
  
  double u = R::runif(0, 1.0);
  int    i = 0;
  double s = p(0);
  
  while((u > s) && (i < N - 1)) {
    i++;
    s += p(i);
  }
  
  return i;
}


//[[Rcpp::export]]
arma::vec zero_bound(arma::vec parvec, double zbound = 1e-04){
  int npar = parvec.n_elem;
  arma::vec trunc_parvec(npar);
  for(int n = 0; n < npar; n++){
    if(parvec(n) < zbound && parvec(n) > -zbound){
      trunc_parvec(n) = 0.0;
    }
    else{
      trunc_parvec(n) = parvec(n);
    }
  }
  return trunc_parvec;
}


//[[Rcpp::export]]
List logisticReg(arma::vec vY, arma::mat vZ, arma::vec SmoothProb,
                 arma::vec vBeta, int maxIter = 5, double tol = 1e-5){
  // needs the vectorized smoothed probabilities
  // Remember to include a column of ones in vZ to include the intercept
  int i,k;
  int iK  = vZ.n_cols;// the CFactor vector + a vector of ones for overall intercept + vZ
  int iN  = vY.n_elem;
  int iter = 0;
  double eps = 1e100;
  arma::mat vI=zeros(iK,iK);
  arma::vec vS(iK);
  arma::vec vW=zeros(iN);
  arma::vec vProb=zeros(iN);
  arma::vec vSEs=zeros(iK);
  arma::vec vEta(iN);
  arma::vec foosbeta(iK);
  
  arma::vec vBetaup(iK);
  arma::vec vBetaold = vBeta;
  // vBeta.zeros();
  // vBetaold.ones();
  arma::mat vH(iK,iK);
  double foostep = 0.0;
  double lk0;
  double lk;
  while(eps > tol && iter<maxIter){
    vS.fill(0.0);
    vH.fill(0.0);
    if(iter == 0){
      vEta = exp(vZ * vBetaold);
      vProb = vEta/(1.0 + vEta);
      vW    = vProb % (1.0 - vProb);
    }
    for(i = 0; i < iN; i++){
      vS = vS + SmoothProb(i) *vZ.row(i).t() * (vY(i) - vProb(i));
      vH = vH + SmoothProb(i)*vProb(i)*(1-vProb(i))*vZ.row(i).t() * vZ.row(i);
    }
    // vH = -vH;
    vBetaup = solve(vH,vS);
    for(k = 0; k < iK; k++){
      foosbeta(k) = abs3(vBetaup(k));
    }
    foostep = max(foosbeta);
    if(foostep > 0.5){
      vBetaup = vBetaup/foostep*0.5;
    }
    vBeta = vBetaold + vBetaup;
    vEta = exp(vZ * vBeta);
    vProb = vEta/(1.0 + vEta);
    vW    = vProb % (1.0 - vProb);
    lk = accu(SmoothProb%vProb);
    vBetaold = vBeta;
    lk0 = lk;
    if(iter > 3){
      eps = abs3(lk -lk0);
    }
    iter += 1;
  }
  
  vI = -inv(vH);
  vSEs = sqrt(diagvec(vI));
  
  List logisticRegout;
  
  logisticRegout["vBeta"]   = vBeta;
  logisticRegout["vSEs"]   = vSEs;
  logisticRegout["vProb"]   = vProb;
  logisticRegout["vW"]   = vW;
  logisticRegout["vI"]   = vI;
  logisticRegout["vS"]   = vS;
  logisticRegout["niter"]   = iter;
  
  return logisticRegout;
}


//[[Rcpp::export]]
List NR_step_cov(arma::mat mX, arma::mat mbeta, arma::mat mU){
  // X should contain a column of ones to include the intercept term
  // mbeta is K-1 x P
  // N is the sample size
  // K is the number of classes
  // P is the number of covariates (including the intercept)
  int N = mX.n_rows;
  int K = mU.n_cols;
  int P = mX.n_cols;
  int j,k,n,p;
  arma::mat X = mX;
  arma::mat w = ones(N,K);
  arma::mat w_i = w;
  arma::mat ta = w;
  arma::cube ibeta(P,P,K-1);
  arma::mat sbeta(K-1,P);
  arma::vec foosbeta(P);
  arma::vec NRstep(P);
  double foostep;
  for(n = 0; n < N; n++){
    for(k = 1; k < K; k++){
      w(n,k) = exp(accu(X.row(n) % mbeta.row(k-1)));
      w_i(n,k) = w(n,k);
      ta(n,k) = w(n,k);
    }
    w_i.row(n) = w_i.row(n)/accu(w_i.row(n));
    ta.row(n) = w.row(n)/pow(accu(w.row(n)),2.0);
  }
  double lk0 = accu(mU%log(w_i));
  for(k = 1; k < K; k++){
    sbeta.row(k-1) = (mU.col(k) - w_i.col(k)).t() * X;
  }
  if(K > 2){
    arma::vec ww(N);
    for(k=1; k < K; k++){
      ww.fill(0.0);
      for(j = 0; j < K; j++){
        if(k!=j){
          ww = ww + w.col(j);
        }
      }
      ibeta.slice(k-1) = X.t() * (repmat(ta.col(k),1,P)%X%repmat(ww,1,P));
      NRstep = solve(ibeta.slice(k-1),sbeta.row(k-1).t());
      for(p = 0; p < P; p++){
        foosbeta(p) = abs3(NRstep(p));
      }
      foostep = max(foosbeta);
      if(foostep > 0.5){
        NRstep = NRstep/foostep*0.5;
      }
      mbeta.row(k-1) = mbeta.row(k-1) + NRstep.t();
    }
  }
  if(K==2){
    double ww;
    ww = accu(w.col(1));
    for(k=1; k < K; k++){
      ibeta.slice(k-1) = X.t() * (repmat(ta.col(k),1,P)%X*ww);
      NRstep = solve(ibeta.slice(k-1),sbeta.row(k-1).t());
      for(p = 0; p < P; p++){
        foosbeta(p) = abs3(NRstep(p));
      }
      foostep = max(foosbeta);
      if(foostep > 0.5){
        NRstep = NRstep/foostep*0.5;
      }
      mbeta.row(k-1) = mbeta.row(k-1) + NRstep.t();
    }
  }
  w.fill(1.0);
  for(n = 0; n < N; n++){
    for(k = 1; k < K; k++){
      w(n,k) = exp(accu(X.row(n) % mbeta.row(k-1)));
    }
    w.row(n) = w.row(n)/accu(w.row(n));
  }
  
  double lk = accu(mU%log(w));
  List NR_out;
  NR_out["beta"] = mbeta;
  NR_out["ibeta"] = ibeta;
  NR_out["sbeta"] = sbeta;
  NR_out["w"] = w;
  NR_out["lk"] = lk;
  NR_out["lk0"] = lk0;
  
  return NR_out;
}


//[[Rcpp::export]]
List NR_step_covIT(arma::mat mX, arma::mat mbeta, arma::mat mU, double tol=1e-06, int maxIt = 100){
  // X should contain a column of ones to include the intercept term
  // mbeta is K-1 x P
  // N is the sample size
  // K is the number of classes
  // P is the number of covariates (including the intercept)
  int N = mX.n_rows;
  int K = mU.n_cols;
  int P = mX.n_cols;
  int j,k,n,p;
  arma::mat X = mX;
  arma::mat w = ones(N,K);
  arma::mat w_i = w;
  arma::mat ta = w;
  arma::cube ibeta(P,P,K-1);
  arma::mat sbeta(K-1,P);
  arma::vec foosbeta(P);
  arma::vec NRstep(P);
  double foostep;
  double eps=1e100;
  double lk0;
  double lk;
  int it=0;
  for(n = 0; n < N; n++){
    for(k = 1; k < K; k++){
      w(n,k) = exp(accu(X.row(n) % mbeta.row(k-1)));
      w_i(n,k) = w(n,k);
      ta(n,k) = w(n,k);
    }
    w_i.row(n) = w_i.row(n)/accu(w_i.row(n));
    ta.row(n) = w.row(n)/pow(accu(w.row(n)),2.0);
  }
  lk0 = accu(mU%log(w_i));
  while(eps > tol && it < maxIt){
    for(k = 1; k < K; k++){
      sbeta.row(k-1) = (mU.col(k) - w_i.col(k)).t() * X;
    }
    if(K > 2){
      arma::vec ww(N);
      for(k=1; k < K; k++){
        ww.fill(0.0);
        for(j = 0; j < K; j++){
          if(k!=j){
            ww = ww + w.col(j);
          }
        }
        ibeta.slice(k-1) = X.t() * (repmat(ta.col(k),1,P)%X%repmat(ww,1,P));
        NRstep = solve(ibeta.slice(k-1),sbeta.row(k-1).t());
        for(p = 0; p < P; p++){
          foosbeta(p) = abs3(NRstep(p));
        }
        foostep = max(foosbeta);
        if(foostep > 0.5){
          NRstep = NRstep/foostep*0.5;
        }
        // needs a trick to make sure it's monotone!!!
        // code below taken from lasso paper
        // begins here!!!!
        // for(c = 0; c < iC; c++){
        //   ctrl = 1;
        //   while(ctrl==1){
        //     alphanext.col(c) = alphastart - NR_con(c)*NR_step;
        //     for(g = 0; g < (G-1);g++){
        //       if(abs3(alphanext(g,c)) > 12.0){
        //         NR_con(c) = NR_con(c)*0.9;
        //       }
        //       else{
        //         ctrl = 0;
        //       }
        //     }
        //   }
        //   for(g = 0; g < (G-1);g++){
        //     pignext(g,c) = exp(alphanext(g,c));
        //   }
        //   pignext.col(c) = pignext.col(c)/accu(pignext.col(c));
        //   par_loglike(c) = (-1.0)*accu(Uplusl%log(pignext.col(c)) - pow(pignext.col(c),gamma)%lambdag%normBg);
        //   
        // }
        // ends here!!!!
        
        mbeta.row(k-1) = mbeta.row(k-1) + NRstep.t();
      }
    }
    if(K==2){
      double ww;
      ww = accu(w.col(1));
      for(k=1; k < K; k++){
        ibeta.slice(k-1) = X.t() * (repmat(ta.col(k),1,P)%X*ww);
        NRstep = solve(ibeta.slice(k-1),sbeta.row(k-1).t());
        for(p = 0; p < P; p++){
          foosbeta(p) = abs3(NRstep(p));
        }
        foostep = max(foosbeta);
        if(foostep > 0.5){
          NRstep = NRstep/foostep*0.5;
        }
        mbeta.row(k-1) = mbeta.row(k-1) + NRstep.t();
      }
    }
    w.fill(1.0);
    w_i.fill(1.0);
    ta.fill(1.0);
    for(n = 0; n < N; n++){
      for(k = 1; k < K; k++){
        w(n,k) = exp(accu(X.row(n) % mbeta.row(k-1)));
        w_i(n,k) = w(n,k);
        ta(n,k) = w(n,k);
      }
      w_i.row(n) = w_i.row(n)/accu(w_i.row(n));
      ta.row(n) = w.row(n)/pow(accu(w.row(n)),2.0);
    }
    
    lk = accu(mU%log(w_i));
    if(it > 3){
      eps = abs3(lk -lk0);
    }
    it = it + 1;
    lk0 = lk;
  }
  
  arma::mat mSbeta(N,(K-1)*P);
  int iter = 0;
  for(k = 1; k < K; k++){
    for(n = 0; n < N; n++){
      mSbeta.row(n).subvec(iter,iter + P - 1) = (mU(n,k) - w_i(n,k)) * X.row(n);
    }
    iter = iter + P;
  }
  iter = 0;
  
  List NR_out;
  NR_out["beta"]   = mbeta;
  NR_out["ibeta"]  = ibeta;
  NR_out["sbeta"]  = sbeta;
  NR_out["w_i"]    = w_i;
  NR_out["mSbeta"] = mSbeta;
  NR_out["lk"] = lk;
  NR_out["lk0"] = lk0;
  
  return NR_out;
}


//[[Rcpp::export]]
List NR_step_covIT_LS(arma::mat mX, arma::mat mbeta, arma::mat mU, double dC, double tol=1e-06, int maxIt = 100){
  // LS stands for line search: preserves monotonicity in the EM algorithm
  // iC is an integer, tells how much the NR step should be decreased to recompute the parameter update
  // X should contain a column of ones to include the intercept term
  // mbeta is K-1 x P
  // N is the sample size
  // K is the number of classes
  // P is the number of covariates (including the intercept)
  int N = mX.n_rows;
  int K = mU.n_cols;
  int P = mX.n_cols;
  int j,k,n,p;
  double NR_LS = 1.0;
  arma::mat mbeta_foo = mbeta;
  arma::mat X = mX;
  arma::mat w = ones(N,K);
  arma::mat w_i = w;
  arma::mat w_i_foo = w_i;
  arma::mat w_i_foonum = w_i;
  arma::mat ta = w;
  arma::cube ibeta(P,P,K-1);
  arma::mat sbeta(K-1,P);
  arma::vec foosbeta(P);
  arma::vec NRstep(P);
  double foostep;
  double eps=1e100;
  double lk0;
  double lk;
  // arma::vec lkfoo = zeros(iC,1);
  double dLKfoo = 0.0;
  double control = 0.0;
  int it=0;
  int fooiter =0;
  for(n = 0; n < N; n++){
    for(k = 1; k < K; k++){
      w(n,k) = exp(accu(X.row(n) % mbeta.row(k-1)));
      w_i(n,k) = w(n,k);
      ta(n,k) = w(n,k);
    }
    w_i.row(n) = w_i.row(n)/accu(w_i.row(n));
    ta.row(n) = w.row(n)/pow(accu(w.row(n)),2.0);
  }
  lk0 = accu(mU%log(w_i));
  while(eps > tol && it < maxIt){
    for(k = 1; k < K; k++){
      sbeta.row(k-1) = (mU.col(k) - w_i.col(k)).t() * X;
    }
    if(K > 2){
      arma::vec ww(N);
      for(k=1; k < K; k++){
        ww.fill(0.0);
        for(j = 0; j < K; j++){
          if(k!=j){
            ww = ww + w.col(j);
          }
        }
        ibeta.slice(k-1) = X.t() * (repmat(ta.col(k),1,P)%X%repmat(ww,1,P));
        NRstep = solve(ibeta.slice(k-1),sbeta.row(k-1).t());
        for(p = 0; p < P; p++){
          foosbeta(p) = abs3(NRstep(p));
        }
        foostep = max(foosbeta);
        if(foostep > 0.5){
          NRstep = NRstep/foostep*0.5;
        }
        // for(c = 0; c < iC; c++){
        //   cbeta.slice(c).row(k-1) = mbeta.row(k-1) + NR_LS(c)*NRstep.t();
        // }
        // mbeta.row(k-1) = mbeta.row(k-1) + NRstep.t();
      }
    }
    if(K==2){
      double ww;
      ww = accu(w.col(1));
      for(k=1; k < K; k++){
        ibeta.slice(k-1) = X.t() * (repmat(ta.col(k),1,P)%X*ww);
        NRstep = solve(ibeta.slice(k-1),sbeta.row(k-1).t());
        for(p = 0; p < P; p++){
          foosbeta(p) = abs3(NRstep(p));
        }
        foostep = max(foosbeta);
        if(foostep > 0.5){
          NRstep = NRstep/foostep*0.5;
        }
        // mbeta.row(k-1) = mbeta.row(k-1) + NRstep.t();
        // for(c = 0; c < iC; c++){
        //   cbeta.slice(c).row(k-1) = mbeta.row(k-1) + NR_LS(c)*NRstep.t();
        // }
        
      }
    }
    NR_LS = 1.0;
    control = 0.0;
    fooiter = 0;
    while(control == 0.0 && fooiter < maxIt){
      w_i_foo.fill(1.0);
      for(k=1; k < K; k++){
        mbeta_foo.row(k-1) = mbeta.row(k-1) + NR_LS*NRstep.t();
      }
      w_i_foo.cols(1,K-1) = exp(X*mbeta_foo.t());
      w_i_foonum = w_i_foo;
      for(n = 0; n < N; n++){
        w_i_foo.row(n) = w_i_foo.row(n)/accu(w_i_foo.row(n));
      }
      dLKfoo = accu(mU%log(w_i_foo));
      if((dLKfoo -lk0) < 0){
        NR_LS = NR_LS*dC;
        fooiter +=1;
      }else{
        control = 1.0;
        mbeta = mbeta_foo;
      }
    }
    w   = w_i_foonum;
    w_i = w_i_foo;
    
    
    lk = dLKfoo;
    
    for(n = 0; n < N; n++){
      ta.row(n) = w.row(n)/pow(accu(w.row(n)),2.0);
    }
    
    eps = abs3(lk -lk0);
    it = it + 1;
    lk0 = lk;
  }
  
  arma::mat mSbeta(N,(K-1)*P);
  int iter = 0;
  for(k = 1; k < K; k++){
    for(n = 0; n < N; n++){
      mSbeta.row(n).subvec(iter,iter + P - 1) = (mU(n,k) - w_i(n,k)) * X.row(n);
    }
    iter = iter + P;
  }
  iter = 0;
  
  List NR_out;
  NR_out["beta"]   = mbeta;
  NR_out["ibeta"]  = ibeta;
  NR_out["sbeta"]  = sbeta;
  NR_out["w_i"]    = w_i;
  NR_out["mSbeta"] = mSbeta;
  NR_out["lk"] = lk;
  NR_out["lk0"] = lk0;
  
  return NR_out;
}



//[[Rcpp::export]]
List NR_step_covIT_wei(arma::mat mX, arma::mat mbeta, arma::mat mU, arma::vec vWei, double tol=1e-06, int maxIt = 100){
  // X should contain a column of ones to include the intercept term
  // mbeta is K-1 x P
  // N is the sample size
  // K is the number of classes
  // P is the number of covariates (including the intercept)
  int N = mX.n_rows;
  int K = mU.n_cols;
  int P = mX.n_cols;
  int j,k,n,p;
  arma::mat X = mX;
  arma::mat w = ones(N,K);
  arma::mat w_i = w;
  arma::mat ta = w;
  arma::cube ibeta(P,P,K-1);
  arma::mat sbeta(K-1,P);
  arma::vec foosbeta(P);
  arma::vec NRstep(P);
  double foostep;
  double eps=1e100;
  double lk0;
  double lk;
  int it=0;
  for(n = 0; n < N; n++){
    for(k = 1; k < K; k++){
      w(n,k) = exp(accu(X.row(n) % mbeta.row(k-1)));
      w_i(n,k) = w(n,k);
      ta(n,k) = w(n,k);
    }
    w_i.row(n) = w_i.row(n)/accu(w_i.row(n));
    ta.row(n) = w.row(n)/pow(accu(w.row(n)),2.0);
  }
  lk0 = accu(mU%log(w_i));
  arma::mat mWei = repmat(vWei,1,P);
  while(eps > tol && it < maxIt){
    for(k = 1; k < K; k++){
      sbeta.row(k-1) = (mU.col(k) - w_i.col(k)).t() *(mWei % X);
    }
    if(K > 2){
      arma::vec ww(N);
      for(k=1; k < K; k++){
        ww.fill(0.0);
        for(j = 0; j < K; j++){
          if(k!=j){
            ww = ww + w.col(j);
          }
        }
        ibeta.slice(k-1) = X.t() * (repmat(ta.col(k),1,P)%X%repmat(ww,1,P)%mWei);
        NRstep = solve(ibeta.slice(k-1),sbeta.row(k-1).t());
        for(p = 0; p < P; p++){
          foosbeta(p) = abs3(NRstep(p));
        }
        foostep = max(foosbeta);
        if(foostep > 0.5){
          NRstep = NRstep/foostep*0.5;
        }
        mbeta.row(k-1) = mbeta.row(k-1) + NRstep.t();
      }
    }
    if(K==2){
      double ww;
      ww = accu(w.col(1));
      for(k=1; k < K; k++){
        ibeta.slice(k-1) = X.t() * (repmat(ta.col(k),1,P)%X%mWei*ww);
        NRstep = solve(ibeta.slice(k-1),sbeta.row(k-1).t());
        for(p = 0; p < P; p++){
          foosbeta(p) = abs3(NRstep(p));
        }
        foostep = max(foosbeta);
        if(foostep > 0.5){
          NRstep = NRstep/foostep*0.5;
        }
        mbeta.row(k-1) = mbeta.row(k-1) + NRstep.t();
      }
    }
    w.fill(1.0);
    w_i.fill(1.0);
    ta.fill(1.0);
    for(n = 0; n < N; n++){
      for(k = 1; k < K; k++){
        w(n,k) = exp(accu(X.row(n) % mbeta.row(k-1)));
        w_i(n,k) = w(n,k);
        ta(n,k) = w(n,k);
      }
      w_i.row(n) = w_i.row(n)/accu(w_i.row(n));
      ta.row(n) = w.row(n)/pow(accu(w.row(n)),2.0);
    }
    
    lk = accu(mU%log(w_i));
    if(it > 3){
      eps = abs3(lk -lk0);
    }
    it = it + 1;
    lk0 = lk;
  }
  
  arma::mat mSbeta(N,(K-1)*P);
  int iter = 0;
  for(k = 1; k < K; k++){
    for(n = 0; n < N; n++){
      mSbeta.row(n).subvec(iter,iter + P - 1) = (mU(n,k) - w_i(n,k)) * (vWei(n) * X.row(n));
    }
    iter = iter + P;
  }
  iter = 0;
  
  List NR_out;
  NR_out["beta"] = mbeta;
  NR_out["ibeta"] = ibeta;
  NR_out["sbeta"] = sbeta;
  NR_out["mSbeta"] = mSbeta;
  NR_out["w_i"] = w_i;
  NR_out["lk"] = lk;
  NR_out["lk0"] = lk0;
  NR_out["niter"] = it;
  
  return NR_out;
}


//[[Rcpp::export]]
List NR_step_covIT_wei_LS(arma::mat mX, arma::mat mbeta, arma::mat mU, double dC, arma::vec vWei, double tol=1e-06, int maxIt = 100){
  // X should contain a column of ones to include the intercept term
  // mbeta is K-1 x P
  // N is the sample size
  // K is the number of classes
  // P is the number of covariates (including the intercept)
  int N = mX.n_rows;
  int K = mU.n_cols;
  int P = mX.n_cols;
  int j,k,n,p;
  arma::mat X = mX;
  arma::mat w = ones(N,K);
  arma::mat w_i = w;
  arma::mat w_i_foo = w_i;
  arma::mat w_i_foonum = w_i;
  arma::mat ta = w;
  arma::cube ibeta(P,P,K-1);
  arma::mat sbeta(K-1,P);
  arma::mat mbeta_foo = mbeta;
  arma::vec foosbeta(P);
  arma::vec NRstep(P);
  double foostep;
  double NR_LS = 1.0;
  double control = 0.0;
  double eps=1e100;
  double lk0;
  double lk;
  double dLKfoo;
  int it=0;
  int  fooiter = 0;
  for(n = 0; n < N; n++){
    for(k = 1; k < K; k++){
      w(n,k) = exp(accu(X.row(n) % mbeta.row(k-1)));
      w_i(n,k) = w(n,k);
      ta(n,k) = w(n,k);
    }
    w_i.row(n) = w_i.row(n)/accu(w_i.row(n));
    ta.row(n) = w.row(n)/pow(accu(w.row(n)),2.0);
  }
  lk0 = accu(mU%log(w_i));
  arma::mat mWei = repmat(vWei,1,P);
  while(eps > tol && it < maxIt){
    for(k = 1; k < K; k++){
      sbeta.row(k-1) = (mU.col(k) - w_i.col(k)).t() *(mWei % X);
    }
    if(K > 2){
      arma::vec ww(N);
      for(k=1; k < K; k++){
        ww.fill(0.0);
        for(j = 0; j < K; j++){
          if(k!=j){
            ww = ww + w.col(j);
          }
        }
        ibeta.slice(k-1) = X.t() * (repmat(ta.col(k),1,P)%X%repmat(ww,1,P)%mWei);
        NRstep = solve(ibeta.slice(k-1),sbeta.row(k-1).t());
        for(p = 0; p < P; p++){
          foosbeta(p) = abs3(NRstep(p));
        }
        foostep = max(foosbeta);
        if(foostep > 0.5){
          NRstep = NRstep/foostep*0.5;
        }
        // mbeta.row(k-1) = mbeta.row(k-1) + NRstep.t();
      }
    }
    if(K==2){
      double ww;
      ww = accu(w.col(1));
      for(k=1; k < K; k++){
        ibeta.slice(k-1) = X.t() * (repmat(ta.col(k),1,P)%X%mWei*ww);
        NRstep = solve(ibeta.slice(k-1),sbeta.row(k-1).t());
        for(p = 0; p < P; p++){
          foosbeta(p) = abs3(NRstep(p));
        }
        foostep = max(foosbeta);
        if(foostep > 0.5){
          NRstep = NRstep/foostep*0.5;
        }
        // mbeta.row(k-1) = mbeta.row(k-1) + NRstep.t();
      }
    }
    // 
    NR_LS = 1.0;
    control = 0.0;
    fooiter = 0;
    while(control == 0.0 && fooiter < maxIt){
      w_i_foo.fill(1.0);
      for(k=1; k < K; k++){
        mbeta_foo.row(k-1) = mbeta.row(k-1) + NR_LS*NRstep.t();
      }
      w_i_foo.cols(1,K-1) = exp(X*mbeta_foo.t());
      w_i_foonum = w_i_foo;
      for(n = 0; n < N; n++){
        w_i_foo.row(n) = w_i_foo.row(n)/accu(w_i_foo.row(n));
      }
      dLKfoo = accu(mU%log(w_i_foo));
      if((dLKfoo -lk0) < 0){
        NR_LS = NR_LS*dC;
        fooiter +=1;
      }else{
        control = 1.0;
        mbeta = mbeta_foo;
      }
    }
    w   = w_i_foonum;
    w_i = w_i_foo;
    // 
    lk = dLKfoo;
    // 
    for(n = 0; n < N; n++){
      ta.row(n) = w.row(n)/pow(accu(w.row(n)),2.0);
    }
    // 
    eps = abs3(lk -lk0);
    it = it + 1;
    lk0 = lk;
  }
  
  arma::mat mSbeta(N,(K-1)*P);
  int iter = 0;
  for(k = 1; k < K; k++){
    for(n = 0; n < N; n++){
      mSbeta.row(n).subvec(iter,iter + P - 1) = (mU(n,k) - w_i(n,k)) * (vWei(n) * X.row(n));
    }
    iter = iter + P;
  }
  iter = 0;
  
  List NR_out;
  NR_out["beta"] = mbeta;
  NR_out["ibeta"] = ibeta;
  NR_out["sbeta"] = sbeta;
  NR_out["mSbeta"] = mSbeta;
  NR_out["w_i"] = w_i;
  NR_out["lk"] = lk;
  NR_out["lk0"] = lk0;
  NR_out["niter"] = it;
  
  return NR_out;
}






// [[Rcpp::export]]
arma::mat vecTomatClass(arma::vec vClass){
  int n;
  int iT = max(vClass);
  int iN = vClass.n_elem;
  vClass = vClass - 1.0;
  arma::mat mClass = zeros(iN,iT);
  for(n = 0; n < iN; n++){
    mClass(n,vClass(n)) = 1.0;
  }
  return mClass;
}


// [[Rcpp::export]]
List AvgMarginalEff(arma::mat beta, arma::mat P_ij, arma::vec weights){
  int J = P_ij.n_cols;
  int N = P_ij.n_rows;
  int K = beta.n_rows;
  int k, j, n;
  arma::cube cDelta = zeros(N,J,K);
  arma::mat betafull = zeros(K,J);
  betafull.cols(1,J-1) = beta;
  double pijfoo;
  for(k = 0; k < K; k++){
    for(n = 0; n < N; n++){
      pijfoo = accu(P_ij.row(n)%betafull.row(k));
      for(j = 0; j < J; j++){
        cDelta(n,j,k) = weights(n)*P_ij(n,j)*(betafull(k,j) - pijfoo);
      }
    }
  }
  arma::mat delta = mean(cDelta,0);
  List Out;
  Out["cDelta"] = cDelta;
  Out["avgMGeff"] = delta;
  return Out;
}

// [[Rcpp::export]]
arma::vec grad_MLTLCA(arma::vec parvec, arma::vec vY, arma::vec vD, arma::vec vPW_N, arma::mat mPX, arma::mat mPMX){
  // vD is a vector of ones if there's no missing on the items
  int k,m,t;
  int iM = mPX.n_cols;
  int iT = mPX.n_rows;
  int iK = vY.n_elem;
  int iP = 1;
  arma::mat mGamma(iT-1,iM);
  arma::mat mPhi(iK,iT);
  arma::vec vOmega = ones(iM,1);
  int foopar = 0;
  vOmega.subvec(1,iM-1) = exp(parvec.subvec(0,iM-2));
  vOmega = vOmega/accu(vOmega);
  foopar = foopar + iM-1;
  for(m = 0; m < iM; m++){
    for(t =1; t< iT; t++){
      mGamma(t-1,m) = parvec(foopar);
      foopar = foopar +1;
    }
  }
  mPhi = reshape(conv_to<mat>::from(parvec(span(foopar,foopar + iK*iT -1))),iK,iT);
  mPhi = exp(mPhi)/(1.0 + exp(mPhi));
  arma::mat mPi     = ones(iT,iM);
  for(m = 0; m < iM; m++){
    for(t = 1; t < iT; t++){
      mPi(t,m) = exp(mGamma(t-1,m));
    }
    mPi.col(m) = mPi.col(m)/accu(mPi.col(m));
  }
  arma::vec vOmega_Score(iM-1);
  for(m =1; m< iM; m++){
    vOmega_Score(m-1) = vPW_N(m)*(1.0 - vOmega(m));
  }
  // 
  arma::vec vGamma_Score((iT-1)*iP*iM);
  int iFoo2 = 0;
  arma::vec vGammaScore_foo = zeros((iT-1)*iP,1);
  for(m = 0; m < iM; m++){
    vGammaScore_foo.fill(0.0);
    for(t = 1; t < iT; t++){
      vGammaScore_foo(t-1) = (mPX(t,m) - mPi(t,m))*(vPW_N(m));
    }
    
    vGamma_Score.subvec(iFoo2, iFoo2 + ((iT-1)*iP) -1) = vGammaScore_foo;
    iFoo2 = iFoo2 + ((iT-1)*iP);
  }
  iFoo2 = 0;
  
  arma::vec vBeta_Score=zeros(iK*iT,1);
  int iroll = 0;
  for(t = 0; t < iT; t++){
    for(k = 0; k < iK; k++){
      if(vD(k)==1){
        for(m = 0; m < iM; m++){
          vBeta_Score(iroll) += mPMX(t,m)*(vY(k) - mPhi(k,t));
        }
        iroll += 1;
      }else{
        for(m = 0; m < iM; m++){
          vBeta_Score(iroll) += 0.0;
        }
        iroll += 1;
      }
    }
  }
  
  
  arma::mat vScore = join_cols(join_cols(vOmega_Score,vGamma_Score),vBeta_Score);
  
  return vScore;
}


// [[Rcpp::export]]
arma::vec grad_MLTLCA_cov(arma::vec parvec, arma::mat mPhi, arma::vec vY, arma::vec vZ, arma::vec vD, arma::vec vPW_N, arma::mat mPX, arma::mat mPMX,
                          arma::vec vPMsumX,
                          arma::ivec ivItemcat, int nstep = 1){
  // vD is a vector of ones if there's no missing on the items
  int l,k,m,p,t,v;
  int iV    = ivItemcat.n_elem;
  int iM = mPX.n_cols;
  int iT = mPX.n_rows;
  int iK = vY.n_elem;
  int iP = vZ.n_elem;
  int foopar = 0;
  arma::cube cGamma(iT-1,iP,iM);
  // arma::mat mPhi(iK,iT);
  arma::vec vOmega = ones(iM,1);
  vOmega.subvec(1,iM-1) = exp(parvec.subvec(0,iM-2));
  vOmega = vOmega/accu(vOmega);
  foopar = foopar + iM-1;
  for(m = 0; m < iM; m++){
    for(p = 0; p < iP; p++){
      for(t =1; t< iT; t++){
        cGamma(t-1,p,m) = parvec(foopar);
        foopar = foopar +1;
      }
    }
  }
  // mPhi = reshape(conv_to<mat>::from(parvec(span(foopar,foopar + iK*iT -1))),iK,iT);
  // mPhi = exp(mPhi)/(1.0 + exp(mPhi));
  arma::mat mPi = ones(iT,iM);
  for(m = 0; m < iM; m++){
    for(t = 1; t < iT; t++){
      mPi(t,m) = exp(accu(vZ.t()%cGamma.slice(m).row(t-1)));
    }
    mPi.col(m) = mPi.col(m)/accu(mPi.col(m));
  }
  arma::vec vOmega_Score(iM-1);
  for(m =1; m< iM; m++){
    vOmega_Score(m-1) = vPW_N(m)*(1.0 - vOmega(m));
  }
  // 
  arma::vec vGamma_Score((iT-1)*iP*iM);
  int iFoo2 = 0;
  int iter  = 0;
  arma::vec vGammaScore_foo = zeros((iT-1)*iP,1);
  for(m = 0; m < iM; m++){
    vGammaScore_foo.fill(0.0);
    iter = 0;
    for(t = 1; t < iT; t++){
      vGammaScore_foo.subvec(iter,iter + iP -1) = (mPX(t,m) - mPi(t,m))*(vPW_N(m)*vZ);
      iter += iP;
    }
    vGamma_Score.subvec(iFoo2, iFoo2 + ((iT-1)*iP) -1) = vGammaScore_foo;
    iFoo2 = iFoo2 + ((iT-1)*iP);
  }
  iFoo2 = 0;
  
  arma::ivec ivItemcat_red = ivItemcat -1;
  int nfreepar_res = sum(ivItemcat_red);
  int iItemfoo = 0;
  arma::vec vBeta_Score=zeros(nfreepar_res*iT);
  // first category as reference
  int iroll = 0;
  if(max(ivItemcat)==2){
    for(t = 0; t < iT; t++){
      for(k = 0; k < iK; k++){
        vBeta_Score(iroll) = vPMsumX(t)*(vY(k) - mPhi(k,t));
        iroll += 1;
      }
    }
  }else{
    for(t = 0; t < iT; t++){
      iItemfoo = 0;
      for(v  = 0; v < iV; v++){
        if(vD(v)==1){
          for(l = 1; l < ivItemcat(v); l++){
            if(v > 0){
              iItemfoo = sum(ivItemcat.subvec(0, v-2))+l;
              vBeta_Score(iroll) = vPMsumX(t)*(vY(iItemfoo) - mPhi(iItemfoo,t));
            }else{
              vBeta_Score(iroll) = vPMsumX(t)*(vY(l) - mPhi(l,t));
            }
            iroll += 1;
          }
        }else{
          for(l = 1; l < ivItemcat(v); l++){
            vBeta_Score(iroll) = 0.0;
            iroll += 1;
          }
        }
      }
    }
  }
  arma::mat vScore = join_cols(vOmega_Score,vGamma_Score);
  if(nstep!=3){
    vScore = join_cols(vScore,vBeta_Score);
  }else{
    arma::mat vScore = join_cols(vOmega_Score,vGamma_Score);
  }
  
  return vScore;
}


