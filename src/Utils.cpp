#include <RcppArmadillo.h>
#include "SafeFunctions.h"

using namespace arma;
using namespace Rcpp;



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
    if(K > 2){
      arma::vec ww = zeros(N,1);
      ww.fill(0.0);
      for(j = 0; j < K; j++){
        if(k!=j){
          ww = ww + w.col(j);
        }
      }
      ibeta.slice(k-1) = X.t() * (repmat(ta.col(k),1,P)%X%repmat(ww,1,P)%mWei);
    }
    if(K==2){
      double ww;
      ww = accu(w.col(1));
      for(k=1; k < K; k++){
        ibeta.slice(k-1) = X.t() * (repmat(ta.col(k),1,P)%X%mWei*ww);
      }
    }
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



