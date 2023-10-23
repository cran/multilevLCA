#include <RcppArmadillo.h>
#include "SafeFunctions.h"
#include "Utils.h"

using namespace arma;
using namespace Rcpp;

//[[Rcpp::export]]
List semisup_LCAcov_includeall(arma::mat mY, arma::mat mZ, arma::mat mDesign, arma::mat mPhi, arma::mat mBeta, arma::vec labels, arma::mat mStep1Var,
                               int fixed = 0, int maxIter = 1e3, double tol = 1e-8,
                               double dC = 3/4, double NRtol = 1e-6, int NRmaxit = 100){
  // mY is equal to the number of observed response patterns x iH
  // mDesign has the same dimentions as mY. The its (i,h)-element is == 1 if the
  // corresponding entry in mY is available, ==0 if this is missing (NA)
  // mZ contain a column of ones
  //
  // //   // NO IMPUTATION IS PERFORMED! If one row has only NA, this should be excluded
  // //   // before running this function
  // //   //
  // //   // Missings should be coded with a numeric entry -- like 99 or 999
  // //   //
  // //   // mU used from input if any fixed class label is available. For unlabelled units, rows can be zeroed.
  int iN = mY.n_rows;
  int iH    = mY.n_cols;
  
  int h,k,j,n;
  int isize = 1;
  double size = 1.0;
  int iK = mBeta.n_cols + 1;
  int iT = mZ.n_cols;
  arma::cube cLogdY = zeros(iN,iK,iH);
  arma::mat mLogdKY = zeros(iN,iK);
  arma::vec vLLK = zeros(iN,1);
  arma::mat mU = zeros(iN,iK);
  arma::mat mPg = ones(iN,iK);
  mPg.cols(1,iK-1) = exp(mZ*mBeta);
  for(n = 0; n < iN; n++){
    mPg.row(n) = mPg.row(n)/accu(mPg.row(n));
    if(labels[n]==1.0){
      mU(n,iK-1) = 1.0;
    }
  }
  if(fixed==1){
    for(n = 0; n < iN; n++){
      for(k = 0; k < iK; k++){
        for(h=0; h< iH;h++){
          if(mDesign(n,h)==1){
            cLogdY(n,k,h) = Rf_dbinom(mY(n,h), size, mPhi(h,k), isize);
            mLogdKY(n,k) += cLogdY(n,k,h);
          }
        }
      }
    }
  }
  
  double eps = 1.0;
  double iter = 0.0;
  arma::vec LLKSeries(maxIter);
  arma::mat mbeta = mBeta.t();
  arma::mat mbeta_score(iN,(iK-1)*iT);
  List NR_step;
  
  while(eps > tol && iter<maxIter){
    vLLK.zeros();
    // step E
    
    for(n = 0; n < iN; n++){
      if(fixed==0){
        mLogdKY.row(n).zeros();
        for(k = 0; k < iK; k++){
          for(h=0; h< iH;h++){
            if(mDesign(n,h)==1){
              cLogdY(n,k,h) = Rf_dbinom(mY(n,h), size, mPhi(h,k), isize);
              mLogdKY(n,k) += cLogdY(n,k,h);
            }
          }
        }
      }
      if(labels[n]==1.0){
        vLLK(n) = accu(mU.row(n).t()%log(mPg.row(n).t())) + accu(mU.row(n).t()%mLogdKY.row(n).t());
      }else{
        vLLK(n) = MixtDensityScale(mPg.row(n).t(), mLogdKY.row(n).t(), iK);
        for(k = 0;k < iK-1; k++){
          mU(n,k) = exp(log(mPg(n,k)) + mLogdKY(n,k) - vLLK(n));
        }
        mU.row(n).subvec(0,iK-2) = mU.row(n).subvec(0,iK-2)/accu(mU.row(n).subvec(0,iK-2));
      }
    }
    // Step M
    // Response probabilities
    for(k = 0; k < iK-1; k++){
      for(h = 0; h < iH; h++){
        mPhi(h,k) = accu(mU.col(k)%mY.col(h)%mDesign.col(h))/accu(mU.col(k)%mDesign.col(h));
        mPhi(h,k) = probcheck(mPhi(h,k));
      }
    }
    
    // //     //  Mixture individual specific weights
    NR_step = NR_step_covIT_LS(mZ, mbeta, mU, dC, NRtol, NRmaxit);
    arma::mat w = NR_step["w_i"];
    mPg = w;
    arma::mat mbeta_next = NR_step["beta"];
    mbeta = mbeta_next;
	arma::mat mbeta_score_curr = NR_step["mSbeta"];
    mbeta_score = mbeta_score_curr;
    
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    iter +=1;
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  
  arma::mat gamma(iH,iK);
  gamma = log(mPhi/(1.0 - mPhi));
  
  arma::vec parvec = join_cols(vectorise(mbeta),vectorise(gamma));
  
  double BIC, AIC;
  int npar = (iK-1)*iT + iH*iK;
  
  BIC = -2.0*LLKSeries(iter-1) + log(iN*1.0)*npar;
  AIC = -2.0*LLKSeries(iter-1) + 2.0*npar;
  arma::vec vPg = mean(mU).t();
  double Terr = accu(vPg%log(vPg));
  arma::mat mlogU = trunc_log(mU);
  double Perr = mean(sum(-mU%mlogU,1));
  double R2entr = (Terr - Perr)/Terr;
  
  /// classification error
  arma::ivec vModalAssnm(iN);
  arma::mat mModalAssnm = zeros(iN,iK);
  for(n=0; n< iN; n++){
    vModalAssnm(n) = WhichMax(mU.row(n).t());
    mModalAssnm(n,vModalAssnm(n)) = 1.0;
  }
  
  arma::mat mClassErr     = zeros(iK,iK);
  arma::mat mClassErrProb = mClassErr;
  for(j = 0; j < iK; j++){
    for(k = 0; k < iK; k++){
      mClassErr(j,k) = accu(mU.col(k)%mModalAssnm.col(j));
    }
    mClassErrProb.row(j)=mClassErr.row(j)/accu(mClassErr.row(j));
  }
  
  double dClassErr_tot =  1.0 - (accu(mClassErr.diag())/iN);
  
  
  // 
  // Computing SEs
  // simultaneous estimation
  
  // Computing the score
  
  
  arma::mat mGamma_Score=zeros(iN,iH*iK);
  int iroll = 0;
  for(k = 0; k < iK; k++){
    for(h = 0; h < iH; h++){
      mGamma_Score.col(iroll) = (mU.col(k)%(mY.col(h) - mPhi(h,k)));
      iroll += 1;
    }
  }
  
  iroll=0;
  
  // Expected information matrix 
  
  arma::mat mScore = join_rows(mbeta_score,mGamma_Score);
  arma::mat Infomat = mScore.t()*mScore/iN;
  arma::mat Varmat = pinv(Infomat)/iN;
  arma::vec SEs_unc =  sqrt(Varmat.diag());
  
  // Matrix formula
  int uncondLatpars   = iK - 1;
  int parsfree        = (iK - 1)*iT;
  arma::mat mSigma11  = mStep1Var.submat(uncondLatpars-1,uncondLatpars-1,uncondLatpars + (iK*iH)-1,uncondLatpars + (iK*iH)-1);
  arma::mat mV2       = Varmat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmat     = Infomat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmatInv  = pinv(mJmat); 
  arma::mat mH        = Infomat.submat(0,parsfree-1,parsfree-1,parsfree + (iK*iH)-1);
  arma::mat mQ        =  mJmatInv*mH*mSigma11*mH.t()*mJmatInv;
  arma::mat mVar_corr = mV2 + mQ;
  arma::vec SEs_cor =  SEs_unc;
  if(fixed == 1){
    SEs_cor.subvec(0,parsfree-1) = sqrt(mVar_corr.diag());
  }
  
  List EMout;
  EMout["mU"]            = mU;
  EMout["mPhi"]          = mPhi;
  EMout["mPg"]           = mPg;
  EMout["mBeta"]         = mbeta.t();
  EMout["mGamma"]        = gamma;
  EMout["SEs_unc"]     = SEs_unc;
  EMout["SEs_cor"]     = SEs_cor;
  EMout["mQ"]          = mQ;
  EMout["mV2"]         = mV2;
  EMout["Varmat_unc"]     = Varmat;
  EMout["Varmat_cor"]     = mVar_corr;
  EMout["LLKSeries"]     = LLKSeries;
  EMout["eps"]           = eps;
  EMout["iter"]          = iter;
  EMout["BIC"]           = BIC;
  EMout["AIC"]           = AIC;
  EMout["R2entr"]        = R2entr;
  EMout["mClassErr"]     = mClassErr;
  EMout["mClassErrProb"] = mClassErrProb;
  EMout["dClassErr_tot"] = dClassErr_tot;
  EMout["mModalAssnm"]   = mModalAssnm;
  EMout["vModalAssnm"]   = vModalAssnm;
  EMout["parvec"]        = parvec;
  
  return EMout;
}





//[[Rcpp::export]]
List LCA_fast_includeall(arma::mat mY, arma::mat mDesign, arma::ivec ivFreq, int iK, arma::mat mU, int maxIter = 1e3, double tol = 1e-8, int reord 
= 0){
  // mY is equal to the number of observed response patterns x iH
  // mDesign has the same dimentions as mY. The its (i,h)-element is == 1 if the
  // corresponding entry in mY is available, ==0 if this is missing (NA)
  // 
  // NO IMPUTATION IS PERFORMED! If one row has only NA, this should be excluded
  // before running this function
  // 
  // Missings should be coded with a numeric entry -- like 99 or 999
  int iNtot    = accu(ivFreq);
  int iN = mY.n_rows;
  int iH    = mY.n_cols;
  
  int h,k,j,n;
  int isize = 1;
  double size = 1.0;
  arma::cube mdY = zeros(iN,iK,iH);
  arma::vec vLLK = zeros(iN,1);
  arma::mat mPhi(iH,iK);
  arma::vec pg(iK);
  arma::vec vHdY(iK);
  double eps = 1.0;
  double iter = 0.0;
  arma::vec LLKSeries(maxIter);
  while(eps > tol && iter<maxIter){
    vLLK.zeros();
    // step M
    for(k = 0; k < iK; k++){
      pg(k) = accu(mU.col(k)%ivFreq)/iNtot;
      for(h = 0; h < iH; h++){
        mPhi(h,k) = accu(mU.col(k)%ivFreq%mY.col(h)%mDesign.col(h))/accu(mU.col(k)%ivFreq%mDesign.col(h));
        mPhi(h,k) = probcheck(mPhi(h,k));
      }
    }
    
    pg = OmegaCheck(pg, iK);
    // step E
    
    for(n = 0; n < iN; n++){
      vHdY.zeros();
      for(k = 0; k < iK; k++){
        for(h=0; h< iH;h++){
          if(mDesign(n,h)==1){
            mdY(n,k,h) = Rf_dbinom(mY(n,h), size, mPhi(h,k), isize);
            vHdY(k) += mdY(n,k,h);
          }
        }
      }
      vLLK(n) = MixtDensityScale(pg, vHdY, iK);
      for(k = 0;k < iK; k++){
        mU(n,k) = exp(log(pg(k)) + vHdY(k) - vLLK(n));
      }
    }
    
    LLKSeries(iter) = accu(vLLK%ivFreq);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    iter +=1;
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  if(reord == 1){
    arma::vec vPhisum = sum(mPhi).t();
    arma::uvec order = sort_index(vPhisum,"descending");
    int ifoo = 0;
    arma::mat mPhi_sorted = mPhi;
    arma::vec pg_sorted = pg;
    arma::mat mU_sorted = mU;
    for(k=0; k< iK; k++){
      ifoo               = order(k);
      mPhi_sorted.col(k) = mPhi.col(ifoo);
      pg_sorted(k)       = pg(ifoo);
      mU_sorted.col(k)   = mU.col(ifoo); 
    }
    mPhi = mPhi_sorted;
    pg = pg_sorted;
    mU = mU_sorted;
  }
  arma::vec alphafoo(iK);
  alphafoo = log(pg/pg(0));
  arma::vec alpha(iK-1);
  alpha = alphafoo.subvec(1,iK-1);
  arma::mat gamma(iH,iK);
  gamma = log(mPhi/(1.0 - mPhi));
  
  arma::vec parvec = join_cols(alpha,vectorise(gamma));
  
  double BIC, AIC;
  BIC = -2.0*LLKSeries(iter-1) + log(iNtot*1.0)*1.0*(iH*iK + (iK - 1));
  AIC = -2.0*LLKSeries(iter-1) + 2.0*(iH*iK + (iK - 1));
  double Terr = accu(-pg%log(pg));
  arma::mat mlogU = trunc_log(mU);
  double Perr = sum(sum(-mU%mlogU,1)%ivFreq)/iNtot;
  double R2entr = (Terr - Perr)/Terr;
  
  // classification error
  arma::ivec vModalAssnm(iN);
  arma::mat mModalAssnm = zeros(iN,iK);
  for(n=0; n< iN; n++){
    vModalAssnm(n) = WhichMax(mU.row(n).t());
    mModalAssnm(n,vModalAssnm(n)) = 1.0;
  }
  
  arma::mat mClassErr     = zeros(iK,iK);
  arma::mat mClassErrProb = mClassErr;
  for(j = 0; j < iK; j++){
    for(k = 0; k < iK; k++){
      mClassErr(j,k) = accu(mU.col(k)%mModalAssnm.col(j)%ivFreq);
    }
    mClassErrProb.row(j)=mClassErr.row(j)/accu(mClassErr.row(j));
  }
  
  double dClassErr_tot =  1.0 - (accu(mClassErr.diag())/iNtot);
  
  
  // Computing the score
  
  arma::mat mPg_Score(iN,iK-1);
  for(k =1; k< iK; k++){
    mPg_Score.col(k-1) = (mU.col(k) - pg(k))%ivFreq;
  }
  
  arma::mat mGamma_Score=zeros(iN,iH*iK);
  int iroll = 0;
  for(k = 0; k < iK; k++){
    for(h = 0; h < iH; h++){
      mGamma_Score.col(iroll) = (mU.col(k)%(mY.col(h) - mPhi(h,k)))%ivFreq%mDesign.col(h);
      iroll += 1;
    }
  }
  
  iroll=0;
  
  arma::mat mScore = join_rows(mPg_Score,mGamma_Score);
  arma::mat Infomat = mScore.t()*mScore/iN;
  arma::mat Varmat = pinv(Infomat)/iN;
  arma::vec SEs =  sqrt(Varmat.diag());
  
  
  
  List EMout;
  EMout["mU"]            = mU;
  EMout["mPhi"]          = mPhi;
  EMout["pg"]            = pg;
  EMout["alphas"]        = alpha;
  EMout["gamma"]         = gamma;
  EMout["LLKSeries"]     = LLKSeries;
  EMout["eps"]           = eps;
  EMout["iter"]          = iter;
  EMout["BIC"]           = BIC;
  EMout["AIC"]           = AIC;
  EMout["R2entr"]        = R2entr;
  EMout["mClassErr"]     = mClassErr;
  EMout["mClassErrProb"] = mClassErrProb;
  EMout["dClassErr_tot"] = dClassErr_tot;
  EMout["mModalAssnm"]   = mModalAssnm;
  EMout["vModalAssnm"]   = vModalAssnm;
  EMout["freq"]          = ivFreq;
  EMout["parvec"]    = parvec;
  EMout["Infomat"]   = Infomat;
  EMout["Varmat"]   = Varmat;
  EMout["SEs"]   = SEs;
  
  return EMout;
}





// [[Rcpp::export]]
List LCAcov_poly(arma::mat mY, arma::mat mZ, int iK, arma::mat mPhi, arma::mat mBeta, arma::mat mStep1Var, arma::ivec ivItemcat, int fixed = 0, int maxIter = 1e3, double tol = 1e-8, double NRtol = 1e-6, int NRmaxit = 100){
  // mU must be n x iK
  // ivItemcat is the vector of number of categories for each item
  // fixed = 1 if measurement model parameters should not be updated
  // mStep1Var to be used only if measurement model parameters are kept fixed at the input value
  // sample = 0 if standard matrix formula for Var(Beta) is to be computed; sample = 1 for the Monte Carlo approach
  // mBeta is the iP+1 x iK-1 matrix of structural regression coefficients
  int iN    = mY.n_rows;
  int iH    = mY.n_cols;
  int NT    = iN;
  int iP    = mBeta.n_rows;
  int iV    = ivItemcat.n_elem;
  int h,j,k,n,p,v;
  int isize = 1;
  double size = 1.0;
  arma::cube mdY = zeros(NT,iK,iV);
  arma::mat mLogDensY = zeros(NT,iK);
  arma::vec vLLK = zeros(NT,1);
  arma::mat mPg = ones(NT,iK);
  arma::mat mU = zeros(NT,iK);
  for(k = 0; k < (iK - 1); k++){
    mPg.col(k + 1) = exp(mZ*mBeta.col(k));
  }
  arma::vec sumPg = sum(mPg,1);
  for(n = 0; n < iN; n++){
    mPg.row(n) = mPg.row(n)/sumPg(n);
  }
  arma::vec vHdY(iK);
  double eps = 1.0;
  double iter = 0.0;
  int ifooDcat = 0;
  // Compute log densities
  if(fixed > 0){
    for(n = 0; n < iN; n++){
      for(k = 0; k < iK; k++){
        ifooDcat = 0;
        for(v = 0; v< iV;v++){
          if(ivItemcat(v)==2){
            mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
            ifooDcat += 1;
          }else{
            for(p = 0; p < ivItemcat(v); p++){
              if(mY(n,ifooDcat) > 0.0){
                mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
              }
              ifooDcat += 1;
            }
          }
          mLogDensY(n,k) += mdY(n,k,v);
        }
      }
    }
  }
  
  // 
  arma::mat mbeta = mBeta.t();
  arma::mat mbeta_score(iN,(iK-1)*iP);
  arma::vec LLKSeries(maxIter);
  List NR_step;
  while(eps > tol && iter<maxIter){
    // 
    // compute log densities
    if(fixed == 0){
      mLogDensY.zeros();
      for(n = 0; n < iN; n++){
        for(k = 0; k < iK; k++){
          ifooDcat = 0;
          for(v = 0; v< iV;v++){
            if(ivItemcat(v)==2){
              mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
              ifooDcat += 1;
            }else{
              for(p = 0; p < ivItemcat(v); p++){
                if(mY(n,ifooDcat) > 0.0){
                  mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
                }
                ifooDcat += 1;
              }
            }
            mLogDensY(n,k) += mdY(n,k,v);
          }
        }
      }
    }
    // E step
    for(n = 0; n < iN; n++){
      vLLK(n) = MixtDensityScale(mPg.row(n).t(), mLogDensY.row(n).t(), iK);
      for(k = 0; k < iK; k++){
        mU(n,k) = exp(log(mPg(n,k)) + mLogDensY(n,k) - vLLK(n));
      } 
    }
    // step M
    if(fixed == 0){
      for(k = 0; k < iK; k++){
        for(h = 0; h < iH; h++){
          mPhi(h,k) = accu(mU.col(k)%mY.col(h))/accu(mU.col(k));
          mPhi(h,k) = probcheck(mPhi(h,k));
        }
      }
    }
    NR_step = NR_step_covIT(mZ, mbeta, mU, NRtol, NRmaxit);
    arma::mat w_i = NR_step["w_i"];
    mPg = w_i;
    arma::mat mbeta_next = NR_step["beta"];
    arma::mat mbeta_score_curr = NR_step["mSbeta"];
    mbeta_score = mbeta_score_curr;
    mbeta = mbeta_next;
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    iter +=1;
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  arma::mat gamma(iH,iK);
  gamma = log(mPhi/(1.0 - mPhi));
  // 
  // mbeta is iK-1 x iP+1
  // 
  arma::vec parvec = join_cols(vectorise(mbeta.t()),vectorise(gamma));
  
  double BIC,AIC;
  BIC = -2.0*LLKSeries(iter-1) + log(NT)*(iH*iK + (iK - 1));
  AIC = -2.0*LLKSeries(iter-1) + 2*(iH*iK + (iK - 1));
  // double Terr = accu(-mPg%log(mPg));
  arma::vec vPg = mean(mU).t();
  double Terr = accu(-vPg%log(vPg));
  arma::mat mlogU = trunc_log(mU);
  double Perr = mean(sum(-mU%mlogU,1));
  double R2entr = (Terr - Perr)/Terr;
  
  /// classification error
  arma::ivec vModalAssnm(iN);
  arma::mat mModalAssnm = zeros(iN,iK);
  for(n=0; n< iN; n++){
    vModalAssnm(n) = WhichMax(mU.row(n).t());
    mModalAssnm(n,vModalAssnm(n)) = 1.0;
  }
  
  arma::mat mClassErr     = zeros(iK,iK);
  arma::mat mClassErrProb = mClassErr;
  for(j = 0; j < iK; j++){
    for(k = 0; k < iK; k++){
      mClassErr(j,k) = accu(mU.col(k)%mModalAssnm.col(j));
    }
    mClassErrProb.row(j)=mClassErr.row(j)/accu(mClassErr.row(j));
  }
  
  double dClassErr_tot =  1.0 - (accu(mClassErr.diag())/iN);
  
  // 
  // Computing SEs
  // simultaneous estimation
  
  // Computing the score
  
  
  arma::mat mGamma_Score=zeros(iN,iH*iK);
  int iroll = 0;
  for(k = 0; k < iK; k++){
    for(h = 0; h < iH; h++){
      mGamma_Score.col(iroll) = (mU.col(k)%(mY.col(h) - mPhi(h,k)));
      iroll += 1;
    }
  }
  
  iroll=0;
  
  // Expected information matrix 
  
  arma::mat mScore = join_rows(mbeta_score,mGamma_Score);
  arma::mat Infomat = mScore.t()*mScore/iN;
  arma::mat Varmat = pinv(Infomat)/iN;
  arma::vec SEs_unc =  sqrt(Varmat.diag());
  
  // Matrix formula
  int uncondLatpars   = iK - 1;
  int parsfree        = (iK - 1)*iP;
  arma::mat mSigma11  = mStep1Var.submat(uncondLatpars-1,uncondLatpars-1,uncondLatpars + (iK*iH)-1,uncondLatpars + (iK*iH)-1);
  arma::mat mV2       = Varmat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmat     = Infomat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmatInv  = pinv(mJmat); 
  arma::mat mH        = Infomat.submat(0,parsfree-1,parsfree-1,parsfree + (iK*iH)-1);
  arma::mat mQ        =  mJmatInv*mH*mSigma11*mH.t()*mJmatInv;
  arma::mat mVar_corr = mV2 + mQ;
  arma::vec SEs_cor =  SEs_unc;
  if(fixed == 1){
    SEs_cor.subvec(0,parsfree-1) = sqrt(mVar_corr.diag());
  }
  
  // if(sample==1){
  //   // this subroutine samples from the approximate sampling distribution of the first step estimator
  //   
  // }
  // 
  List EMout;
  EMout["mU"]        = mU;
  EMout["mPhi"]      = mPhi;
  EMout["mPg"]        = mPg;
  EMout["beta"]    = mbeta.t();
  EMout["gamma"]     = gamma;
  EMout["SEs_unc"]     = SEs_unc;
  EMout["SEs_cor"]     = SEs_cor;
  EMout["mQ"]          = mQ;
  EMout["mV2"]         = mV2;
  EMout["Varmat_unc"]     = Varmat;
  EMout["Varmat_cor"]     = mVar_corr;
  EMout["LLKSeries"] = LLKSeries;
  EMout["eps"]       = eps;
  EMout["iter"]      = iter;
  EMout["BIC"]       = BIC;
  EMout["AIC"]       = AIC;
  EMout["R2entr"]    = R2entr;
  EMout["mClassErr"]     = mClassErr;
  EMout["mClassErrProb"] = mClassErrProb;
  EMout["dClassErr_tot"] = dClassErr_tot;
  EMout["vModalAssnm"]   = vModalAssnm;
  EMout["mModalAssnm"]   = mModalAssnm;
  EMout["parvec"]    = parvec;
  
  return EMout;
}  

//[[Rcpp::export]]
List LCA_poly(arma::mat mY, int iK, arma::mat mU, arma::ivec ivItemcat, int maxIter = 1e3, double tol = 1e-8, int reord = 0){
  // mU must be n*Ti x iK
  // ivItemcat is the vector of number of categories for each item
  //
  int iN    = mY.n_rows;
  int iH    = mY.n_cols;
  int iV    = ivItemcat.n_elem;
  int NT    = iN;
  int h,j,k,n,p,v;
  int isize = 1;
  double size = 1.0;
  arma::cube mdY = zeros(NT,iK,iV);
  arma::vec vLLK = zeros(NT,1);
  arma::mat mPhi(iH,iK);
  arma::vec pg(iK);
  arma::vec vHdY(iK);
  double eps = 1.0;
  double iter = 0.0;
  int ifooDcat = 0;
  arma::vec LLKSeries(maxIter);
  while(eps > tol && iter<maxIter){
    vLLK.zeros();
    // step M
    for(k = 0; k < iK; k++){
      pg(k) = mean(mU.col(k));
      for(h = 0; h < iH; h++){
        mPhi(h,k) = accu(mU.col(k)%mY.col(h))/accu(mU.col(k));
        mPhi(h,k) = probcheck(mPhi(h,k));
      }
    }
    
    pg = OmegaCheck(pg, iK);
    // step E
    
    for(n = 0; n < iN; n++){
      vHdY.zeros();
      for(k = 0; k < iK; k++){
        ifooDcat = 0;
        for(v = 0; v< iV;v++){
          if(ivItemcat(v)==2){
            mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
            ifooDcat += 1;
          }else{
            for(p = 0; p < ivItemcat(v); p++){
              if(mY(n,ifooDcat) > 0.0){
                mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
              }
              ifooDcat += 1;
            }
          }
          vHdY(k) += mdY(n,k,v);
        }
      }
      vLLK(n) = MixtDensityScale(pg, vHdY, iK);
      for(k = 0;k < iK; k++){
        mU(n,k) = exp(log(pg(k)) + vHdY(k) - vLLK(n));
      }
    }
    
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    iter +=1;
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  if(reord == 1){
    arma::vec vPhisum = sum(mPhi).t();
    arma::uvec order = sort_index(vPhisum,"descending");
    int ifoo = 0;
    arma::mat mPhi_sorted = mPhi;
    arma::vec pg_sorted = pg;
    arma::mat mU_sorted = mU;
    for(k=0; k< iK; k++){
      ifoo               = order(k);
      mPhi_sorted.col(k) = mPhi.col(ifoo);
      pg_sorted(k)       = pg(ifoo);
      mU_sorted.col(k)   = mU.col(ifoo); 
    }
    mPhi = mPhi_sorted;
    pg = pg_sorted;
    mU = mU_sorted;
  }
  arma::vec alphafoo(iK);
  alphafoo = log(pg/pg(0));
  arma::vec alpha(iK-1);
  alpha = alphafoo.subvec(1,iK-1);
  arma::mat gamma(iH,iK);
  gamma = log(mPhi/(1.0 - mPhi));
  
  arma::vec parvec = join_cols(alpha,vectorise(gamma));
  
  double BIC, AIC;
  BIC = -2.0*LLKSeries(iter-1) + log(NT)*(iH*iK + (iK - 1));
  AIC = -2.0*LLKSeries(iter-1) + 2*(iH*iK + (iK - 1));
  double Terr = accu(-pg%log(pg));
  arma::mat mlogU = trunc_log(mU);
  double Perr = mean(sum(-mU%mlogU,1));
  double R2entr = (Terr - Perr)/Terr;
  
  // classification error
  arma::ivec vModalAssnm(iN);
  arma::mat mModalAssnm = zeros(iN,iK);
  for(n=0; n< iN; n++){
    vModalAssnm(n) = WhichMax(mU.row(n).t());
    mModalAssnm(n,vModalAssnm(n)) = 1.0;
  }
  
  arma::mat mClassErr     = zeros(iK,iK);
  arma::mat mClassErrProb = mClassErr;
  for(j = 0; j < iK; j++){
    for(k = 0; k < iK; k++){
      mClassErr(j,k) = accu(mU.col(k)%mModalAssnm.col(j));
    }
    mClassErrProb.row(j)=mClassErr.row(j)/accu(mClassErr.row(j));
  }
  
  double dClassErr_tot =  1.0 - (accu(mClassErr.diag())/iN);
  
  
  List EMout;
  EMout["mU"]            = mU;
  EMout["mPhi"]          = mPhi;
  EMout["pg"]            = pg;
  EMout["alphas"]        = alpha;
  EMout["gamma"]         = gamma;
  EMout["LLKSeries"]     = LLKSeries;
  EMout["eps"]           = eps;
  EMout["iter"]          = iter;
  EMout["BIC"]           = BIC;
  EMout["AIC"]           = AIC;
  EMout["R2entr"]        = R2entr;
  EMout["mClassErr"]     = mClassErr;
  EMout["mClassErrProb"] = mClassErrProb;
  EMout["dClassErr_tot"] = dClassErr_tot;
  EMout["vModalAssnm"]   = vModalAssnm;
  EMout["mModalAssnm"]   = mModalAssnm;
  EMout["parvec"]    = parvec;
  
  return EMout;
}

//[[Rcpp::export]]
List LCA_fast_poly(arma::mat mY, arma::ivec ivFreq, int iK, arma::mat mU, arma::ivec ivItemcat,  int maxIter = 1e3, double tol = 1e-8, int reord = 0){
  // mY is equal to the number of observed response patterns x iH
  // ivItempos tells the number of categories for each item 
  //
  int iNtot    = accu(ivFreq);
  int iN = mY.n_rows;
  int iH    = mY.n_cols;
  int iV    = ivItemcat.n_elem;
  int NT    = iN;
  int h,j,k,n,p,v;
  int isize = 1;
  double size = 1.0;
  arma::cube mdY = zeros(NT,iK,iV);
  arma::vec vLLK = zeros(NT,1);
  arma::mat mPhi(iH,iK);
  arma::vec pg(iK);
  arma::vec vHdY(iK);
  double eps = 1.0;
  double iter = 0.0;
  int ifooDcat = 0;
  arma::vec LLKSeries(maxIter);
  while(eps > tol && iter<maxIter){
    vLLK.zeros();
    // step M
    for(k = 0; k < iK; k++){
      pg(k) = accu(mU.col(k)%ivFreq)/iNtot;
      for(h = 0; h < iH; h++){
        mPhi(h,k) = accu(mU.col(k)%ivFreq%mY.col(h))/accu(mU.col(k)%ivFreq);
        mPhi(h,k) = probcheck(mPhi(h,k));
      }
    }
    
    pg = OmegaCheck(pg, iK);
    // step E
    for(n = 0; n < iN; n++){
      vHdY.zeros();
      for(k = 0; k < iK; k++){
        ifooDcat = 0;
        for(v = 0; v< iV;v++){
          if(ivItemcat(v)==2){
            mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
            ifooDcat += 1;
          }else{
            for(p = 0; p < ivItemcat(v); p++){
              if(mY(n,ifooDcat) > 0.0){
                mdY(n,k,v) = Rf_dbinom(mY(n,ifooDcat), size, mPhi(ifooDcat,k), isize);
              }
              ifooDcat += 1;
            }
          }
          vHdY(k) += mdY(n,k,v);
        }
      }
      vLLK(n) = MixtDensityScale(pg, vHdY, iK);
      for(k = 0;k < iK; k++){
        mU(n,k) = exp(log(pg(k)) + vHdY(k) - vLLK(n));
      }
    }
    LLKSeries(iter) = accu(vLLK%ivFreq);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    iter +=1;
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  if(reord == 1){
    arma::vec vPhisum = sum(mPhi).t();
    arma::uvec order = sort_index(vPhisum,"descending");
    int ifoo = 0;
    arma::mat mPhi_sorted = mPhi;
    arma::vec pg_sorted = pg;
    arma::mat mU_sorted = mU;
    for(k=0; k< iK; k++){
      ifoo               = order(k);
      mPhi_sorted.col(k) = mPhi.col(ifoo);
      pg_sorted(k)       = pg(ifoo);
      mU_sorted.col(k)   = mU.col(ifoo); 
    }
    mPhi = mPhi_sorted;
    pg = pg_sorted;
    mU = mU_sorted;
  }
  arma::vec alphafoo(iK);
  alphafoo = log(pg/pg(0));
  arma::vec alpha(iK-1);
  alpha = alphafoo.subvec(1,iK-1);
  arma::mat gamma(iH,iK);
  gamma = log(mPhi/(1.0 - mPhi));
  
  arma::vec parvec = join_cols(alpha,vectorise(gamma));
  
  double BIC, AIC;
  BIC = -2.0*LLKSeries(iter-1) + log(iNtot)*(iH*iK + (iK - 1));
  AIC = -2.0*LLKSeries(iter-1) + 2*(iH*iK + (iK - 1));
  double Terr = accu(-pg%log(pg));
  arma::mat mlogU = trunc_log(mU);
  double Perr = sum(sum(-mU%mlogU,1)%ivFreq)/iNtot;
  double R2entr = (Terr - Perr)/Terr;
  
  // classification error
  arma::ivec vModalAssnm(iN);
  arma::mat mModalAssnm = zeros(iN,iK);
  for(n=0; n< iN; n++){
    vModalAssnm(n) = WhichMax(mU.row(n).t());
    mModalAssnm(n,vModalAssnm(n)) = 1.0;
  }
  
  arma::mat mClassErr     = zeros(iK,iK);
  arma::mat mClassErrProb = mClassErr;
  for(j = 0; j < iK; j++){
    for(k = 0; k < iK; k++){
      mClassErr(j,k) = accu(mU.col(k)%mModalAssnm.col(j)%ivFreq);
    }
    mClassErrProb.row(j)=mClassErr.row(j)/accu(mClassErr.row(j));
  }
  
  double dClassErr_tot =  1.0 - (accu(mClassErr.diag())/iNtot);
  
  
  // Computing the score
  
  arma::mat mPg_Score(iN,iK-1);
  for(k =1; k< iK; k++){
    mPg_Score.col(k-1) = (mU.col(k) - pg(k))%ivFreq;
  }
  
  arma::mat mGamma_Score=zeros(iN,iH*iK);
  int iroll = 0;
  for(k = 0; k < iK; k++){
    for(h = 0; h < iH; h++){
      mGamma_Score.col(iroll) = (mU.col(k)%(mY.col(h) - mPhi(h,k)))%ivFreq;
      iroll += 1;
    }
  }
  
  iroll=0;
  
  arma::mat mScore = join_rows(mPg_Score,mGamma_Score);
  arma::mat Infomat = mScore.t()*mScore/iN;
  arma::mat Varmat = pinv(Infomat)/iN;
  arma::vec SEs =  sqrt(Varmat.diag());
  
  List EMout;
  EMout["mU"]            = mU;
  EMout["mPhi"]          = mPhi;
  EMout["pg"]            = pg;
  EMout["alphas"]        = alpha;
  EMout["gamma"]         = gamma;
  EMout["mScore"] =mScore;
  EMout["Varmat"] =Varmat;
  EMout["SEs"] =SEs;
  EMout["LLKSeries"]     = LLKSeries;
  EMout["eps"]           = eps;
  EMout["iter"]          = iter;
  EMout["BIC"]           = BIC;
  EMout["AIC"]           = AIC;
  EMout["R2entr"]        = R2entr;
  EMout["mClassErr"]     = mClassErr;
  EMout["mClassErrProb"] = mClassErrProb;
  EMout["dClassErr_tot"] = dClassErr_tot;
  EMout["mModalAssnm"]   = mModalAssnm;
  EMout["vModalAssnm"]   = vModalAssnm;
  EMout["freq"]          = ivFreq;
  EMout["parvec"]    = parvec;
  
  return EMout;
}

//[[Rcpp::export]]
double LCA_LLK(arma::vec parvec, arma::mat mY, int iK){
  // mU must be n*Ti x iK
  // works only with dichotomous items
  int iN    = mY.n_rows;
  int iH    = mY.n_cols;
  int NT    = iN;
  int h;
  int k;
  int n;
  int isize = 1;
  double size = 1.0;
  arma::cube mdY = zeros(NT,iK,iH);
  arma::vec vLLK = zeros(NT,1);
  arma::vec pg=ones(iK,1);
  pg.subvec(1,iK-1) = exp(parvec.subvec(0,iK-2));
  pg = pg/accu(pg);
  
  int foopar = iK-1;
  arma::mat mPhi = reshape(conv_to<mat>::from(parvec(span(foopar,foopar + iH*iK -1))),iH,iK);
  mPhi = exp(mPhi)/(1.0 + exp(mPhi));
  
  arma::vec vHdY(iK);
  
  for(n = 0; n < iN; n++){
    vHdY.zeros();
    for(k = 0; k < iK; k++){
      for(h=0; h< iH;h++){
        mdY(n,k,h) = Rf_dbinom(mY(n,h), size, mPhi(h,k), isize);
        vHdY(k) += mdY(n,k,h);
      }
    }
    vLLK(n) = MixtDensityScale(pg, vHdY, iK);
  }
  
  double log_like = accu(vLLK);
  
  return log_like;
}

//[[Rcpp::export]]
arma::vec LCA_LLK_j(arma::vec parvec, arma::mat mY, int iK){
  // mU must be n*Ti x iK
  // works only with dichotomous items
  int iN    = mY.n_rows;
  int iH    = mY.n_cols;
  int NT    = iN;
  int h;
  int k;
  int n;
  int isize = 1;
  double size = 1.0;
  arma::cube mdY = zeros(NT,iK,iH);
  arma::vec vLLK = zeros(NT,1);
  arma::vec pg=ones(iK,1);
  pg.subvec(1,iK-1) = exp(parvec.subvec(0,iK-2));
  pg = pg/accu(pg);
  
  int foopar = iK-1;
  arma::mat mPhi = reshape(conv_to<mat>::from(parvec(span(foopar,foopar + iH*iK -1))),iH,iK);
  mPhi = exp(mPhi)/(1.0 + exp(mPhi));
  
  arma::vec vHdY(iK);
  
  for(n = 0; n < iN; n++){
    vHdY.zeros();
    for(k = 0; k < iK; k++){
      for(h=0; h< iH;h++){
        mdY(n,k,h) = Rf_dbinom(mY(n,h), size, mPhi(h,k), isize);
        vHdY(k) += mdY(n,k,h);
      }
    }
    vLLK(n) = MixtDensityScale(pg, vHdY, iK);
  }
  
  return vLLK;
}

// [[Rcpp::export]]
List LCAcov(arma::mat mY, arma::mat mZ, int iK, arma::mat mPhi, arma::mat mBeta, arma::mat mStep1Var, int fixed = 0, int maxIter = 1e3, double tol = 1e-8, double NRtol = 1e-6, int NRmaxit = 100){
  // mU must be n x iK
  //
  // fixed = 1 if measurement model parameters should not be updated
  // mStep1Var to be used only if measurement model parameters are kept fixed at the input value
  // sample = 0 if standard matrix formula for Var(Beta) is to be computed; sample = 1 for the Monte Carlo approach
  // mBeta is the iP+1 x iK-1 matrix of structural regression coefficients
  int iN    = mY.n_rows;
  int iH    = mY.n_cols;
  int NT    = iN;
  int iP    = mBeta.n_rows;
  int h,j,k,n;
  int isize = 1;
  double size = 1.0;
  arma::cube mdY = zeros(NT,iK,iH);
  arma::vec vLLK = zeros(NT,1);
  arma::mat mPg = ones(NT,iK);
  arma::mat mU = zeros(NT,iK);
  for(k = 0; k < (iK - 1); k++){
    mPg.col(k + 1) = exp(mZ*mBeta.col(k));
  }
  arma::vec sumPg = sum(mPg,1);
  for(n = 0; n < iN; n++){
    mPg.row(n) = mPg.row(n)/sumPg(n);
  }
  arma::vec vHdY(iK);
  double eps = 1.0;
  double iter = 0.0;
  // arma::mat mbeta = zeros(iK-1,iP);
  arma::mat mbeta = mBeta.t();
  arma::mat mbeta_score(iN,(iK-1)*iP);
  arma::vec LLKSeries(maxIter);
  List NR_step;
  while(eps > tol && iter<maxIter){
    vLLK.zeros();
    
    // step E
    
    for(n = 0; n < iN; n++){
      vHdY.zeros();
      for(k = 0; k < iK; k++){
        for(h=0; h< iH;h++){
          mdY(n,k,h) = Rf_dbinom(mY(n,h), size, mPhi(h,k), isize);
          vHdY(k) += mdY(n,k,h);
        }
      }
      vLLK(n) = MixtDensityScale(mPg.row(n).t(), vHdY, iK);
      for(k = 0;k < iK; k++){
        mU(n,k) = exp(log(mPg(n,k)) + vHdY(k) - vLLK(n));
      }
    }
    
    // step M
    if(fixed==0){
      for(k = 0; k < iK; k++){
        // pg(k) = mean(mU.col(k));
        for(h = 0; h < iH; h++){
          mPhi(h,k) = accu(mU.col(k)%mY.col(h))/accu(mU.col(k));
          mPhi(h,k) = probcheck(mPhi(h,k));
        }
      }
    }
    if(fixed==1){
      for(k = 0; k < iK; k++){
        // pg(k) = mean(mU.col(k));
        for(h = 0; h < iH; h++){
          mPhi(h,k) = probcheck(mPhi(h,k));
        }
      }
    }
    NR_step = NR_step_covIT(mZ, mbeta, mU, NRtol, NRmaxit);
    arma::mat w_i = NR_step["w_i"];
    mPg = w_i;
    arma::mat mbeta_next = NR_step["beta"];
    arma::mat mbeta_score_curr = NR_step["mSbeta"];
    mbeta_score = mbeta_score_curr;
    mbeta = mbeta_next;
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    iter +=1;
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  arma::mat gamma(iH,iK);
  gamma = log(mPhi/(1.0 - mPhi));
  // 
  // mbeta is iK-1 x iP+1
  // 
  arma::vec parvec = join_cols(vectorise(mbeta.t()),vectorise(gamma));
  
  double BIC,AIC;
  BIC = -2.0*LLKSeries(iter-1) + log(NT)*(iH*iK + (iK - 1));
  AIC = -2.0*LLKSeries(iter-1) + 2*(iH*iK + (iK - 1));
  // double Terr = accu(-mPg%log(mPg));
  arma::vec vPg = mean(mU).t();
  double Terr = accu(-vPg%log(vPg));
  arma::mat mlogU = trunc_log(mU);
  double Perr = mean(sum(-mU%mlogU,1));
  double R2entr = (Terr - Perr)/Terr;
  
  /// classification error
  arma::ivec vModalAssnm(iN);
  arma::mat mModalAssnm = zeros(iN,iK);
  for(n=0; n< iN; n++){
    vModalAssnm(n) = WhichMax(mU.row(n).t());
    mModalAssnm(n,vModalAssnm(n)) = 1.0;
  }
  
  arma::mat mClassErr     = zeros(iK,iK);
  arma::mat mClassErrProb = mClassErr;
  for(j = 0; j < iK; j++){
    for(k = 0; k < iK; k++){
      mClassErr(j,k) = accu(mU.col(k)%mModalAssnm.col(j));
    }
    mClassErrProb.row(j)=mClassErr.row(j)/accu(mClassErr.row(j));
  }
  
  double dClassErr_tot =  1.0 - (accu(mClassErr.diag())/iN);
  
  // 
  // Computing SEs
  // simultaneous estimation
  
  // Computing the score
  
  
  arma::mat mGamma_Score=zeros(iN,iH*iK);
  int iroll = 0;
  for(k = 0; k < iK; k++){
    for(h = 0; h < iH; h++){
      mGamma_Score.col(iroll) = (mU.col(k)%(mY.col(h) - mPhi(h,k)));
      iroll += 1;
    }
  }
  
  iroll=0;
  
  // Expected information matrix 
  
  arma::mat mScore = join_rows(mbeta_score,mGamma_Score);
  arma::mat Infomat = mScore.t()*mScore/iN;
  arma::mat Varmat = pinv(Infomat)/iN;
  arma::vec SEs_unc =  sqrt(Varmat.diag());
  
  // Matrix formula
  int uncondLatpars   = iK - 1;
  int parsfree        = (iK - 1)*iP;
  arma::mat mSigma11  = mStep1Var.submat(uncondLatpars-1,uncondLatpars-1,uncondLatpars + (iK*iH)-1,uncondLatpars + (iK*iH)-1);
  arma::mat mV2       = Varmat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmat     = Infomat.submat(0,0,parsfree-1,parsfree-1);
  arma::mat mJmatInv  = pinv(mJmat); 
  arma::mat mH        = Infomat.submat(0,parsfree-1,parsfree-1,parsfree + (iK*iH)-1);
  arma::mat mQ        =  mJmatInv*mH*mSigma11*mH.t()*mJmatInv;
  arma::mat mVar_corr = mV2 + mQ;
  arma::vec SEs_cor =  SEs_unc;
  if(fixed == 1){
    SEs_cor.subvec(0,parsfree-1) = sqrt(mVar_corr.diag());
  }
  
  // if(sample==1){
  //   // this subroutine samples from the approximate sampling distribution of the first step estimator
  //   
  // }
  // 
  List EMout;
  EMout["mU"]        = mU;
  EMout["mPhi"]      = mPhi;
  EMout["mPg"]        = mPg;
  EMout["beta"]    = mbeta.t();
  EMout["gamma"]     = gamma;
  EMout["SEs_unc"]     = SEs_unc;
  EMout["SEs_cor"]     = SEs_cor;
  EMout["mQ"]          = mQ;
  EMout["mV2"]         = mV2;
  EMout["Varmat_unc"]     = Varmat;
  EMout["Varmat_cor"]     = mVar_corr;
  EMout["LLKSeries"] = LLKSeries;
  EMout["eps"]       = eps;
  EMout["iter"]      = iter;
  EMout["BIC"]       = BIC;
  EMout["AIC"]       = AIC;
  EMout["R2entr"]    = R2entr;
  EMout["mClassErr"]     = mClassErr;
  EMout["mClassErrProb"] = mClassErrProb;
  EMout["dClassErr_tot"] = dClassErr_tot;
  EMout["vModalAssnm"]   = vModalAssnm;
  EMout["mModalAssnm"]   = mModalAssnm;
  EMout["parvec"]    = parvec;
  
  return EMout;
}

//[[Rcpp::export]]
List LCA(arma::mat mY, int iK, arma::mat mU, int maxIter = 1e3, double tol = 1e-8, int reord = 0){
  // mU must be n*Ti x iK
  //
  int iN    = mY.n_rows;
  int iH    = mY.n_cols;
  int NT    = iN;
  int h,k,j,n;
  int isize = 1;
  double size = 1.0;
  arma::cube mdY = zeros(NT,iK,iH);
  arma::vec vLLK = zeros(NT,1);
  arma::mat mPhi(iH,iK);
  arma::vec pg(iK);
  arma::vec vHdY(iK);
  double eps = 1.0;
  double iter = 0.0;
  arma::vec LLKSeries(maxIter);
  while(eps > tol && iter<maxIter){
    vLLK.zeros();
    // step M
    for(k = 0; k < iK; k++){
      pg(k) = mean(mU.col(k));
      for(h = 0; h < iH; h++){
        mPhi(h,k) = accu(mU.col(k)%mY.col(h))/accu(mU.col(k));
        mPhi(h,k) = probcheck(mPhi(h,k));
      }
    }
    
    pg = OmegaCheck(pg, iK);
    // step E
    
    for(n = 0; n < iN; n++){
      vHdY.zeros();
      for(k = 0; k < iK; k++){
        for(h=0; h< iH;h++){
          mdY(n,k,h) = Rf_dbinom(mY(n,h), size, mPhi(h,k), isize);
          vHdY(k) += mdY(n,k,h);
        }
      }
      vLLK(n) = MixtDensityScale(pg, vHdY, iK);
      for(k = 0;k < iK; k++){
        mU(n,k) = exp(log(pg(k)) + vHdY(k) - vLLK(n));
      }
    }
    
    LLKSeries(iter) = accu(vLLK);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    iter +=1;
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  if(reord == 1){
    arma::vec vPhisum = sum(mPhi).t();
    arma::uvec order = sort_index(vPhisum,"descending");
    int ifoo = 0;
    arma::mat mPhi_sorted = mPhi;
    arma::vec pg_sorted = pg;
    arma::mat mU_sorted = mU;
    for(k=0; k< iK; k++){
      ifoo               = order(k);
      mPhi_sorted.col(k) = mPhi.col(ifoo);
      pg_sorted(k)       = pg(ifoo);
      mU_sorted.col(k)   = mU.col(ifoo); 
    }
    mPhi = mPhi_sorted;
    pg = pg_sorted;
    mU = mU_sorted;
  }
  arma::vec alphafoo(iK);
  alphafoo = log(pg/pg(0));
  arma::vec alpha(iK-1);
  alpha = alphafoo.subvec(1,iK-1);
  arma::mat gamma(iH,iK);
  gamma = log(mPhi/(1.0 - mPhi));
  
  arma::vec parvec = join_cols(alpha,vectorise(gamma));
  
  double BIC, AIC;
  BIC = -2.0*LLKSeries(iter-1) + log(NT)*(iH*iK + (iK - 1));
  AIC = -2.0*LLKSeries(iter-1) + 2*(iH*iK + (iK - 1));
  double Terr = accu(-pg%log(pg));
  arma::mat mlogU = trunc_log(mU);
  double Perr = mean(sum(-mU%mlogU,1));
  double R2entr = (Terr - Perr)/Terr;
  
  // classification error
  arma::ivec vModalAssnm(iN);
  arma::mat mModalAssnm = zeros(iN,iK);
  for(n=0; n< iN; n++){
    vModalAssnm(n) = WhichMax(mU.row(n).t());
    mModalAssnm(n,vModalAssnm(n)) = 1.0;
  }
  
  arma::mat mClassErr     = zeros(iK,iK);
  arma::mat mClassErrProb = mClassErr;
  for(j = 0; j < iK; j++){
    for(k = 0; k < iK; k++){
      mClassErr(j,k) = accu(mU.col(k)%mModalAssnm.col(j));
    }
    mClassErrProb.row(j)=mClassErr.row(j)/accu(mClassErr.row(j));
  }
  
  double dClassErr_tot =  1.0 - (accu(mClassErr.diag())/iN);
  
  
  List EMout;
  EMout["mU"]            = mU;
  EMout["mPhi"]          = mPhi;
  EMout["pg"]            = pg;
  EMout["alphas"]        = alpha;
  EMout["gamma"]         = gamma;
  EMout["LLKSeries"]     = LLKSeries;
  EMout["eps"]           = eps;
  EMout["iter"]          = iter;
  EMout["BIC"]           = BIC;
  EMout["AIC"]           = AIC;
  EMout["R2entr"]        = R2entr;
  EMout["mClassErr"]     = mClassErr;
  EMout["mClassErrProb"] = mClassErrProb;
  EMout["dClassErr_tot"] = dClassErr_tot;
  EMout["vModalAssnm"]   = vModalAssnm;
  EMout["mModalAssnm"]   = mModalAssnm;
  EMout["parvec"]    = parvec;
  
  return EMout;
}

//[[Rcpp::export]]
List LCA_fast(arma::mat mY, arma::ivec ivFreq, int iK, arma::mat mU, int maxIter = 1e3, double tol = 1e-8, int reord = 0){
  // mY is equal to the number of observed response patterns x iH
  //
  int iNtot    = accu(ivFreq);
  int iN = mY.n_rows;
  int iH    = mY.n_cols;
  int NT    = iN;
  int h,k,j,n;
  int isize = 1;
  double size = 1.0;
  arma::cube mdY = zeros(NT,iK,iH);
  arma::vec vLLK = zeros(NT,1);
  arma::mat mPhi(iH,iK);
  arma::vec pg(iK);
  arma::vec vHdY(iK);
  double eps = 1.0;
  double iter = 0.0;
  arma::vec LLKSeries(maxIter);
  while(eps > tol && iter<maxIter){
    vLLK.zeros();
    // step M
    for(k = 0; k < iK; k++){
      pg(k) = accu(mU.col(k)%ivFreq)/iNtot;
      for(h = 0; h < iH; h++){
        mPhi(h,k) = accu(mU.col(k)%ivFreq%mY.col(h))/accu(mU.col(k)%ivFreq);
        mPhi(h,k) = probcheck(mPhi(h,k));
      }
    }
    
    pg = OmegaCheck(pg, iK);
    // step E
    
    for(n = 0; n < iN; n++){
      vHdY.zeros();
      for(k = 0; k < iK; k++){
        for(h=0; h< iH;h++){
          mdY(n,k,h) = Rf_dbinom(mY(n,h), size, mPhi(h,k), isize);
          vHdY(k) += mdY(n,k,h);
        }
      }
      vLLK(n) = MixtDensityScale(pg, vHdY, iK);
      for(k = 0;k < iK; k++){
        mU(n,k) = exp(log(pg(k)) + vHdY(k) - vLLK(n));
      }
    }
    
    LLKSeries(iter) = accu(vLLK%ivFreq);
    if(iter > 10){
      eps = abs3(LLKSeries(iter) - LLKSeries(iter-1));
    }
    iter +=1;
  }
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  if(reord == 1){
    arma::vec vPhisum = sum(mPhi).t();
    arma::uvec order = sort_index(vPhisum,"descending");
    int ifoo = 0;
    arma::mat mPhi_sorted = mPhi;
    arma::vec pg_sorted = pg;
    arma::mat mU_sorted = mU;
    for(k=0; k< iK; k++){
      ifoo               = order(k);
      mPhi_sorted.col(k) = mPhi.col(ifoo);
      pg_sorted(k)       = pg(ifoo);
      mU_sorted.col(k)   = mU.col(ifoo); 
    }
    mPhi = mPhi_sorted;
    pg = pg_sorted;
    mU = mU_sorted;
  }
  arma::vec alphafoo(iK);
  alphafoo = log(pg/pg(0));
  arma::vec alpha(iK-1);
  alpha = alphafoo.subvec(1,iK-1);
  arma::mat gamma(iH,iK);
  gamma = log(mPhi/(1.0 - mPhi));
  
  arma::vec parvec = join_cols(alpha,vectorise(gamma));
  
  double BIC, AIC;
  BIC = -2.0*LLKSeries(iter-1) + log(iNtot)*(iH*iK + (iK - 1));
  AIC = -2.0*LLKSeries(iter-1) + 2*(iH*iK + (iK - 1));
  double Terr = accu(-pg%log(pg));
  arma::mat mlogU = trunc_log(mU);
  double Perr = sum(sum(-mU%mlogU,1)%ivFreq)/iNtot;
  double R2entr = (Terr - Perr)/Terr;
  
  // classification error
  arma::ivec vModalAssnm(iN);
  arma::mat mModalAssnm = zeros(iN,iK);
  for(n=0; n< iN; n++){
    vModalAssnm(n) = WhichMax(mU.row(n).t());
    mModalAssnm(n,vModalAssnm(n)) = 1.0;
  }
  
  arma::mat mClassErr     = zeros(iK,iK);
  arma::mat mClassErrProb = mClassErr;
  for(j = 0; j < iK; j++){
    for(k = 0; k < iK; k++){
      mClassErr(j,k) = accu(mU.col(k)%mModalAssnm.col(j)%ivFreq);
    }
    mClassErrProb.row(j)=mClassErr.row(j)/accu(mClassErr.row(j));
  }
  
  double dClassErr_tot =  1.0 - (accu(mClassErr.diag())/iNtot);
  
  
  // Computing the score
  
  arma::mat mPg_Score(iN,iK-1);
  for(k =1; k< iK; k++){
    mPg_Score.col(k-1) = (mU.col(k) - pg(k))%ivFreq;
  }
  
  arma::mat mGamma_Score=zeros(iN,iH*iK);
  int iroll = 0;
  for(k = 0; k < iK; k++){
    for(h = 0; h < iH; h++){
      mGamma_Score.col(iroll) = (mU.col(k)%(mY.col(h) - mPhi(h,k)))%ivFreq;
      iroll += 1;
    }
  }
  
  iroll=0;
  
  arma::mat mScore = join_rows(mPg_Score,mGamma_Score);
  arma::mat Infomat = mScore.t()*mScore/iN;
  arma::mat Varmat = pinv(Infomat)/iN;
  arma::vec SEs =  sqrt(Varmat.diag());
  
  List EMout;
  EMout["mU"]            = mU;
  EMout["mPhi"]          = mPhi;
  EMout["pg"]            = pg;
  EMout["alphas"]        = alpha;
  EMout["gamma"]         = gamma;
  EMout["mScore"] =mScore;
  EMout["Varmat"] =Varmat;
  EMout["SEs"] =SEs;
  EMout["LLKSeries"]     = LLKSeries;
  EMout["eps"]           = eps;
  EMout["iter"]          = iter;
  EMout["BIC"]           = BIC;
  EMout["AIC"]           = AIC;
  EMout["R2entr"]        = R2entr;
  EMout["mClassErr"]     = mClassErr;
  EMout["mClassErrProb"] = mClassErrProb;
  EMout["dClassErr_tot"] = dClassErr_tot;
  EMout["mModalAssnm"]   = mModalAssnm;
  EMout["vModalAssnm"]   = vModalAssnm;
  EMout["freq"]          = ivFreq;
  EMout["parvec"]    = parvec;
  
  return EMout;
}

