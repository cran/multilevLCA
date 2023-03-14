library(MASS)
library(dplyr)
library(tidyr)
library(klaR)
library(magrittr)


LCA_fast_init = function(mY, iK, kmea=T, npc=NULL, maxIter = 1e3, tol = 1e-8, reord = 1){
  mY_df     = data.frame(mY)
  group_by_all = NULL
  mY_aggr   = as.matrix(mY_df  %>% group_by_all %>% count)
  iHf       = dim(mY_aggr)[2]
  freq      = mY_aggr[,iHf]
  mY_unique = mY_aggr[,-iHf]
  if(kmea == F){
    clusfoo   = klaR::kmodes(mY_unique,modes=iK)$cluster
  }else{
    if(is.null(npc)){
      prscores = prcomp(mY_unique)$x[,1:2]
    }else{
      prscores = prcomp(mY_unique)$x[,1:npc]
    }
    spectclust = kmeans(prscores,iK)
    clusfoo = spectclust$cluster
  }
  mU = vecTomatClass(clusfoo)
  out = LCA_fast(mY_unique, freq, iK, mU, maxIter, tol, reord)
  return(out)
}
#
LCA_fast_init_wcov = function(mY, mZ, iK, kmea=T, npc=NULL, maxIter = 1e3, tol = 1e-8, fixed = 0, reord = 1,
                              NRtol = 1e-6, NRmaxit = 100){
  # mZ must include a column of ones!
  group_by_all = NULL
  mY_df     = data.frame(mY)
  mY_aggr   = as.matrix(mY_df  %>% group_by_all %>% count)
  iHf       = dim(mY_aggr)[2]
  freq      = mY_aggr[,iHf]
  mY_unique = mY_aggr[,-iHf]
  if(kmea == F){
    clusfoo   = klaR::kmodes(mY_unique,modes=iK)$cluster
  }else{
    if(is.null(npc)){
      prscores = prcomp(mY_unique)$x[,1:2]
    }else{
      prscores = prcomp(mY_unique)$x[,1:npc]
    }
    spectclust = kmeans(prscores,iK)
    clusfoo = spectclust$cluster
    # clusfoo = kmeans(mY_unique,iK)$cluster
  }
  mU = vecTomatClass(clusfoo)
  out        = LCA_fast(mY_unique, freq, iK, mU, maxIter, tol, reord)
  P = ncol(mZ)
  mBeta_init = matrix(0,P,iK-1)
  mBeta_init[1,] = out$alphas
  Step1Var = out$Varmat
  outcov = LCAcov(mY,mZ,iK,out$mPhi,mBeta_init,Step1Var,fixed, maxIter, 
                  tol, NRtol, NRmaxit)
  return(list(out=out,outcov=outcov))
}

# 
LCA_fast_init_whigh = function(mY, id_high, iK, kmea=T, npc=NULL, maxIter = 1e3, tol = 1e-8, reord = 1, debug=F){
  mY_df     = data.frame(mY,id_high)
  group_by_all = NULL
  mY_aggr   = as.matrix(mY_df  %>% group_by_all %>% count)
  iHf       = dim(mY_aggr)[2]
  mY_aggr   = mY_aggr[order(mY_aggr[,iHf-1]),]
  freq      = mY_aggr[,iHf]
  mY_unique = mY_aggr[,-c(iHf-1,iHf)]
  # 
  # mY_unique_df = data.frame(mY_unique)
  # for(var in 1:(dim(mY_unique_df)[2])){
  #   mY_unique_df[,var] = factor(mY_unique_df[,var])
  # }
  # clusfoo = kproto(mY_unique_df,k=iK,nstart=10)
  # # 
  if(kmea==F){
    clusfoo   = try(kmodes(mY_unique,modes=iK),silent = T)
    check = 0
    while((inherits(clusfoo,"try-error"))){ 
      clusfoo   = try(kmodes(mY_unique,modes=iK),silent = T)
      check = check + 1
      if(check > 2 & debug==T) print(check)
    }
    clusfoo = clusfoo$cluster
  }else{
    if(is.null(npc)){
      prscores = prcomp(mY_unique)$x[,1:2]
    }else{
      prscores = prcomp(mY_unique)$x[,1:npc]
    }
    spectclust = kmeans(prscores,iK)
    clusfoo = spectclust$cluster
  }
  
  mU = vecTomatClass(clusfoo)
  out = LCA_fast(mY_unique, freq, iK, mU, maxIter, tol, reord)
  return(out)
}
# 
LCA_fast_init_norecode = function(mY_unique,freq,iK, kmea=T, npc=NULL, maxIter = 1e3, tol = 1e-8, reord = 1){
  if(kmea==F){
    clusfoo = klaR::kmodes(mY_unique,modes=iK)$cluster
  }else{
    if(is.null(npc)){
      prscores = prcomp(mY_unique)$x[,1:2]
    }else{
      prscores = prcomp(mY_unique)$x[,1:npc]
    }
    spectclust = kmeans(prscores,iK)
    clusfoo = spectclust$cluster
  }
  mU = vecTomatClass(clusfoo)
  out = LCA_fast(mY_unique, freq, iK, mU, maxIter, tol, reord)
  return(out)
}
# 
# 
meas_Init = function(mY, id_high, vNj, iM, iT,kmea=T,npc = NULL){
  # fixed number of low- (iT) and high- (iM) level classes
  iJ = length(vNj)
  iN = dim(mY)[1]
  iK = dim(mY)[2]
  # 
  # working out starting values at low level first
  # 
  out_LCA  = try(LCA_fast_init_whigh(mY,id_high = id_high,iK = iT,kmea=kmea),silent = T)
  check = 0
  while((inherits(out_LCA,"try-error")) | any(out_LCA$mPhi<1e-3)){
    out_LCA  = try(LCA_fast_init_whigh(mY,id_high = id_high,iK = iT,kmea=kmea),silent=T)
    check = check + 1
    if(check > 2) break
  }
  # 
  # now turning to higher level
  # 
  lowlev_relclassprop = matrix(0,iJ,iT)
  foo = 0
  for(j in 1:iJ){
    lowlev_relclassprop[j,] = colMeans(out_LCA$mU[foo+(1:vNj[j]),])
    foo = foo +vNj[j]
  }
  high_out = kmeans(lowlev_relclassprop,centers = iM,nstart = 50)
  vOmega_start = high_out$size/iJ 
  Wmodal_mat = vecTomatClass(high_out$cluster)
  # reordering in decreasing order
  highreord = order(vOmega_start,decreasing = T);
  Wmodal_mat_reord = Wmodal_mat[,highreord]
  Wmodal = apply(Wmodal_mat_reord,1,which.max)
  vOmega_start = vOmega_start[highreord]
  index_indiv_high = rep(Wmodal,times=vNj)
  # 
  mPhi_start = out_LCA$mPhi
  mPi_fast = table(index_indiv_high,rep(out_LCA$vModalAssnm+1,times=out_LCA$freq))
  mPi_start = t(mPi_fast/rowSums(mPi_fast))
  # 
  # 
  return(list(vOmega_start=vOmega_start, mPi_start = mPi_start, mPhi_start = mPhi_start))
} 


simultsel_fun = function(data,iM,iT,kmea=T, npc = NULL,
                         maxIter = 1e3, tol = 1e-8, fixedpars=0, reord = 1){
  LCAout=NULL
  mY       = data$Y
  iH       = ncol(mY)
  K        = iH
  vNj      = table(data$idhigh_long-1)
  iN       = length(vNj)
  if(iM ==1 & iT ==1){
    ll = sum(apply(mY,2,function(x){dbinom(x,1,mean(x),log=T)}))
    BIClow  = -2*ll + K*log(sum(vNj))
    BIChigh = -2*ll + K*log(iN)
    AIC     = -2*ll + 2*K
    ICL_BIClow = Inf
    ICL_BIChigh = Inf
    outmuLCA=list()
  }else if(iM > 1 & iT == 1){
    ll = -Inf
    BIClow  = Inf
    BIChigh = Inf
    AIC     = Inf
    ICL_BIClow = Inf
    ICL_BIChigh = Inf
    outmuLCA=list()
  }else if(iM==1 & iT>1){
    LCAout = LCA_fast_init(mY, iK=iT)
    ll      = tail(LCAout$LLKSeries,1)
    npar    = iT*iH + iT - 1
    BIClow  = LCAout$BIC
    BIChigh = -2*ll + npar*log(iN)
    AIC     = -2*ll + 2*npar
    ICL_BIClow = BIClow + 2*sum(LCAout$freq*apply(-LCAout$mU*log(LCAout$mU),1,sum))
    ICL_BIChigh = Inf
    outmuLCA=LCAout
  }else{
    # actual multilevel LC fit
    start    = try(meas_Init(mY,id_high=data$idhigh_long,vNj,iM=iM,iT=iT,kmea=kmea,npc=npc),silent=T)
    if(inherits(start,"try-error")){
      ll      = -Inf 
      BIClow  = Inf
      BIChigh = Inf
      AIC     = Inf
      ICL_BIClow = Inf
      ICL_BIChigh = Inf
      outmuLCA=list()
    }else{
      vOmegast = start$vOmega_start
      mPi = start$mPi_start
      outmuLCA = try(MLTLCA(mY, vNj, vOmegast, mPi, start$mPhi_start,
                            maxIter = maxIter, tol = tol, reord = reord),silent=T)
      if(inherits(outmuLCA,"try-error")){
        ll      = -Inf 
        BIClow  = Inf
        BIChigh = Inf
        AIC     = Inf
        ICL_BIClow = Inf
        ICL_BIChigh = Inf
      }else{
        ll      = tail(outmuLCA$LLKSeries,1) 
        BIClow  = outmuLCA$BIClow
        BIChigh = outmuLCA$BIChigh
        AIC     = outmuLCA$AIC
        ICL_BIClow = outmuLCA$ICL_BIClow
        ICL_BIChigh = outmuLCA$ICL_BIChigh
      }
    }
    
  }
  return(list(modFIT = outmuLCA,ll=ll, BIClow = BIClow, BIChigh = BIChigh, AIC = AIC,
              ICL_BIClow = ICL_BIClow, ICL_BIChigh = ICL_BIChigh))
}
# 
lukosel_fun = function(data,iM_max,iT_max,kmea=T,
                       maxIter = 1e3, tol = 1e-8, fixedpars=0, reord = 1){
  mY       = data$Y
  iH       = ncol(mY)
  K        = iH
  vNj      = table(data$idhigh_long-1)
  iN       = length(vNj)
  # 
  # step 1 - lower level
  # 
  LCAout = list()
  BIClow_step1  = rep(0,iT_max)
  BIChigh_step1 = rep(0,iT_max)
  AIC_step1     = rep(0,iT_max)
  ICL_BIClow_step1    = rep(0,iT_max)
  ICL_BIChigh_step1   = rep(0,iT_max)
  
  llfoo = sum(apply(mY,2,function(x){dbinom(x,1,mean(x),log=T)}))
  BIClow_step1[1]      = -2*llfoo + K*log(sum(vNj))
  BIChigh_step1[1]     = -2*llfoo + K*log(iN)
  AIC_step1[1]         = -2*llfoo + 2*K
  ICL_BIClow_step1[1]  = Inf
  ICL_BIChigh_step1[1] = Inf
  for(iT in 2:iT_max){
    # 
    LCAout = LCA_fast_init(mY, iK=iT,kmea=kmea)
    ll      = tail(LCAout$LLKSeries,1)
    npar    = iT*iH + iT - 1
    BIClow_step1[iT]     = LCAout$BIC
    BIChigh_step1[iT]    = -2*ll + npar*log(iN)
    AIC_step1[iT]        = -2*ll + 2*npar
    ICL_BIClow_step1[iT] = BIClow_step1[iT] + 2*sum(LCAout$freq*apply(-LCAout$mU*log(LCAout$mU),1,sum))
    ICL_BIChigh_step1[iT] = Inf
    if(BIClow_step1[iT] < BIClow_step1[iT-1] ) LCAout_store = LCAout
  }
  iT_currbest = which.min(BIClow_step1)
  # step 2 - higher level
  # for iM = 1
  BIClow_step2  = rep(NA,iM_max)
  BIChigh_step2 = rep(NA,iM_max)
  AIC_step2     = rep(NA,iM_max)
  ICL_BIClow_step2    = rep(NA,iM_max)
  ICL_BIChigh_step2   = rep(NA,iM_max)
  BIClow_step2[1] = BIClow_step1[iT_currbest]
  BIChigh_step2[1] = BIChigh_step1[iT_currbest]
  AIC_step2[1]     = AIC_step1[iT_currbest]
  ICL_BIClow_step2[1]    = ICL_BIClow_step1[iT_currbest]
  ICL_BIChigh_step2[1]   = ICL_BIClow_step1[iT_currbest]
  # iM > 1
  if(iT_currbest > 1){
    outmuLCA = list()
    for(iM in 2:iM_max){
      start    = try(meas_Init(mY,id_high=data$idhigh_long,vNj,iM=iM,iT=iT_currbest,kmea=kmea))
      vOmegast = start$vOmega_start
      mPi = start$mPi_start
      
      
      outmuLCA[[iM-1]] = try(MLTLCA(mY, vNj, vOmegast, mPi, start$mPhi_start),silent=T)
      if(inherits(outmuLCA[[iM-1]],"try-error")){
        ll      = -Inf 
        BIClow_step2[iM] = Inf
        BIChigh_step2[iM] = Inf
        AIC_step2[iM]     = Inf
        ICL_BIClow_step2[iM]    = Inf
        ICL_BIChigh_step2[iM]   = Inf
      }else{
        ll      = tail(outmuLCA[[iM-1]]$LLKSeries,1) 
        BIClow_step2[iM] = outmuLCA[[iM-1]]$BIClow
        BIChigh_step2[iM] = outmuLCA[[iM-1]]$BIChigh
        AIC_step2[iM]     = outmuLCA[[iM-1]]$AIC
        ICL_BIClow_step2[iM]    = outmuLCA[[iM-1]]$ICL_BIClow
        ICL_BIChigh_step2[iM]   = outmuLCA[[iM-1]]$ICL_BIChigh
      }
    }
  }
  iM_currbest = which(BIChigh_step2==min(BIChigh_step2,na.rm=T))[1]
  outmuLCA_step2 = outmuLCA[[iM_currbest-1]]
  # step 3 - revisiting lower level
  BIClow_step3  = rep(NA,iT_max-1)
  BIChigh_step3 = rep(NA,iT_max-1)
  AIC_step3     = rep(NA,iT_max-1)
  ICL_BIClow_step3    = rep(NA,iT_max-1)
  ICL_BIChigh_step3   = rep(NA,iT_max-1)
  if(iM_currbest > 1){
    outmuLCA3 = list()
    for(iT in 1:(iT_max-1)){
      start    = try(meas_Init(mY,id_high=data$idhigh_long,vNj,iM=iM_currbest,iT=(iT+1),kmea=kmea))
      vOmegast = start$vOmega_start
      mPi = start$mPi_start
      
      outmuLCA3[[iT]] = try(MLTLCA(mY, vNj, vOmegast, mPi, start$mPhi_start,
                                   maxIter = maxIter, tol = tol, reord = reord),silent=T)
      if(inherits(outmuLCA3[[iT]],"try-error")){
        BIClow_step3[iT] = Inf
        BIChigh_step3[iT] = Inf
        AIC_step3[iT]     = Inf
        ICL_BIClow_step3[iT]    = Inf
        ICL_BIChigh_step3[iT]   = Inf  
      }else{
        ll      = tail(outmuLCA3[[iT]]$LLKSeries,1) 
        BIClow_step3[iT] = outmuLCA3[[iT]]$BIClow
        BIChigh_step3[iT] = outmuLCA3[[iT]]$BIChigh
        AIC_step3[iT]     = outmuLCA3[[iT]]$AIC
        ICL_BIClow_step3[iT]    = outmuLCA3[[iT]]$ICL_BIClow
        ICL_BIChigh_step3[iT]   = outmuLCA3[[iT]]$ICL_BIChigh
      }
    }
    iT_currbest = which(BIClow_step3==min(BIClow_step3,na.rm=T))+1
    outmuLCA_step3 = outmuLCA3[[iT_currbest-1]]
  }
  return(list(outmuLCA_step3=outmuLCA_step3,iT_opt=iT_currbest, iM_opt=iM_currbest, BIClow_step1 = BIClow_step1, BIChigh_step1 = BIChigh_step1, AIC_step1 = AIC_step1,
              ICL_BIClow_step1 = ICL_BIClow_step1, ICL_BIChigh_step1 = ICL_BIChigh_step1,
              BIClow_step2 = BIClow_step2, BIChigh_step2 = BIChigh_step2, AIC_step2 = AIC_step2,
              ICL_BIClow_step2 = ICL_BIClow_step2, ICL_BIChigh_step2 = ICL_BIChigh_step2,
              BIClow_step3 = BIClow_step3, BIChigh_step3 = BIChigh_step3, AIC_step3 = AIC_step3,
              ICL_BIClow_step3 = ICL_BIClow_step3, ICL_BIChigh_step3 = ICL_BIChigh_step3))
}

####################
### JOHAN'S CODE ###
####################

simultsel = function(data,Y,iT,id_high,iM){
  
  if(iT==1&iM==1){
    
    iH = ncol(data[,Y])
    iN = nrow(data)
    iJ = length(table(data[,id_high])[table(data[,id_high])>0])
    
    ll = sum(apply(data[,Y],2,function(x){dbinom(x,1,mean(x),log=TRUE)}))
    
    AIC         = -2*ll+2*iH
    BIClow      = -2*ll+iH*log(iN)
    BIChigh     = -2*ll+iH*log(iJ)
    ICL_BIClow  = NA
    ICL_BIChigh = NA
    
    modFIT = list()
    
  } else if(iM>1&iT==1){
    
    ll          = NA
    AIC         = NA
    BIClow      = NA
    BIChigh     = NA
    ICL_BIClow  = NA
    ICL_BIChigh = NA
    
    modFIT = list()
    
  } else if(iM==1&iT>1){
    
    est = estlca(data=data,Y=Y,iT=iT,id_high=NULL,iM=NULL,Z=NULL,Zh=NULL,
                 extout=TRUE,dataout=FALSE,kmea=TRUE,
                 maxIter=1e3,tol=1e-8,fixed=TRUE,reord=TRUE,
                 fixedpars=1,NRtol=1e-6,NRmaxit=100,freqout=TRUE)
    
    iH    = ncol(data[,Y])
    iJ    = length(table(data[,id_high])[table(data[,id_high])>0])
    npar  = iT*iH+iT-1
    
    ll = tail(est$LLKSeries,1)
    
    AIC         = -2*ll+2*npar
    BIClow      = est$BIC
    BIChigh     = -2*ll+npar*log(iJ)
    ICL_BIClow  = est$BIC+2*sum(est$freq*apply(-est$mU*log(est$mU),1,sum))
    ICL_BIChigh = NA
    
    modFIT = est
    
  } else{
    
    est = estlca(data=data,Y=Y,iT=iT,id_high=id_high,iM=iM,Z=NULL,Zh=NULL,
                 extout=FALSE,dataout=FALSE,kmea=TRUE,
                 maxIter=1e3,tol=1e-8,fixed=TRUE,reord=TRUE,
                 fixedpars=1,NRtol=1e-6,NRmaxit=100)
    
    ll          = tail(est$LLKSeries,1)
    AIC         = est$AIC
    BIClow      = est$BIClow
    BIChigh     = est$BIChigh
    ICL_BIClow  = est$ICL_BIClow
    ICL_BIChigh = est$ICL_BIChigh
    
    modFIT = est
    
  }
  
  return(list(modFIT=modFIT,ll=ll,AIC=AIC,BIClow=BIClow,BIChigh=BIChigh,
              ICL_BIClow=ICL_BIClow,ICL_BIChigh=ICL_BIChigh))
  
}

estlca = function(data,Y,iT,id_high,iM,Z,Zh,extout,dataout,kmea,maxIter,tol,fixed,reord,fixedpars,NRtol,NRmaxit,freqout=FALSE){
  
  if(is.null(id_high) & is.null(iM) & is.null(Z) & is.null(Zh)){
    
    ######################################
    ### Single-level measurement model ###
    ######################################
    
    mY = data[,Y]
    iH = ncol(mY)
    
    # Initialization
    mY_df     = data.frame(mY)
    mY_aggr   = as.matrix(mY_df%>%group_by_all%>%count)
    iHf       = dim(mY_aggr)[2]
    freq      = mY_aggr[,iHf]
    mY_unique = mY_aggr[,-iHf]
    if(kmea==FALSE){
      
      clusfoo     = klaR::kmodes(mY_unique,modes=iT)$cluster
      
    } else{
      
      prscores    = prcomp(mY_unique)
      num_pc      = max(round(ncol(mY_unique)/2),(sum(cumsum(prscores$sdev^2/sum(prscores$sdev^2))<0.85)+1))
      prscores    = prscores$x[,1:num_pc]
      spectclust  = kmeans(prscores,centers=iT,iter.max=100,nstart=100)
      clusfoo     = spectclust$cluster
      
    }
    mU        = vecTomatClass(clusfoo)
    
    # Estimation
    est = LCA_fast(mY_unique,freq,iT,mU,maxIter,tol,reord)
    
    # Clean output
    if(extout == FALSE){
      
      # Create empty list for output
      out = list()
      
      # Add vector of class proportions
      out$vPg           = est$pg
      rownames(out$vPg) = paste0("P(C",1:iT,")")
      colnames(out$vPg) = ""
      
      # Add matrix of conditional response probabilities
      out$mPhi            = est$mPhi
      rownames(out$mPhi)  = paste0("P(",Y,"|C)")
      colnames(out$mPhi)  = paste0("C",1:iT)
      
      # Add Akaike information criterion
      out$AIC = est$AIC
      
      # Add Bayesian information criterion
      out$BIC = est$BIC
      
      # Add average proportion of classification errors
      out$AvgClassErrProb = est$dClassErr_tot
      
      # Add entropy R-sqr
      out$R2entr = est$R2entr
      
      # Add number of iterations
      out$iter = est$iter
      
      # Add log-likelihood series
      out$LLKSeries = est$LLKSeries
      
      # Add specification
      out$spec            = as.matrix("Single-level LC measurement model")
      rownames(out$spec)  = colnames(out$spec) = ""
      
    } else{
      
      # Create empty list for output
      out = list()
      
      # Add vector of class proportions
      out$vPg           = est$pg
      rownames(out$vPg) = paste0("P(C",1:iT,")")
      colnames(out$vPg) = ""
      
      # Add matrix of conditional response probabilities
      out$mPhi            = est$mPhi
      rownames(out$mPhi)  = paste0("P(",Y,"|C)")
      colnames(out$mPhi)  = paste0("C",1:iT)
      
      # Add alphas
      out$alphas            = est$alphas
      rownames(out$alphas)  = paste0("alpha(C",2:iT,")")
      colnames(out$alphas)  = ""
      
      # Add gammas
      out$gammas            = est$gamma
      rownames(out$gammas)  = paste0("gamma(",Y,"|C)")
      colnames(out$gammas)  = colnames(out$mPhi)
      
      # Add vector of model parameters
      out$parvec            = est$parvec
      rownames(out$parvec)  = c(rownames(out$alphas),paste0(rep(substr(rownames(out$gammas),1,nchar(rownames(out$gammas))-2),iT),rep(colnames(out$gammas),rep(iH,iT)),")"))
      colnames(out$parvec)  = ""
      
      # Add variance-covariance matrix
      out$Varmat            = est$Varmat
      rownames(out$Varmat)  = colnames(out$Varmat) = rownames(out$parvec)
      
      # Add vector of standard errors
      out$SEs           = est$SEs
      rownames(out$SEs) = rownames(out$parvec)
      colnames(out$SEs) = ""
      
      # Add epsilon
      out$eps = est$eps
      
      # Add matrix of posterior class membership probabilities
      out$mU            = est$mU
      colnames(out$mU)  = colnames(out$mPhi)
      
      # Add matrix of modal class assignment
      out$mU_modal            = est$mModalAssnm
      colnames(out$mU_modal)  = colnames(out$mPhi)
      
      # Add vector of modal class assignment
      out$vU_modal = est$vModalAssnm
      
      # Add matrix of classification errors
      out$mClassErr           = est$mClassErr
      rownames(out$mClassErr) = paste0("C",1:iT,"_true")
      colnames(out$mClassErr) = paste0("C",1:iT,"_post")
      
      # Add matrix of proportion of classification errors
      out$mClassErrProb           = est$mClassErrProb
      rownames(out$mClassErrProb) = rownames(out$mClassErr)
      colnames(out$mClassErrProb) = colnames(out$mClassErr)
      
      # Add average proportion of classification errors
      out$AvgClassErrProb = est$dClassErr_tot
      
      # Add entropy R-sqr
      out$R2entr = est$R2entr
      
      # Add Akaike information criterion
      out$AIC = est$AIC
      
      # Add Bayesian information criterion
      out$BIC = est$BIC
      
      # Add number of iterations
      out$iter = est$iter
      
      # Add log-likelihood series
      out$LLKSeries = est$LLKSeries
      
      # Add matrix of model parameter contributions to log-likelihood score
      out$mScore            = est$mScore
      colnames(out$mScore)  = rownames(out$parvec)
      
      # Add specification
      out$spec            = as.matrix("Single-level LC measurement model")
      rownames(out$spec)  = colnames(out$spec) = ""
      
      if(freqout==TRUE){
        
        out$freq = est$freq
        
      }
      
    }
    if(dataout == TRUE){
      
      out$data            = mY_unique
      rownames(out$data)  = NULL
      
    }
    
    # Return
    return(out)
    
  } else if(is.null(id_high) & is.null(iM) & !is.null(Z) & is.null(Zh)){
    
    #####################################
    ### Single-level structural model ###
    #####################################
    
    mY        = data[,Y]
    mY        = as.matrix(mY)
    mZ        = as.data.frame(data[,Z])
    names(mZ) = Z
    if(FALSE%in%unlist(lapply(as.data.frame(mZ),is.numeric))){
      
      mZkeep            = as.matrix(mZ[,unlist(lapply(as.data.frame(mZ),is.numeric))])
      colnames(mZkeep)  = colnames(mZ)[unlist(lapply(as.data.frame(mZ),is.numeric))]
      mZclean           = as.matrix(mZ[,!unlist(lapply(as.data.frame(mZ),is.numeric))])
      mZclean_names     = colnames(mZ)[!unlist(lapply(as.data.frame(mZ),is.numeric))]
      if(ncol(mZclean) == 1){
        
        mZclean = as.matrix(as.character(mZclean))
        
      } else{
        
        mZclean = apply(mZclean,2,function(x){as.character(x)})
        
      }
      mZclean[mZclean==""] = NA
      for(i in 1:ncol(mZclean)){
        
        for(j in 2:length(unique(na.omit(mZclean[,i])))){
          
          lvl               = sort(unique(na.omit(mZclean[,i])))[j]
          cleaned           = as.matrix(ifelse(mZclean[,i] == lvl,1,0))
          colnames(cleaned) = paste0(mZclean_names[i],".",lvl)
          mZkeep            = cbind(mZkeep, cleaned)
          
        }
        
      }
      mZ = as.matrix(apply(mZkeep,2,function(x){as.numeric(x)}))
      
      
    }
    mZ        = cbind(1,mZ)
    mZ        = as.matrix(mZ)
    Z         = c("Intercept",colnames(mZ)[-1])
    iH        = ncol(mY)
    
    # Initialization
    mY_df           = data.frame(mY)
    mY_aggr         = as.matrix(mY_df%>%group_by_all%>%count)
    iHf             = dim(mY_aggr)[2]
    freq            = mY_aggr[,iHf]
    mY_unique       = mY_aggr[,-iHf]
    if(kmea==FALSE){
      
      clusfoo     = klaR::kmodes(mY_unique,modes=iT)$cluster
      
    } else{
      
      prscores    = prcomp(mY_unique)
      num_pc      = max(round(ncol(mY_unique)/2),(sum(cumsum(prscores$sdev^2/sum(prscores$sdev^2))<0.85)+1))
      prscores    = prscores$x[,1:num_pc]
      spectclust  = kmeans(prscores,centers=iT,iter.max=100,nstart=100)
      clusfoo     = spectclust$cluster
      
    }
    mU              = vecTomatClass(clusfoo)
    outinit         = LCA_fast(mY_unique,freq,iT,mU,maxIter,tol,reord)
    P               = ncol(mZ)
    mBeta_init      = matrix(0,P,iT-1)
    mBeta_init[1,]  = outinit$alphas
    Step1Var        = outinit$Varmat
    mPhi_start      = outinit$mPhi
    
    # Remove missing values on mZ, update mY
    nomissing = complete.cases(mZ)
    mY        = mY[nomissing,]
    mZ        = mZ[nomissing,]
    
    # Estimation
    est = LCAcov(mY,mZ,iT,mPhi_start,mBeta_init,Step1Var,fixed,maxIter,tol,NRtol,NRmaxit)
    
    # Clean output
    if(extout == FALSE){
      
      # Create empty list for output
      out = list()
      
      # Add vector of sample means of class proportions
      out$vPg_avg           = as.matrix(apply(est$mPg,2,mean))
      rownames(out$vPg_avg) = paste0("P(C",1:iT,")")
      colnames(out$vPg_avg) = ""
      
      # Add matrix of conditional response probabilities
      out$mPhi            = est$mPhi
      rownames(out$mPhi)  = paste0("P(",Y,"|C)")
      colnames(out$mPhi)  = paste0("C",1:iT)
      
      # Add betas
      out$betas           = est$beta
      rownames(out$betas) = paste0("beta(",Z,"|C)")
      colnames(out$betas) = paste0("C",2:iT)
      
      # Add corrected standard errors for beta
      out$SEs_cor_beta            = matrix(est$SEs_cor[1:((iT-1)*P),],P,iT-1)
      rownames(out$SEs_cor_beta)  = rownames(out$betas)
      colnames(out$SEs_cor_beta)  = colnames(out$betas)
      
      # Add Akaike information criterion
      out$AIC = est$AIC
      
      # Add Bayesian information criterion
      out$BIC = est$BIC
      
      # Add average proportion of classification errors
      out$AvgClassErrProb = est$dClassErr_tot
      
      # Add entropy R-sqr
      out$R2entr = est$R2entr
      
      # Add number of iterations
      out$iter = est$iter
      
      # Add log-likelihood series
      out$LLKSeries = est$LLKSeries
      
      # Add specification
      out$spec            = as.matrix("Single-level LC structural model")
      rownames(out$spec)  = colnames(out$spec) = ""
      
    } else{
      
      # Create empty list for output
      out = list()
      
      # Add matrix of class proportions
      out$mPg           = est$mPg
      colnames(out$mPg) = paste0("P(C",1:iT,")")
      
      # Add vector of sample means of class proportions
      out$vPg_avg           = as.matrix(apply(est$mPg,2,mean))
      rownames(out$vPg_avg) = colnames(out$mPg)
      colnames(out$vPg_avg) = ""
      
      # Add matrix of conditional response probabilities
      out$mPhi            = est$mPhi
      rownames(out$mPhi)  = paste0("P(",Y,"|C)")
      colnames(out$mPhi)  = paste0("C",1:iT)
      
      # Add gammas
      out$gammas            = est$gamma
      rownames(out$gammas)  = paste0("gamma(",Y,"|C)")
      colnames(out$gammas)  = colnames(out$mPhi)
      
      # Add betas
      out$betas           = est$beta
      rownames(out$betas) = paste0("beta(",Z,"|C)")
      colnames(out$betas) = paste0("C",2:iT)
      
      # Add vector of model parameters
      out$parvec            = est$parvec
      rownames(out$parvec)  = c(paste0(rep(substr(rownames(out$betas),1,nchar(rownames(out$betas))-2),iT-1),rep(colnames(out$betas),rep(P,iT-1)),")"),paste0(rep(substr(rownames(out$gammas),1,nchar(rownames(out$gammas))-2),iT),rep(colnames(out$gammas),rep(iH,iT)),")"))
      colnames(out$parvec)  = ""
      
      # Add inverse of the information matrix from the second step
      out$mV2           = est$mV2
      rownames(out$mV2) = colnames(out$mV2) = paste0(rep(substr(rownames(out$betas),1,nchar(rownames(out$betas))-2),iT-1),rep(colnames(out$betas),rep(P,iT-1)),")")
      
      # Add mQ
      out$mQ            = est$mQ
      rownames(out$mQ)  = colnames(out$mQ) = rownames(out$mV2)
      
      # Add uncorrected variance-covariance matrix
      out$Varmat_unc            = est$Varmat_unc[1:((iT-1)*P),1:((iT-1)*P)]
      rownames(out$Varmat_unc)  = colnames(out$Varmat_unc) = rownames(out$mV2)
      
      # Add corrected variance-covariance matrix
      out$Varmat_cor            = est$Varmat_cor
      rownames(out$Varmat_cor)  = colnames(out$Varmat_cor) = rownames(out$mV2)
      
      # Add vector of uncorrected standard errors
      out$SEs_unc           = as.matrix(est$SEs_unc[1:((iT-1)*P),])
      rownames(out$SEs_unc) = rownames(out$mV2)
      colnames(out$SEs_unc) = ""
      
      # Add vector of corrected standard errors
      out$SEs_cor           = as.matrix(est$SEs_cor[1:((iT-1)*P),])
      rownames(out$SEs_cor) = rownames(out$mV2)
      colnames(out$SEs_cor) = ""
      
      # Add corrected standard errors for beta
      out$SEs_cor_beta            = matrix(est$SEs_cor[1:((iT-1)*P),],P,iT-1)
      rownames(out$SEs_cor_beta)  = rownames(out$betas)
      colnames(out$SEs_cor_beta)  = colnames(out$betas)
      
      # Add epsilon
      out$eps = est$eps
      
      # Add matrix of posterior class membership probabilities
      out$mU            = est$mU
      colnames(out$mU)  = colnames(out$mPhi)
      
      # Add matrix of classification errors
      out$mClassErr           = est$mClassErr
      rownames(out$mClassErr) = paste0("C",1:iT,"_true")
      colnames(out$mClassErr) = paste0("C",1:iT,"_post")
      
      # Add matrix of proportion of classification errors
      out$mClassErrProb           = est$mClassErrProb
      rownames(out$mClassErrProb) = paste0("C",1:iT,"_true")
      colnames(out$mClassErrProb) = paste0("C",1:iT,"_post")
      
      # Add average proportion of classification errors
      out$AvgClassErrProb = est$dClassErr_tot
      
      # Add entropy R-sqr
      out$R2entr = est$R2entr
      
      # Add Akaike information criterion
      out$AIC = est$AIC
      
      # Add Bayesian information criterion
      out$BIC = est$BIC
      
      # Add number of iterations
      out$iter = est$iter
      
      # Add log-likelihood series
      out$LLKSeries = est$LLKSeries
      
      # Add specification
      out$spec            = as.matrix("Single-level LC structural model")
      rownames(out$spec)  = colnames(out$spec) = ""
      
    }
    if(dataout == TRUE){
      
      out$data            = cbind(mY,mZ)
      rownames(out$data)  = NULL
      colnames(out$data)  = c(Y,Z)
      
    }
    
    # Return
    return(out)
    
  } else if(!is.null(id_high) & !is.null(iM) & is.null(Z) & is.null(Zh)){
    
    ####################################
    ### Multilevel measurement model ###
    ####################################
    
    data    = data[order(data[,id_high]),]
    mY      = data[,Y]
    mY      = as.matrix(mY)
    iH      = ncol(mY)
    id_high = data[,id_high]
    idnames = as.character(sort(unique(id_high)))
    id_high = as.numeric(factor(id_high))
    vNj     = as.vector(table(id_high))
    iJ      = length(vNj)
    
    # Initialization
    mY_df             = data.frame(mY,id_high)
    mY_aggr           = as.matrix(mY_df%>%group_by_all%>%count)
    iHf               = dim(mY_aggr)[2]
    mY_aggr           = mY_aggr[order(mY_aggr[,iHf-1]),]
    freq              = mY_aggr[,iHf]
    mY_unique         = mY_aggr[,-c(iHf-1,iHf)]
    if(kmea==FALSE){
      
      clusfoo     = klaR::kmodes(mY_unique,modes=iT)$cluster
      
    } else{
      
      prscores    = prcomp(mY_unique)
      num_pc      = max(round(ncol(mY_unique)/2),(sum(cumsum(prscores$sdev^2/sum(prscores$sdev^2))<0.85)+1))
      prscores    = prscores$x[,1:num_pc]
      spectclust  = kmeans(prscores,centers=iT,iter.max=100,nstart=100)
      clusfoo     = spectclust$cluster
      
    }
    mU                = vecTomatClass(clusfoo)
    outinit           = LCA_fast(mY_unique,freq,iT,mU,maxIter,tol,reord)
    mU_full           = vecTomatClass(rep(outinit$vModalAssnm+1,outinit$freq))
    vUmodal           = rep(0,iJ)
    mUmodal           = matrix(0,iJ,iT)
    foo               = 0
    for(j in 1:iJ){
      
      if(vNj[j]>1){
        
        mUmodal[j,] = colSums(mU_full[foo+(1:vNj[j]),])/vNj[j]
        foo         = foo+vNj[j]
        
      } else{
        
        mUmodal[j,] = colSums(t(mU_full[foo+1,]))/vNj[j]
        foo         = foo+vNj[j]
        
      }
      
      
    }
    high_out          = kmeans(mUmodal,centers=iM,iter.max=100,nstart=100)
    vOmega_start      = high_out$size/iJ
    mIndex_indiv_high = vecTomatClass(rep(high_out$cluster,vNj))
    highreord         = order(vOmega_start,decreasing=TRUE)
    vOmega_start      = vOmega_start[highreord]
    index_indiv_high  = apply(mIndex_indiv_high[,highreord],1,which.max)
    mPhi_start        = outinit$mPhi
    mPi_fast          = table(index_indiv_high,rep(outinit$vModalAssnm+1,outinit$freq))
    mPi_start         = t(mPi_fast/rowSums(mPi_fast))
    
    # Estimation
    est = MLTLCA(mY,vNj,vOmega_start,mPi_start,mPhi_start,maxIter,tol,reord)
    
    # Clean output
    if(extout == FALSE){
      
      # Create empty list for output
      out = list()
      
      # Add matrix of high-level class proportions
      out$vOmega           = est$vOmega
      rownames(out$vOmega) = paste0("P(G",1:iM,")")
      colnames(out$vOmega) = ""
      
      # Add matrix of conditional low-level class proportions
      out$mPi           = est$mPi
      rownames(out$mPi) = paste0("P(C",1:iT,"|G)")
      colnames(out$mPi) = paste0("G",1:iM)
      
      # Add matrix of conditional response probabilities
      out$mPhi            = est$mPhi
      rownames(out$mPhi)  = paste0("P(",Y,"|C)")
      colnames(out$mPhi)  = paste0("C",1:iT)
      
      # Add Akaike information criterion
      out$AIC = est$AIC
      
      # Add low-level Bayesian information criterion
      out$BIClow = est$BIClow
      
      # Add high-level Bayesian information criterion
      out$BIChigh = est$BIChigh
      
      # Add low-level integrated completed likelihood Bayesian information criterion
      out$ICL_BIClow = est$ICL_BIClow
      
      # Add high-level integrated completed likelihood Bayesian information criterion
      out$ICL_BIChigh = est$ICL_BIChigh
      
      # Add low-level entropy R-sqr
      out$R2entr_low = est$R2entr_low
      
      # Add high-level entropy R-sqr
      out$R2entr_high = est$R2entr_high
      
      # Add number of iterations
      out$iter = est$iter
      
      # Add log-likelihood series
      out$LLKSeries = est$LLKSeries
      
      # Add specification
      out$spec            = as.matrix("Multilevel LC measurement model")
      rownames(out$spec)  = colnames(out$spec) = ""
      
    } else{
      
      # Create empty list for output
      out = list()
      
      # Add matrix of high-level class proportions
      out$vOmega           = est$vOmega
      rownames(out$vOmega) = paste0("P(G",1:iM,")")
      colnames(out$vOmega) = ""
      
      # Add matrix of conditional low-level class proportions
      out$mPi           = est$mPi
      rownames(out$mPi) = paste0("P(C",1:iT,"|G)")
      colnames(out$mPi) = paste0("G",1:iM)
      
      # Add matrix of conditional response probabilities
      out$mPhi            = est$mPhi
      rownames(out$mPhi)  = paste0("P(",Y,"|C)")
      colnames(out$mPhi)  = paste0("C",1:iT)
      
      # Add deltas
      out$vDelta            = est$vDelta
      rownames(out$vDelta)  = paste0("delta(G",2:iM,")")
      colnames(out$vDelta)  = ""
      
      # Add gammas
      out$mGamma            = est$mGamma
      rownames(out$mGamma)  = paste0("gamma(C",2:iT,"|G)")
      colnames(out$mGamma)  = colnames(out$mPi)
      
      # Add betas
      out$mBeta           = est$mBeta
      rownames(out$mBeta) = paste0("beta(",Y,"|C)")
      colnames(out$mBeta) = colnames(out$mPhi)
      
      # Add vector of model parameters
      out$parvec            = est$parvec
      rownames(out$parvec)  = c(rownames(out$vDelta),paste0(rep(substr(rownames(out$mGamma),1,nchar(rownames(out$mGamma))-2),iM),rep(colnames(out$mGamma),rep(iT-1,iM)),")"),paste0(rep(substr(rownames(out$mBeta),1,nchar(rownames(out$mBeta))-2),iT),rep(colnames(out$mBeta),rep(iH,iT)),")"))
      colnames(out$parvec)  = ""
      
      # Add Fisher information matrix
      out$Infomat           = est$Infomat
      rownames(out$Infomat) = colnames(out$Infomat) = rownames(out$parvec)
      
      # Add variance-covariance matrix
      out$Varmat            = est$Varmat
      rownames(out$Varmat)  = colnames(out$Varmat) = rownames(out$parvec)
      
      # Add vector of standard errors
      out$SEs           = est$SEs
      rownames(out$SEs) = rownames(out$parvec)
      colnames(out$SEs) = ""
      
      # Add epsilon
      out$eps = est$eps
      
      # Add cube of joint posterior class membership probabilities
      out$cPMX = array(est$cPMX,dim(est$cPMX),dimnames=list(NULL, paste0("C",1:iT,",G"), colnames(out$mPi)))
      
      # Add cube of log of joint posterior class membership probabilities
      out$cLogPMX = array(est$cLogPMX,dim(est$cLogPMX),dimnames=dimnames(out$cPMX))
      
      # Add cube of conditional posterior low-level class membership probabilities
      out$cPX = array(est$cPX,dim(est$cPX),dimnames=list(NULL, paste0("C",1:iT,"|G"),colnames(out$mPi)))
      
      # Add cube of log of conditional posterior low-level class membership probabilities
      out$cLogPX = array(est$cLogPX, dim(est$cLogPX),dimnames=dimnames(out$cPX))
      
      # Add matrix of posterior high-level class membership probabilities for low-level units after marginalizing over low-level classes
      out$mSumPX = est$mSumPX
      colnames(out$mSumPX) = colnames(out$mPi)
      
      # Add matrix of posterior high-level class membership probabilities for high-level units
      out$mPW           = est$mPW
      rownames(out$mPW) = idnames
      colnames(out$mPW) = colnames(out$mPi)
      
      # Add matrix of log of posterior high-level class membership probabilities for high-level units
      out$mlogPW            = est$mlogPW
      rownames(out$mlogPW)  = idnames
      colnames(out$mlogPW)  = colnames(out$mPi)
      
      # Add matrix of posterior high-level class membership probabilities for low-level units
      out$mPW_N           = est$mPW_N
      colnames(out$mPW_N) = colnames(out$mPi)
      
      # Add matrix of posterior low-level class membership probabilities for low-level units after marginalizing over high-level classes
      out$mPMsumX           = est$mPMsumX
      colnames(out$mPMsumX) = colnames(out$mPhi)
      
      # Add low-level entropy R-sqr
      out$R2entr_low = est$R2entr_low
      
      # Add high-level entropy R-sqr
      out$R2entr_high = est$R2entr_high
      
      # Add Akaike information criterion
      out$AIC = est$AIC
      
      # Add low-level Bayesian information criterion
      out$BIClow = est$BIClow
      
      # Add high-level Bayesian information criterion
      out$BIChigh = est$BIChigh
      
      # Add low-level integrated completed likelihood Bayesian information criterion
      out$ICL_BIClow = est$ICL_BIClow
      
      # Add high-level integrated completed likelihood Bayesian information criterion
      out$ICL_BIChigh = est$ICL_BIChigh
      
      # Add number of iterations
      out$iter = est$iter
      
      # Add log-likelihood series
      out$LLKSeries = est$LLKSeries
      
      # Add current log-likelihood for high-level units
      out$vLLK            = est$vLLK
      rownames(out$vLLK)  = idnames
      colnames(out$vLLK)  = ""
      
      # Add matrix of model parameter contributions to log-likelihood score
      out$mScore            = est$mScore
      colnames(out$mScore)  = rownames(out$parvec)
      
      # Add specification
      out$spec            = as.matrix("Multilevel LC measurement model")
      rownames(out$spec)  = colnames(out$spec) = ""
      
    }
    if(dataout == TRUE){
      
      out$data            = cbind(mY,id_high)
      rownames(out$data)  = NULL
      colnames(out$data)  = c(Y,colnames(data)[!colnames(data)%in%Y])
      
    }
    
    # Return
    return(out)
    
  } else if(!is.null(id_high) & !is.null(iM) & !is.null(Z) & is.null(Zh)){
    
    #############################################################
    ### Multilevel structural model with low-level covariates ###
    #############################################################
    
    data      = data[order(data[,id_high]),]
    mY        = data[,Y]
    mY        = as.matrix(mY)
    iH        = ncol(mY)
    id_high   = data[,id_high]
    idnames   = as.character(sort(unique(id_high)))
    id_high   = as.numeric(factor(id_high))
    vNj       = as.vector(table(id_high))
    iJ        = length(vNj)
    mZ        = as.data.frame(data[,Z])
    names(mZ) = Z
    if(FALSE%in%unlist(lapply(as.data.frame(mZ),is.numeric))){
      
      mZkeep            = as.matrix(mZ[,unlist(lapply(as.data.frame(mZ),is.numeric))])
      colnames(mZkeep)  = colnames(mZ)[unlist(lapply(as.data.frame(mZ),is.numeric))]
      mZclean           = as.matrix(mZ[,!unlist(lapply(as.data.frame(mZ),is.numeric))])
      mZclean_names     = colnames(mZ)[!unlist(lapply(as.data.frame(mZ),is.numeric))]
      if(ncol(mZclean) == 1){
        
        mZclean = as.matrix(as.character(mZclean))
        
      } else{
        
        mZclean = apply(mZclean,2,function(x){as.character(x)})
        
      }
      mZclean[mZclean==""] = NA
      for(i in 1:ncol(mZclean)){
        
        for(j in 2:length(unique(na.omit(mZclean[,i])))){
          
          lvl               = sort(unique(na.omit(mZclean[,i])))[j]
          cleaned           = as.matrix(ifelse(mZclean[,i] == lvl,1,0))
          colnames(cleaned) = paste0(mZclean_names[i],".",lvl)
          mZkeep            = cbind(mZkeep, cleaned)
          
        }
        
      }
      mZ = as.matrix(apply(mZkeep,2,function(x){as.numeric(x)}))
      
      
    }
    mZ      = cbind(1,mZ)
    mZ      = as.matrix(mZ)
    Z       = c("Intercept",colnames(mZ)[-1])
    P       = ncol(mZ)
    
    # Initialization
    mY_df             = data.frame(mY,id_high)
    mY_aggr           = as.matrix(mY_df%>%group_by_all%>%count)
    iHf               = dim(mY_aggr)[2]
    mY_aggr           = mY_aggr[order(mY_aggr[,iHf-1]),]
    freq              = mY_aggr[,iHf]
    mY_unique         = mY_aggr[,-c(iHf-1,iHf)]
    if(kmea==FALSE){
      
      clusfoo     = klaR::kmodes(mY_unique,modes=iT)$cluster
      
    } else{
      
      prscores    = prcomp(mY_unique)
      num_pc      = max(round(ncol(mY_unique)/2),(sum(cumsum(prscores$sdev^2/sum(prscores$sdev^2))<0.85)+1))
      prscores    = prscores$x[,1:num_pc]
      spectclust  = kmeans(prscores,centers=iT,iter.max=100,nstart=100)
      clusfoo     = spectclust$cluster
      
    }
    mU                = vecTomatClass(clusfoo)
    outinit1          = LCA_fast(mY_unique,freq,iT,mU,maxIter,tol,reord)
    mU_full           = vecTomatClass(rep(outinit1$vModalAssnm+1,outinit1$freq))
    vUmodal           = rep(0,iJ)
    mUmodal           = matrix(0,iJ,iT)
    foo               = 0
    for(j in 1:iJ){
      
      if(vNj[j]>1){
        
        mUmodal[j,] = colSums(mU_full[foo+(1:vNj[j]),])/vNj[j]
        foo         = foo+vNj[j]
        
      } else{
        
        mUmodal[j,] = colSums(t(mU_full[foo+1,]))/vNj[j]
        foo         = foo+vNj[j]
        
      }
      
    }
    high_out          = kmeans(mUmodal,centers=iM,iter.max=100,nstart=100)
    vOmega_start      = high_out$size/iJ
    mIndex_indiv_high = vecTomatClass(rep(high_out$cluster,vNj))
    highreord         = order(vOmega_start,decreasing=TRUE)
    vOmega_start      = vOmega_start[highreord]
    index_indiv_high  = apply(mIndex_indiv_high[,highreord],1,which.max)
    mPhi_start        = outinit1$mPhi
    mPi_fast          = table(index_indiv_high,rep(outinit1$vModalAssnm+1,outinit1$freq))
    mPi_start         = t(mPi_fast/rowSums(mPi_fast))
    outinit2          = MLTLCA(mY,vNj,vOmega_start,mPi_start,mPhi_start,maxIter,tol,reord)
    vOmega_start      = outinit2$vOmega
    cGamma_start      = array(c(rbind(outinit2$mGamma,matrix(0,(iT-1)*(P-1),iM))),c(iT-1,P,iM))
    mPhi_start        = outinit2$mPhi
    mStep1Var         = outinit2$Varmat
    
    # Remove missing values on mZ, update mY and vNj
    nomissing = complete.cases(mZ)
    mY        = mY[nomissing,]
    idnames   = idnames[sort(unique(id_high[nomissing]))]
    id_high   = id_high[nomissing]
    vNj       = as.vector(table(id_high))
    mZ        = mZ[nomissing,]
    
    # Estimation
    est = MLTLCA_cov(mY,mZ,vNj,vOmega_start,cGamma_start,mPhi_start,mStep1Var,maxIter,tol,fixedpars,NRtol,NRmaxit)
    
    # Return
    if(extout == FALSE){
      
      # Create empty list for output
      out = list()
      
      # Add matrix of high-level class proportions
      out$vOmega            = est$vOmega
      rownames(out$vOmega)  = paste0("P(G",1:iM,")")
      colnames(out$vOmega)  = ""
      
      # Add vector of sample means of conditional low-level class proportions
      out$mPi_avg           = apply(est$cPi,3,function(x){apply(x,2,mean)})
      rownames(out$mPi_avg) = paste0("P(C",1:iT,"|G)")
      colnames(out$mPi_avg) = paste0("G",1:iM)
      
      # Add matrix of conditional response probabilities
      out$mPhi            = est$mPhi
      rownames(out$mPhi)  = paste0("P(",Y,"|C)")
      colnames(out$mPhi)  = paste0("C",1:iT)
      
      # Add gammas
      out$cGamma = array(apply(est$cGamma,3,t),c(P,iT-1,iM),dimnames=list(paste0("gamma(",Z,"|C)"),paste0("C",2:iT,"|G"),colnames(out$mPi_avg)))
      
      # Add corrected standard errors for gamma
      out$SEs_cor_gamma = array(apply(array(est$SEs_cor[iM:(iM-1+(P*(iT-1)*iM))],c(iT-1,P,iM)),3,t),dim(out$cGamma),dimnames=dimnames(out$cGamma))
      
      # Add Akaike information criterion
      out$AIC = est$AIC
      
      # Add low-level Bayesian information criterion
      out$BIClow = est$BIClow
      
      # Add high-level Bayesian information criterion
      out$BIChigh = est$BIChigh
      
      # Add low-level integrated completed likelihood Bayesian information criterion
      out$ICL_BIClow = est$ICL_BIClow
      
      # Add high-level integrated completed likelihood Bayesian information criterion
      out$ICL_BIChigh = est$ICL_BIChigh
      
      # Add low-level entropy R-sqr
      out$R2entr_low = est$R2entr_low
      
      # Add high-level entropy R-sqr
      out$R2entr_high = est$R2entr_high
      
      # Add number of iterations
      out$iter = est$iter
      
      # Add log-likelihood series
      out$LLKSeries = est$LLKSeries
      
      # Add specification
      out$spec            = as.matrix("Multilevel structural LC model with low-level covariates")
      rownames(out$spec)  = colnames(out$spec) = ""
      
    } else{
      
      # Create empty list for output
      out = list()
      
      # Add matrix of high-level class proportions
      out$vOmega            = est$vOmega
      rownames(out$vOmega)  = paste0("P(G",1:iM,")")
      colnames(out$vOmega)  = ""
      
      # Add matrix of conditional low-level class proportions
      out$mPi           = matrix(est$cPi,nrow(est$cPi),iT*iM)
      colnames(out$mPi) = paste0(rep(paste0("P(C",1:iT,"|"), iM),rep(paste0("G",1:iM,")"),rep(iT,iM)))
      
      # Add vector of sample means of conditional low-level class proportions
      out$mPi_avg           = apply(est$cPi,3,function(x){apply(x,2,mean)})
      rownames(out$mPi_avg) = paste0("P(C",1:iT,"|G)")
      colnames(out$mPi_avg) = paste0("G",1:iM)
      
      # Add matrix of conditional response probabilities
      out$mPhi            = est$mPhi
      rownames(out$mPhi)  = paste0("P(",Y,"|C)")
      colnames(out$mPhi)  = paste0("C",1:iT)
      
      # Add deltas
      out$vDelta            = est$vDelta
      rownames(out$vDelta)  = paste0("delta(G",2:iM,")")
      colnames(out$vDelta)  = ""
      
      # Add gammas
      out$cGamma = array(apply(est$cGamma,3,t),c(P,iT-1,iM),dimnames=list(paste0("gamma(",Z,"|C)"),paste0("C",2:iT,"|G"),paste0("G",1:iM)))
      
      # Add betas
      out$mBeta           = est$mBeta
      rownames(out$mBeta) = paste0("beta(",Y,"|C)")
      colnames(out$mBeta) = colnames(out$mPhi)
      
      # Add vector of model parameters
      out$parvec            = est$parvec
      rownames(out$parvec)  = c(rownames(out$vDelta),paste0(apply(out$cGamma,3,function(x){paste0(rep(substr(rownames(x),1,nchar(rownames(x))-2),rep(iT-1,P)),rep(substr(colnames(x),1,nchar(colnames(x))-2),P),"|")}),rep(unlist(dimnames(out$cGamma)[3]),rep(P*(iT-1),iM)),")"),paste0(rep(substr(rownames(out$mBeta),1,nchar(rownames(out$mBeta))-2),iT),rep(colnames(out$mBeta),rep(iH,iT)),")"))
      colnames(out$parvec)  = ""
      
      # Add Fisher information matrix
      out$Infomat           = est$Infomat
      rownames(out$Infomat) = rownames(out$parvec)
      
      # Add cube of Fisher information matrix for gamma
      out$cGamma_Info = array(est$cGamma_Info,dim(est$cGamma_Info),dimnames=list(paste0("gamma(",Z,"|C|G)"),paste0("gamma(",Z,"|C|G)"),paste0("C", 2:iT,rep(paste0("|G", 1:iM),rep(iT-1,iM)))))
      
      # Add inverse of the information matrix from the second step
      out$mV2           = est$mV2
      rownames(out$mV2) = colnames(out$mV2) = c(rownames(out$vDelta),paste0(apply(out$cGamma,3,function(x){paste0(rep(substr(rownames(x),1,nchar(rownames(x))-2),rep(iT-1,P)),rep(substr(colnames(x),1,nchar(colnames(x))-2),P),"|")}),rep(unlist(dimnames(out$cGamma)[3]),rep(P*(iT-1),iM)),")"))
      
      # Add mQ
      out$mQ            = est$mQ
      rownames(out$mQ)  = colnames(out$mQ) = rownames(out$mV2)
      
      # Add uncorrected variance-covariance matrix
      out$Varmat_unc            = est$Varmat[1:(iM-1+P*(iT-1)*iM),1:(iM-1+P*(iT-1)*iM)]
      rownames(out$Varmat_unc)  = colnames(out$Varmat_unc) = rownames(out$mV2)
      
      # Add corrected variance-covariance matrix
      out$Varmat_cor            = est$mVar_corr
      rownames(out$Varmat_cor)  = colnames(out$Varmat_cor) = rownames(out$mV2)
      
      # Add vector of uncorrected standard errors
      out$SEs_unc           = as.matrix(est$SEs_unc[1:(iM-1+P*(iT-1)*iM),])
      rownames(out$SEs_unc) = rownames(out$mV2)
      colnames(out$SEs_unc) = ""
      
      # Add vector of corrected standard errors
      out$SEs_cor           = as.matrix(est$SEs_cor[1:(iM-1+P*(iT-1)*iM),])
      rownames(out$SEs_cor) = rownames(out$mV2)
      colnames(out$SEs_cor) = ""
      
      # Add corrected standard errors for gamma
      out$SEs_cor_gamma = array(apply(array(est$SEs_cor[iM:(iM-1+(P*(iT-1)*iM))],c(iT-1,P,iM)),3,t),dim(out$cGamma), dimnames=dimnames(out$cGamma))
      
      # Add epsilon
      out$eps = est$eps
      
      # Add cube of joint posterior class membership probabilities
      out$cPMX = array(est$cPMX,dim(est$cPMX),dimnames=list(NULL,paste0("C",1:iT,",G"),colnames(out$mPi_avg)))
      
      # Add cube of log of joint posterior class membership probabilities
      out$cLogPMX = array(est$cLogPMX, dim(est$cLogPMX), dimnames=dimnames(out$cPMX))
      
      # Add cube of conditional posterior low-level class membership probabilities
      out$cPX = array(est$cPX, dim(est$cPX), dimnames=list(NULL,paste0("C",1:iT,"|G"),colnames(out$mPi_avg)))
      
      # Add cube of log of conditional posterior low-level class membership probabilities
      out$cLogPX = array(est$cLogPX,dim(est$cLogPX),dimnames=dimnames(out$cPX))
      
      # Add matrix of posterior high-level class membership probabilities for low-level units after marginalizing over low-level classes
      out$mSumPX            = est$mSumPX
      colnames(out$mSumPX)  = colnames(out$mPi_avg)
      
      # Add matrix of posterior high-level class membership probabilities for high-level units
      out$mPW           = est$mPW
      rownames(out$mPW) = idnames
      colnames(out$mPW) = colnames(out$mPi_avg)
      
      # Add matrix of log of posterior high-level class membership probabilities for high-level units
      out$mlogPW            = est$mlogPW
      rownames(out$mlogPW)  = idnames
      colnames(out$mlogPW)  = colnames(out$mPi_avg)
      
      # Add matrix of posterior high-level class membership probabilities for low-level units
      out$mPW_N           = est$mPW_N
      colnames(out$mPW_N) = colnames(out$mPi_avg)
      
      # Add matrix of posterior low-level class membership probabilities for low-level units after marginalizing over high-level classes
      out$mPMsumX           = est$mPMsumX
      colnames(out$mPMsumX) = colnames(out$mPhi)
      
      # Add low-level entropy R-sqr
      out$R2entr_low = est$R2entr_low
      
      # Add high-level entropy R-sqr
      out$R2entr_high = est$R2entr_high
      
      # Add Akaike information criterion
      out$AIC = est$AIC
      
      # Add low-level Bayesian information criterion
      out$BIClow = est$BIClow
      
      # Add high-level Bayesian information criterion
      out$BIChigh = est$BIChigh
      
      # Add low-level integrated completed likelihood Bayesian information criterion
      out$ICL_BIClow = est$ICL_BIClow
      
      # Add high-level integrated completed likelihood Bayesian information criterion
      out$ICL_BIChigh = est$ICL_BIChigh
      
      # Add number of iterations
      out$iter = est$iter
      
      # Add log-likelihood series
      out$LLKSeries = est$LLKSeries
      
      # Add current log-likelihood for high-level units
      out$vLLK            = est$vLLK
      rownames(out$vLLK)  = idnames
      colnames(out$vLLK)  = ""
      
      # Add matrix of model parameter contributions to log-likelihood score
      out$mScore            = est$mScore
      colnames(out$mScore)  = rownames(out$parvec)
      
      # Add subset of matrix of model parameter contributions to log-likelihood score for gamma
      out$mGamma_Score            = est$mGamma_Score
      colnames(out$mGamma_Score)  = paste0(apply(out$cGamma,3,function(x){paste0(rep(substr(rownames(x),1,nchar(rownames(x))-2),rep(iT-1,P)),rep(substr(colnames(x),1,nchar(colnames(x))-2),P),"|")}),rep(unlist(dimnames(out$cGamma)[3]),rep(P*(iT-1),iM)),")")
      
      # Add specification
      out$spec            = as.matrix("Multilevel structural LC model with low-level covariates")
      rownames(out$spec)  = colnames(out$spec) = ""
      
    }
    if(dataout == TRUE){
      
      out$data            = cbind(mY,id_high,mZ)
      rownames(out$data)  = NULL
      colnames(out$data)  = c(Y,colnames(data)[!colnames(data)%in%c(Y,Z)],Z)
      
    }
    return(out)
    
  } else if(!is.null(id_high)&!is.null(iM)&!is.null(Z)&!is.null(Zh)){
    
    #######################################################################
    ### Multilevel structural model with low- and high-level covariates ###
    #######################################################################
    
    data    = data[order(data[,id_high]),]
    mY      = data[,Y]
    mY      = as.matrix(mY)
    iH      = ncol(mY)
    id_high = data[,id_high]
    idnames = as.character(sort(unique(id_high)))
    id_high = as.numeric(factor(id_high))
    vNj     = as.vector(table(id_high))
    iJ      = length(vNj)
    mZ        = as.data.frame(data[,Z])
    names(mZ) = Z
    if(FALSE%in%unlist(lapply(as.data.frame(mZ),is.numeric))){
      
      mZkeep            = as.matrix(mZ[,unlist(lapply(as.data.frame(mZ),is.numeric))])
      colnames(mZkeep)  = colnames(mZ)[unlist(lapply(as.data.frame(mZ),is.numeric))]
      mZclean           = as.matrix(mZ[,!unlist(lapply(as.data.frame(mZ),is.numeric))])
      mZclean_names     = colnames(mZ)[!unlist(lapply(as.data.frame(mZ),is.numeric))]
      if(ncol(mZclean) == 1){
        
        mZclean = as.matrix(as.character(mZclean))
        
      } else{
        
        mZclean = apply(mZclean,2,function(x){as.character(x)})
        
      }
      mZclean[mZclean==""] = NA
      for(i in 1:ncol(mZclean)){
        
        for(j in 2:length(unique(na.omit(mZclean[,i])))){
          
          lvl               = sort(unique(na.omit(mZclean[,i])))[j]
          cleaned           = as.matrix(ifelse(mZclean[,i] == lvl,1,0))
          colnames(cleaned) = paste0(mZclean_names[i],".",lvl)
          mZkeep            = cbind(mZkeep, cleaned)
          
        }
        
      }
      mZ = as.matrix(apply(mZkeep,2,function(x){as.numeric(x)}))
      
      
    }
    mZ      = cbind(1,mZ)
    mZ      = as.matrix(mZ)
    Z       = c("Intercept",colnames(mZ)[-1])
    P       = ncol(mZ)
    mZh        = as.data.frame(data[,Zh])
    names(mZh) = Zh
    if(FALSE%in%unlist(lapply(as.data.frame(mZh),is.numeric))){
      
      mZhkeep            = as.matrix(mZh[,unlist(lapply(as.data.frame(mZh),is.numeric))])
      colnames(mZhkeep)  = colnames(mZh)[unlist(lapply(as.data.frame(mZh),is.numeric))]
      mZhclean           = as.matrix(mZh[,!unlist(lapply(as.data.frame(mZh),is.numeric))])
      mZhclean_names     = colnames(mZh)[!unlist(lapply(as.data.frame(mZh),is.numeric))]
      if(ncol(mZhclean) == 1){
        
        mZhclean = as.matrix(as.character(mZhclean))
        
      } else{
        
        mZhclean = apply(mZhclean,2,function(x){as.character(x)})
        
      }
      mZhclean[mZhclean==""] = NA
      for(i in 1:ncol(mZhclean)){
        
        for(j in 2:length(unique(na.omit(mZhclean[,i])))){
          
          lvl               = sort(unique(na.omit(mZhclean[,i])))[j]
          cleaned           = as.matrix(ifelse(mZhclean[,i] == lvl,1,0))
          colnames(cleaned) = paste0(mZhclean_names[i],".",lvl)
          mZhkeep            = cbind(mZhkeep, cleaned)
          
        }
        
      }
      mZh = as.matrix(apply(mZhkeep,2,function(x){as.numeric(x)}))
      
      
    }
    mZh     = cbind(1,mZh)
    mZh     = as.matrix(mZh)
    Zh      = c("Intercept",colnames(mZh)[-1])
    P_high  = ncol(mZh)
    
    # Initialization
    mY_df             = data.frame(mY,id_high)
    mY_aggr           = as.matrix(mY_df%>%group_by_all%>%count)
    iHf               = dim(mY_aggr)[2]
    mY_aggr           = mY_aggr[order(mY_aggr[,iHf-1]),]
    freq              = mY_aggr[,iHf]
    mY_unique         = mY_aggr[,-c(iHf-1,iHf)]
    if(kmea==FALSE){
      
      clusfoo     = klaR::kmodes(mY_unique,modes=iT)$cluster
      
    } else{
      
      prscores    = prcomp(mY_unique)
      num_pc      = max(round(ncol(mY_unique)/2),(sum(cumsum(prscores$sdev^2/sum(prscores$sdev^2))<0.85)+1))
      prscores    = prscores$x[,1:num_pc]
      spectclust  = kmeans(prscores,centers=iT,iter.max=100,nstart=100)
      clusfoo     = spectclust$cluster
      
    }
    mU                = vecTomatClass(clusfoo)
    outinit1          = LCA_fast(mY_unique,freq,iT,mU,maxIter,tol,reord)
    mU_full           = vecTomatClass(rep(outinit1$vModalAssnm+1,outinit1$freq))
    vUmodal           = rep(0,iJ)
    mUmodal           = matrix(0,iJ,iT)
    foo               = 0
    for(j in 1:iJ){
      
      if(vNj[j]>1){
        
        mUmodal[j,] = colSums(mU_full[foo+(1:vNj[j]),])/vNj[j]
        foo         = foo+vNj[j]
        
      } else{
        
        mUmodal[j,] = colSums(t(mU_full[foo+1,]))/vNj[j]
        foo         = foo+vNj[j]
        
      }
      
    }
    high_out          = kmeans(mUmodal,centers=iM,iter.max=100,nstart=100)
    vOmega_start      = high_out$size/iJ
    mIndex_indiv_high = vecTomatClass(rep(high_out$cluster,vNj))
    highreord         = order(vOmega_start,decreasing=TRUE)
    vOmega_start      = vOmega_start[highreord]
    index_indiv_high  = apply(mIndex_indiv_high[,highreord],1,which.max)
    mPhi_start        = outinit1$mPhi
    mPi_fast          = table(index_indiv_high,rep(outinit1$vModalAssnm+1,outinit1$freq))
    mPi_start         = t(mPi_fast/rowSums(mPi_fast))
    outinit2          = MLTLCA(mY,vNj,vOmega_start,mPi_start,mPhi_start,maxIter,tol,reord)
    cGamma_start      = array(c(rbind(outinit2$mGamma,matrix(0,(iT-1)*(P-1),iM))),c(iT-1,P,iM))
    mPhi_start        = outinit2$mPhi
    mStep1Var         = outinit2$Varmat
    mDelta_start      = matrix(0,iM-1,P_high)
    mDelta_start[,1]  = outinit2$vDelta
    
    # Remove missing values on mZ, update mY and vNj
    nomissing = complete.cases(mZ)&complete.cases(mZh)
    mY        = mY[nomissing,]
    idnames   = idnames[sort(unique(id_high[nomissing]))]
    id_high   = id_high[nomissing]
    vNj       = as.vector(table(id_high))
    mZ        = mZ[nomissing,]
    mZh       = mZh[nomissing,]
    mZh       = mZh[!duplicated(id_high),]
    
    # Estimation
    est = MLTLCA_covlowhigh(mY,mZ,mZh,vNj,mDelta_start,cGamma_start,mPhi_start,mStep1Var,maxIter,tol,fixedpars,NRtol,NRmaxit)
    
    # Return
    if(extout == FALSE){
      
      # Create empty list for output
      out = list()
      
      # Add matrix of sample means of high-level class proportions
      out$vOmega_avg            = as.matrix(apply(est$mOmega,2,mean))
      rownames(out$vOmega_avg)  = paste0("P(G",1:iM,")")
      colnames(out$vOmega_avg)  = ""
      
      # Add matrix of sample means of conditional low-level class proportions
      out$mPi_avg           = apply(est$cPi,3,function(x){apply(x,2,mean)})
      rownames(out$mPi_avg) = paste0("P(C",1:iT,"|G)")
      colnames(out$mPi_avg) = paste0("G",1:iM)
      
      # Add matrix of conditional response probabilities
      out$mPhi            = est$mPhi
      rownames(out$mPhi)  = paste0("P(",Y,"|C)")
      colnames(out$mPhi)  = paste0("C",1:iT)
      
      # Add deltas
      out$mDelta            = t(est$mDelta)
      rownames(out$mDelta)  = paste0("delta(",Zh,"|G)")
      colnames(out$mDelta)  = paste0("G",2:iM)
      
      # Add gammas
      out$cGamma = array(apply(est$cGamma,3,t),c(P,iT-1,iM),dimnames=list(paste0("gamma(",Z,"|C)"),paste0("C",2:iT,"|G"),colnames(out$mPi_avg)))
      
      # Add corrected standard errors for delta
      out$SEs_cor_delta = t(matrix(est$SEs_cor[1:(P_high*(iM-1))],iM-1,P_high))
      rownames(out$SEs_cor_delta) = rownames(out$mDelta)
      colnames(out$SEs_cor_delta) = colnames(out$mDelta)
      
      # Add corrected standard errors for gamma
      out$SEs_cor_gamma = array(apply(array(est$SEs_cor[(1+P_high*(iM-1)):(P_high*(iM-1)+P*(iT-1)*iM)],c(iT-1,P,iM)),3,t),dim(out$cGamma),dimnames=dimnames(out$cGamma))
      
      # Add Akaike information criterion
      out$AIC = est$AIC
      
      # Add low-level Bayesian information criterion
      out$BIClow = est$BIClow
      
      # Add high-level Bayesian information criterion
      out$BIChigh = est$BIChigh
      
      # Add low-level integrated completed likelihood Bayesian information criterion
      out$ICL_BIClow = est$ICL_BIClow
      
      # Add high-level integrated completed likelihood Bayesian information criterion
      out$ICL_BIChigh = est$ICL_BIChigh
      
      # Add low-level entropy R-sqr
      out$R2entr_low = est$R2entr_low
      
      # Add high-level entropy R-sqr
      out$R2entr_high = est$R2entr_high
      
      # Add number of iterations
      out$iter = est$iter
      
      # Add log-likelihood series
      out$LLKSeries = est$LLKSeries
      
      # Add specification
      out$spec            = as.matrix("Multilevel structural LC model with low- and high-level covariates")
      rownames(out$spec)  = colnames(out$spec) = ""
      
    } else{
      
      # Create empty list for output
      out = list()
      
      # Add matrix of high-level class proportions
      out$mOmega            = est$mOmega
      colnames(out$mOmega)  = paste0("P(G",1:iM,")")
      
      # Add matrix of sample means of high-level class proportions
      out$vOmega_avg            = as.matrix(apply(est$mOmega,2,mean))
      rownames(out$vOmega_avg)  = colnames(out$mOmega)
      colnames(out$vOmega_avg)  = ""
      
      # Add matrix of conditional low-level class proportions
      out$mPi           = matrix(est$cPi,nrow(est$cPi),iT*iM)
      colnames(out$mPi) = paste0(rep(paste0("P(C",1:iT,"|"),iM),rep(paste0("G",1:iM,")"),rep(iT,iM)))
      
      # Add matrix of sample means of conditional low-level class proportions
      out$mPi_avg           = apply(est$cPi,3,function(x){apply(x,2,mean)})
      rownames(out$mPi_avg) = paste0("P(C",1:iT,"|G)")
      colnames(out$mPi_avg) = paste0("G",1:iM)
      
      # Add matrix of conditional response probabilities
      out$mPhi            = est$mPhi
      rownames(out$mPhi)  = paste0("P(",Y,"|C)")
      colnames(out$mPhi)  = paste0("C",1:iT)
      
      # Add deltas
      out$mDelta            = t(est$mDelta)
      rownames(out$mDelta)  = paste0("delta(",Zh,"|G)")
      colnames(out$mDelta)  = paste0("G",2:iM)
      
      # Add gammas
      out$cGamma = array(apply(est$cGamma,3,t),c(P,iT-1,iM),dimnames=list(paste0("gamma(",Z,"|C)"),paste0("C",2:iT,"|G"),colnames(out$mPi_avg)))
      
      # Add betas
      out$mBeta           = est$mBeta
      rownames(out$mBeta) = paste0("beta(",Y,"|C)")
      colnames(out$mBeta) = colnames(out$mPhi)
      
      # Add vector of model parameters
      out$parvec            = est$parvec
      rownames(out$parvec)  = c(paste0(rep(substr(rownames(out$mDelta),1,nchar(rownames(out$mDelta))-2),rep(iM-1,P_high)),rep(colnames(out$mDelta),P_high),")"),paste0(apply(out$cGamma,3,function(x){paste0(rep(substr(rownames(x),1,nchar(rownames(x))-2),rep(iT-1,P)),rep(substr(colnames(x),1,nchar(colnames(x))-2),P),"|")}),rep(unlist(dimnames(out$cGamma)[3]),rep(P*(iT-1),iM)),")"),paste0(rep(substr(rownames(out$mBeta),1,nchar(rownames(out$mBeta))-2),iT),rep(colnames(out$mBeta),rep(iH,iT)),")"))
      colnames(out$parvec)  = ""
      
      # Add Fisher information matrix
      out$Infomat           = est$Infomat
      rownames(out$Infomat) = colnames(out$Infomat) = rownames(out$parvec)
      
      # Add cube of Fisher information matrix for delta
      out$cDelta_Info = array(est$cDelta_Info, dim(est$cDelta_Info), dimnames=list(paste0("delta(",Zh,"|G)"),paste0("delta(",Zh,"|G)"),colnames(out$mDelta)))
      
      # Add cube of Fisher information matrix for gamma
      out$cGamma_Info = array(est$cGamma_Info, dim(est$cGamma_Info), dimnames=list(paste0("gamma(",Z,"|C|G)"),paste0("gamma(",Z,"|C|G)"),paste0("C",2:iT,rep(paste0("|G",1:iM),rep(iT-1,iM)))))
      
      # Add inverse of the information matrix from the second step
      out$mV2           = est$mV2
      rownames(out$mV2) = colnames(out$mV2) = c(paste0(rep(substr(rownames(out$mDelta),1,nchar(rownames(out$mDelta))-2),rep(iM-1,P_high)),rep(colnames(out$mDelta),P_high),")"),paste0(apply(out$cGamma,3,function(x){paste0(rep(substr(rownames(x),1,nchar(rownames(x))-2),rep(iT-1,P)),rep(substr(colnames(x),1,nchar(colnames(x))-2),P),"|")}),rep(unlist(dimnames(out$cGamma)[3]),rep(P*(iT-1),iM)),")"))
      
      # Add mQ
      out$mQ            = est$mQ
      rownames(out$mQ)  = colnames(out$mV2) = rownames(out$mV2)
      
      # Add uncorrected variance-covariance matrix
      out$Varmat_unc            = est$Varmat[1:(P_high*(iM-1)+P*(iT-1)*iM),1:(P_high*(iM-1)+P*(iT-1)*iM)]
      rownames(out$Varmat_unc)  = colnames(out$Varmat_unc) = rownames(out$mV2)
      
      # Add corrected variance-covariance matrix
      out$Varmat_cor            = est$mVar_corr
      rownames(out$Varmat_cor)  = colnames(out$Varmat_cor) = rownames(out$Varmat_unc)
      
      # Add vector of uncorrected standard errors
      out$SEs_unc           = as.matrix(est$SEs_unc[1:(P_high*(iM-1)+P*(iT-1)*iM),])
      rownames(out$SEs_unc) = rownames(out$mV2)
      colnames(out$SEs_unc) = ""
      
      # Add vector of corrected standard errors
      out$SEs_cor           = as.matrix(est$SEs_cor[1:(P_high*(iM-1)+P*(iT-1)*iM),])
      rownames(out$SEs_cor) = rownames(out$mV2)
      colnames(out$SEs_cor) = ""
      
      # Add corrected standard errors for delta
      out$SEs_cor_delta = t(matrix(est$SEs_cor[1:(P_high*(iM-1))],iM-1,P_high))
      rownames(out$SEs_cor_delta) = rownames(out$mDelta)
      colnames(out$SEs_cor_delta) = colnames(out$mDelta)
      
      # Add corrected standard errors for gamma
      out$SEs_cor_gamma = array(apply(array(est$SEs_cor[(1+P_high*(iM-1)):(P_high*(iM-1)+P*(iT-1)*iM)],c(iT-1,P,iM)),3,t),dim(out$cGamma),dimnames=dimnames(out$cGamma))
      
      # Add epsilon
      out$eps = est$eps
      
      # Add cube of joint posterior class membership probabilities
      out$cPMX = array(est$cPMX,dim(est$cPMX),dimnames=list(NULL,paste0("C",1:iT,",G"), paste0("G", 1:iM)))
      
      # Add cube of log of joint posterior class membership probabilities
      out$cLogPMX = array(est$cLogPMX,dim(est$cLogPMX),dimnames=dimnames(out$cPMX))
      
      # Add cube of conditional posterior low-level class membership probabilities
      out$cPX = array(est$cPX,dim(est$cPX),dimnames=list(NULL,paste0("C", 1:iT,"|G"),paste0("G", 1:iM)))
      
      # Add cube of log of conditional posterior low-level class membership probabilities
      out$cLogPX = array(est$cLogPX,dim(est$cLogPX),dimnames=dimnames(out$cPX))
      
      # Add matrix of posterior high-level class membership probabilities for low-level units after marginalizing over low-level classes
      out$mSumPX            = est$mSumPX
      colnames(out$mSumPX)  = colnames(out$mPi_avg)
      
      # Add matrix of posterior high-level class membership probabilities for high-level units
      out$mPW           = est$mPW
      rownames(out$mPW) = idnames
      colnames(out$mPW) = colnames(out$mPi_avg)
      
      # Add matrix of log of posterior high-level class membership probabilities for high-level units
      out$mlogPW            = est$mlogPW
      rownames(out$mlogPW)  = idnames
      colnames(out$mlogPW)  = colnames(out$mPi_avg)
      
      # Add matrix of posterior high-level class membership probabilities for low-level units
      out$mPW_N           = est$mPW_N
      colnames(out$mPW_N) = colnames(out$mPi_avg)
      
      # Add matrix of posterior low-level class membership probabilities for low-level units after marginalizing over high-level classes
      out$mPMsumX           = est$mPMsumX
      colnames(out$mPMsumX) = colnames(out$mPhi)
      
      # Add low-level entropy R-sqr
      out$R2entr_low = est$R2entr_low
      
      # Add high-level entropy R-sqr
      out$R2entr_high = est$R2entr_high
      
      # Add Akaike information criterion
      out$AIC = est$AIC
      
      # Add low-level Bayesian information criterion
      out$BIClow = est$BIClow
      
      # Add high-level Bayesian information criterion
      out$BIChigh = est$BIChigh
      
      # Add low-level integrated completed likelihood Bayesian information criterion
      out$ICL_BIClow = est$ICL_BIClow
      
      # Add high-level integrated completed likelihood Bayesian information criterion
      out$ICL_BIChigh = est$ICL_BIChigh
      
      # Add number of iterations
      out$iter = est$iter
      
      # Add log-likelihood series
      out$LLKSeries = est$LLKSeries
      
      # Add current log-likelihood for high-level units
      out$vLLK            = est$vLLK
      rownames(out$vLLK)  = idnames
      colnames(out$vLLK)  = ""
      
      # Add matrix of model parameter contributions to log-likelihood score
      out$mScore            = est$mScore
      colnames(out$mScore)  = rownames(out$parvec)
      
      # Add subset of matrix of model parameter contributions to log-likelihood score for delta
      out$mDelta_Score            = est$mDelta_Score
      rownames(out$mDelta_Score)  = idnames
      colnames(out$mDelta_Score)  = paste0(rep(substr(rownames(out$mDelta),1,nchar(rownames(out$mDelta))-2),rep(iM-1,P_high)),rep(colnames(out$mDelta),P_high),")")
      
      # Add subset of matrix of model parameter contributions to log-likelihood score for gamma
      out$mGamma_Score            = est$mGamma_Score
      colnames(out$mGamma_Score)  = paste0(apply(out$cGamma,3,function(x){paste0(rep(substr(rownames(x),1,nchar(rownames(x))-2),rep(iT-1,P)),rep(substr(colnames(x),1,nchar(colnames(x))-2),P),"|")}),rep(unlist(dimnames(out$cGamma)[3]),rep(P*(iT-1),iM)),")")
      
      # Add specification
      out$spec            = as.matrix("Multilevel structural LC model with low- and high-level covariates")
      rownames(out$spec)  = colnames(out$spec) = ""
      
    }
    if(dataout == TRUE){
      
      out$data_low            = cbind(mY,id_high,mZ)
      rownames(out$data_low)  = NULL
      colnames(out$data_low)  = c(Y,colnames(data)[!colnames(data)%in%c(Y,Z,Zh)],Z)
      out$data_high           = mZh
      rownames(out$data_high) = idnames
      colnames(out$data_high) = Zh
      
    }
    return(out)
    
  }
  
}

multiLCA = function(data, Y, iT, id_high = NULL, iM = NULL, Z = NULL, Zh = NULL,
                    extout = FALSE, dataout = FALSE, kmea = TRUE,
                    sequential = TRUE, numFreeCores = 2,
                    maxIter = 1e3, tol = 1e-8, fixed = TRUE, reord = TRUE,
                    fixedpars = 1, NRtol = 1e-6, NRmaxit = 100, verbose = TRUE){
  
  #########################
  ### 1. Control inputs ###
  #########################
  
  # <data> is matrix or dataframe
  if(!is.data.frame(data)){
    
    stop("data must be matrix or dataframe.",call.=FALSE)
    
  }
  
  # <Y>, <id_high>, <Z> and <Zh> are column names in <data>
  if(FALSE%in%(c(Y,id_high,Z,Zh)%in%colnames(data))){
    
    if(FALSE%in%(Y%in%colnames(data))){
      
      stop("Not all Y in data.",call.=FALSE)
      
    }
    if(!is.null(id_high)&FALSE%in%(id_high%in%colnames(data))){
      
      stop("id_high not in data.",call.=FALSE)
      
    }
    if(!is.null(Z)&FALSE%in%(Z%in%colnames(data))){
      
      stop("Not all Z in data.",call.=FALSE)
      
    }
    if(!is.null(Zh)&FALSE%in%(Zh%in%colnames(data))){
      
      stop("Not all Zh in data.",call.=FALSE)
      
    }
    
  }
  
  # Single <id_high>
  if(!is.null(id_high)){
    
    if(length(id_high)!=1){
      
      stop("Invalid id_high.",call.=FALSE)
      
    }
    
  }
  
  # <iT> positive integer or sequential positive integers
  if(!length(iT)>=1){
    
    stop("Invalid iT.",call.=FALSE)
    
  } else{
    
    if(length(iT)==1){
      
      if(!is.numeric(iT)){
        
        stop("Invalid iT.",call.=FALSE)
        
      } else if(iT!=round(iT)){
        
        stop("Invalid iT.",call.=FALSE)
        
      } else if(iT<1){
        
        stop("Invalid iT.",call.=FALSE)
        
      }
      
    } else{
      
      if(!is.vector(iT)|is.list(iT)){
        
        stop("Invalid iT.",call.=FALSE)
        
      } else if(!is.numeric(iT)){
        
        stop("Invalid iT.",call.=FALSE)
        
      } else if(FALSE%in%(iT==round(iT))){
        
        stop("Invalid iT.",call.=FALSE)
        
      } else if(min(iT)<1){
        
        stop("Invalid iT.",call.=FALSE)
        
      } else if(length(iT)!=length(min(iT):max(iT))){
        
        stop("Invalid iT.",call.=FALSE)
        
      } else if(FALSE%in%(iT==min(iT):max(iT))){
        
        stop("Invalid iT.",call.=FALSE)
        
      }
      
    }
    
  }
  
  # <iM> positive integer or sequential positive integers
  if(!is.null(iM)&!length(iM)>=1){
    
    stop("Invalid iM.",call.=FALSE)
    
  } else if(!is.null(iM)){
    
    if(length(iM)==1){
      
      if(!is.numeric(iM)){
        
        stop("Invalid iM.",call.=FALSE)
        
      } else if(iM!=round(iM)){
        
        stop("Invalid iM.",call.=FALSE)
        
      } else if(iM<1){
        
        stop("Invalid iM.",call.=FALSE)
        
      }
      
    } else{
      
      if(!is.vector(iM)|is.list(iM)){
        
        stop("Invalid iM.",call.=FALSE)
        
      } else if(!is.numeric(iM)){
        
        stop("Invalid iM.",call.=FALSE)
        
      } else if(FALSE%in%(iM==round(iM))){
        
        stop("Invalid iM.",call.=FALSE)
        
      } else if(min(iM)<1){
        
        stop("Invalid iM.",call.=FALSE)
        
      } else if(length(iM)!=length(min(iM):max(iM))){
        
        stop("Invalid iM.",call.=FALSE)
        
      } else if(FALSE%in%(iM==min(iM):max(iM))){
        
        stop("Invalid iM.",call.=FALSE)
        
      }
      
    }
    
  }
  
  ##############################################
  ### 2. Identify specification and approach ###
  ##############################################
  
  # Identify specification
  if(is.null(id_high)&is.null(iM)){
    
    specification = "single-level"
    
  } else if(!is.null(id_high)&!is.null(iM)){
    
    specification = "multilevel"
    
  } else{
    
    if(is.null(id_high)){
      
      stop("iM specified, id_high not.",call.=FALSE)
      
    }
    if(is.null(iM)){
      
      stop("id_high specified, iM not.",call.=FALSE)
      
    }
    
  }
  if(is.null(Z)){
    
    specification = paste(specification,"measurement")
    
    if(!is.null(Zh)){
      
      stop("Zh specified, Z not.",call.=FALSE)
      
    }
    
  } else if(!is.null(Z)){
    
    specification = paste(specification,"structural")
    
    if(specification=="single-level structural"&!is.null(Zh)){
      
      stop("Zh specified, iM and id_high not.",call.=FALSE)
      
    }
    
  }
  
  # Identify approach
  if(length(iT)==1&is.null(iM)){
    
    approach = "direct"
    
  } else if(length(iT)==1&length(iM)==1){
    
    approach = "direct"
    
  } else if(length(iT)>1&is.null(iM)){
    
    approach = "model selection on low"
    
    if("structural"%in%unlist(strsplit(specification," "))){
      
      warning("Estimating measurement model in model selection.",call.=FALSE)
      
    }
    
  } else if(length(iT)>1&length(iM)==1){
    
    approach = "model selection on low"
    
    if("structural"%in%unlist(strsplit(specification," "))){
      
      warning("Estimating measurement model in model selection.",call.=FALSE)
      
    }
    
  } else if(length(iT)==1&length(iM)>1){
    
    approach = "model selection on high"
    
    if("structural"%in%unlist(strsplit(specification," "))){
      
      warning("Estimating measurement model in model selection.",call.=FALSE)
      
    }
    
  } else if(length(iT)>1&length(iM)>1){
    
    approach = "model selection on low and high"
    
    if("structural"%in%unlist(strsplit(specification," "))){
      
      warning("Estimating measurement model in model selection.",call.=FALSE)
      
    }
    
  }
  
  #####################
  ### 3. Clean data ###
  #####################
  
  # Remove irrelevant columns
  data = data[,c(Y,id_high,Z,Zh)]
  
  # Remove missing values with respect to <Y> and <id_high>
  if(nrow(data)>nrow(na.omit(data[,c(Y,id_high)]))){
    
    warning("Missing values removed.",call.=FALSE)
    
    data = data[complete.cases(data[,Y]),]
    
  }
  
  # Replace <fixed> and <reord> with numeric
  if(fixed==TRUE){
    
    fixed = 1
    
  } else{
    
    fixed = 0
    
  }
  if(reord==TRUE){
    
    reord = 1
    
  } else{
    
    reord = 0
    
  }
  
  #####################################
  ### 4. Control contents of <data> ###
  #####################################
  
  # <Y> 0-1 coded
  if(TRUE%in%(apply(data[,Y],2,function(x){length(unique(x))})>2)){
    
    stop("Items not 0-1 coded.",call.=FALSE)
    
  } else if(FALSE%in%apply(data[,Y],2,function(x){sort(unique(x))==c(0,1)})){
    
    stop("Items not 0-1 coded.",call.=FALSE)
    
  }
  
  # No missings in <id_high>
  if(TRUE%in%is.na(data[,id_high])){
    
    stop("id_high cannot contain missing values.",call.=FALSE)
    
  }
  
  # <id_high> not labelled "Intercept"
  if(!is.null(id_high)){
    
    if(id_high=="Intercept"){
      
      stop("id_high may not be labelled 'Intercept'.",call.=FALSE)
      
    }
    
  }
  
  # No duplicate <Y>, <Z> or >Zh>
  if(TRUE%in%duplicated(as.list(data[,Y]))){
    
    stop("Duplicate items.",call.=FALSE)
    
  }
  if(!is.null(Z)){
    
    if(length(Z)>2){
      
      if(TRUE%in%duplicated(as.list(data[,Z]))){
        
        stop("Duplicate low-level covariates.",call.=FALSE)
        
      }
      
    }
    
  }
  if(!is.null(Zh)){
    
    if(length(Zh)>2){
      
      if(TRUE%in%duplicated(as.list(data[,Zh]))){
        
        stop("Duplicate high-level covariates.",call.=FALSE)
        
      }
      
    }
    
  }
  
  # No constant <Y>, <Z> or <Zh>
  if(FALSE%in%apply(data[,Y],2,function(x){length(unique(x))>1})){
    
    stop("Constant in items.",call.=FALSE)
    
  }
  if(!is.null(Z)){
    
    if(FALSE%in%apply(data[,Y],2,function(x){length(unique(x))>1})){
      
      stop("Constant in low-level covariates.",call.=FALSE)
      
    }
    
  }
  if(!is.null(Zh)){
    
    if(FALSE%in%apply(data[,Y],2,function(x){length(unique(x))>1})){
      
      stop("Constant in high-level covariates.",call.=FALSE)
      
    }
    
  }
  
  #############################
  ### 5. Perform estimation ###
  #############################
  
  # Direct to estimation without model selection
  if(approach=="direct"){
    
    # Non-mixture specification: one-class single-level model
    if(("single-level"%in%unlist(strsplit(specification," ")))&iT==1){
      
      iH = ncol(data[,Y])
      iN = nrow(data)
      
      ll = sum(apply(data[,Y],2,function(x){dbinom(x,1,mean(x),log=TRUE)}))
      
      AIC     = -2*ll+2*iH
      BIC     = -2*ll+iH*log(iN)
      avgCE   = 0
      R2entr  = 1
      
      if(verbose)cat("\nCALL:\n")
      if(verbose)print(match.call())
      if(verbose)cat("\nSPECIFICATION:\n")
      if(verbose)print("iT = 1")
      if(verbose)cat("\n------------------------------------\n")
      if(verbose)cat("\nMODEL AND CLASSIFICATION STATISTICS:\n")
      stat            = t(round(c(ll,AIC,BIC,avgCE,R2entr),2))
      rownames(stat)  = ""
      colnames(stat)  = c("LL","AIC","BIC","ClassErr","EntR-sqr")
      if(verbose)print(stat,right=TRUE)
      if(verbose)cat("\n")
      
      stop("Non-mixture specification.",call.=FALSE)
      
    }
    
    # Non-mixture specification: one-low one-high multilevel model
    if(("multilevel"%in%unlist(strsplit(specification," ")))){
      
      if(iT==1&iM==1){
        
        iH = ncol(data[,Y])
        iN = nrow(data)
        iJ = length(table(data[,id_high])[table(data[,id_high]) > 0])
        
        ll = sum(apply(data[,Y],2,function(x){dbinom(x,1,mean(x),log=TRUE)}))
        
        AIC         = -2*ll+2*iH
        BIClow      = -2*ll+iH*log(iN)
        BIChigh     = -2*ll+iH*log(iJ)
        ICL_BIClow  = NA
        ICL_BIChigh = NA
        
        if(verbose)cat("\nCALL:\n")
        if(verbose)print(match.call())
        if(verbose)cat("\nSPECIFICATION:\n")
        if(verbose)print("iT = 1, iM = 1")
        if(verbose)cat("\n---------------------------\n")
        if(verbose)cat("\nGOODNESS-OF-FIT STATISTICS:\n")
        stat              = t(round(c(ll,AIC,BIClow,BIChigh,ICL_BIClow,ICL_BIChigh),2))
        stat[is.na(stat)] = "-"
        rownames(stat)    = ""
        colnames(stat)    = c("LL","AIC","BIClow","BIChigh","ICL_BIClow","ICL_BIChigh")
        if(verbose)print(noquote(stat),right=TRUE)
        if(verbose)cat("\n")
        
        stop("Non-mixture specification.",call.=FALSE)
        
      }
      
    }
    
    # Non-mixture specification: one-low multilevel model
    if(("multilevel"%in%unlist(strsplit(specification," ")))){
      
      if(iT==1&iM>1){
        
        if(verbose)cat("\nCALL:\n")
        if(verbose)print(match.call())
        if(verbose)cat("\nSPECIFICATION:\n")
        if(verbose)print("iT = 1, iM > 1")
        if(verbose)cat("\n---------------------------\n")
        if(verbose)cat("\nGOODNESS-OF-FIT STATISTICS:\n")
        stat              = t(rep("-",6))
        rownames(stat)    = ""
        colnames(stat)    = c("LL","AIC","BIClow","BIChigh","ICL_BIClow","ICL_BIChigh")
        if(verbose)print(noquote(stat),right=TRUE)
        if(verbose)cat("\n")
        
        stop("Non-mixture specification.",call.=FALSE)
        
      }
      
    }
    
    # Non-mixture specification: one-high multilevel model
    if(("multilevel"%in%unlist(strsplit(specification," ")))){
      
      if(iT>1&iM==1){
        
        warning("Non-mixture specification. Estimating single-level model.",call.=FALSE)
        
        id_high = NULL
        iM      = NULL
        
        if(!is.null(Zh)){
          
          Z = c(Z,Zh)
          
          Zh = NULL
          
        }
        
      }
      
    }
    
    # Mixture specification
    est       = estlca(data=data,Y=Y,iT=iT,id_high=id_high,iM=iM,Z=Z,Zh=Zh,extout=extout,dataout=dataout,kmea=kmea,maxIter=maxIter,tol=tol,fixed=fixed,reord=reord,fixedpars=fixedpars,NRtol=NRtol,NRmaxit=NRmaxit)
    est$call  = match.call()
    
  }
  
  # Model selection with respect to <iT>
  if(approach=="model selection on low"){
    
    # Single-level model
    if("single-level"%in%unlist(strsplit(specification," "))){
      
      iT_max = max(iT)
      iT_min = min(iT)
      num_iT = iT_max-iT_min+1
      
      ll_seq      = rep(NA,num_iT)
      AIC_seq     = rep(NA,num_iT)
      BIC_seq     = rep(NA,num_iT)
      avgCE_seq   = rep(NA,num_iT)
      R2entr_seq  = rep(NA,num_iT)
      
      LCAseq = list()
      
      for(i in iT_min:iT_max){
        
        if(i==1){
          
          LCAseq[[1]] = list()
          
          iH = ncol(data[,Y])
          iN = nrow(data)
          
          ll = sum(apply(data[,Y],2,function(x){dbinom(x,1,mean(x),log=TRUE)}))
          
          ll_seq[1]     = ll
          AIC_seq[1]    = -2*ll+2*iH
          BIC_seq[1]    = -2*ll+iH*log(iN)
          avgCE_seq[1]  = 0
          R2entr_seq[1] = 1
          
        } else{
          
          LCAseq[[1+i-iT_min]] = estlca(data=data,Y=Y,iT=i,id_high=NULL,iM=NULL,
                                        Z=NULL,Zh=NULL,
                                        extout=FALSE,dataout=FALSE,kmea=TRUE,
                                        maxIter=1e3,tol=1e-8,fixed=TRUE,
                                        reord=TRUE,fixedpars=1,
                                        NRtol=1e-6,NRmaxit=100)
          ll_seq[1+i-iT_min]      = tail(LCAseq[[1+i-iT_min]]$LLKSeries,1)
          AIC_seq[1+i-iT_min]     = LCAseq[[1+i-iT_min]]$AIC
          BIC_seq[1+i-iT_min]     = LCAseq[[1+i-iT_min]]$BIC
          avgCE_seq[1+i-iT_min]   = LCAseq[[1+i-iT_min]]$AvgClassErrProb
          R2entr_seq[1+i-iT_min]  = LCAseq[[1+i-iT_min]]$R2entr
          
        }
        
      }
      
      iT_best   = which.min(BIC_seq)+iT_min-1
      est       = LCAseq[[1+iT_best-iT_min]]
      est$call  = match.call()
      
      if(verbose)cat("\nCALL:\n")
      if(verbose)print(match.call())
      if(verbose)cat("\n------------------------------------\n")
      if(verbose)cat("\nMODEL AND CLASSIFICATION STATISTICS:\n")
      stat            = round(cbind(iT_min:iT_max,ll_seq,AIC_seq,BIC_seq,avgCE_seq,R2entr_seq),2)
      rownames(stat)  = rep("",num_iT)
      colnames(stat)  = c("iT","LL","AIC","BIC","ClassErr","EntR-sqr")
      if(verbose)print(stat,right=TRUE)
      if(verbose)cat("\nBest iT based on BIC:",iT_best,"\n")
      
    }
    
    # Multilevel model
    if("multilevel"%in%unlist(strsplit(specification," "))){
      
      # Non-mixture specification: one-high multilevel model
      if(iM==1){
        
        warning("Non-mixture specification. Estimating single-level models.",call.=FALSE)
        
        iT_max = max(iT)
        iT_min = min(iT)
        num_iT = iT_max-iT_min+1
        
        ll_seq      = rep(NA,num_iT)
        AIC_seq     = rep(NA,num_iT)
        BIC_seq     = rep(NA,num_iT)
        avgCE_seq   = rep(NA,num_iT)
        R2entr_seq  = rep(NA,num_iT)
        
        LCAseq = list()
        
        for(i in iT_min:iT_max){
          
          if(i==1){
            
            LCAseq[[1]] = list()
            
            iH = ncol(data[,Y])
            iN = nrow(data)
            
            ll = sum(apply(data[,Y],2,function(x){dbinom(x,1,mean(x),log=TRUE)}))
            
            ll_seq[1]     = ll
            AIC_seq[1]    = -2*ll+2*iH
            BIC_seq[1]    = -2*ll+iH*log(iN)
            avgCE_seq[1]  = 0
            R2entr_seq[1] = 1
            
          } else{
            
            LCAseq[[1+i-iT_min]] = estlca(data=data,Y=Y,iT=i,id_high=NULL,iM=NULL,
                                          Z=NULL,Zh=NULL,
                                          extout=FALSE,dataout=FALSE,kmea=TRUE,
                                          maxIter=1e3,tol=1e-8,fixed=TRUE,
                                          reord=TRUE,fixedpars=1,
                                          NRtol=1e-6,NRmaxit=100)
            ll_seq[1+i-iT_min]      = tail(LCAseq[[1+i-iT_min]]$LLKSeries,1)
            AIC_seq[1+i-iT_min]     = LCAseq[[1+i-iT_min]]$AIC
            BIC_seq[1+i-iT_min]     = LCAseq[[1+i-iT_min]]$BIC
            avgCE_seq[1+i-iT_min]   = LCAseq[[1+i-iT_min]]$AvgClassErrProb
            R2entr_seq[1+i-iT_min]  = LCAseq[[1+i-iT_min]]$R2entr
            
          }
          
        }
        
        iT_best   = which.min(BIC_seq)+iT_min-1
        est       = LCAseq[[1+iT_best-iT_min]]
        est$call  = match.call()
        
        if(verbose)cat("\nCALL:\n")
        if(verbose)print(match.call())
        if(verbose)cat("\n------------------------------------\n")
        if(verbose)cat("\nMODEL AND CLASSIFICATION STATISTICS:\n")
        stat            = round(cbind(iT_min:iT_max,ll_seq,AIC_seq,BIC_seq,avgCE_seq,R2entr_seq),2)
        rownames(stat)  = rep("",num_iT)
        colnames(stat)  = c("iT","LL","AIC","BIC","ClassErr","EntR-sqr")
        if(verbose)print(stat,right=TRUE)
        if(verbose)cat("\nBest iT based on BIC:",iT_best,"\n")
        
      }
      
      # Mixture specification
      if(iM>1){
        
        iT_max = max(iT)
        iT_min = min(iT)
        num_iT = iT_max-iT_min+1
        
        ll_seq          = rep(NA,num_iT)
        AIC_seq         = rep(NA,num_iT)
        BIClow_seq      = rep(NA,num_iT)
        BIChigh_seq     = rep(NA,num_iT)
        ICL_BIClow_seq  = rep(NA,num_iT)
        ICL_BIChigh_seq = rep(NA,num_iT)
        
        LCAseq = list()
        
        for(i in iT_min:iT_max){
          
          if(i==1){
            
            LCAseq[[1]] = list()
            
          } else{
            
            LCAseq[[1+i-iT_min]] = estlca(data=data,Y=Y,iT=i,id_high=id_high,iM=iM,
                                          Z=NULL,Zh=NULL,
                                          extout=FALSE,dataout=FALSE,kmea=TRUE,
                                          maxIter=1e3,tol=1e-8,fixed=TRUE,
                                          reord=TRUE,fixedpars=1,
                                          NRtol=1e-6,NRmaxit=100)
            
            ll_seq[1+i-iT_min]          = tail(LCAseq[[1+i-iT_min]]$LLKSeries,1)
            AIC_seq[1+i-iT_min]         = LCAseq[[1+i-iT_min]]$AIC
            BIClow_seq[1+i-iT_min]      = LCAseq[[1+i-iT_min]]$BIClow
            BIChigh_seq[1+i-iT_min]     = LCAseq[[1+i-iT_min]]$BIChigh
            ICL_BIClow_seq[1+i-iT_min]  = LCAseq[[1+i-iT_min]]$ICL_BIClow
            ICL_BIChigh_seq[1+i-iT_min] = LCAseq[[1+i-iT_min]]$ICL_BIChigh
            
          }
          
        }
        
        iT_best   = which.min(BIClow_seq)+iT_min-1
        est       = LCAseq[[1+iT_best-iT_min]]
        est$call  = match.call()
        
        if(verbose)cat("\nCALL:\n")
        if(verbose)print(match.call())
        if(verbose)cat("\n---------------------------\n")
        if(verbose)cat("\nGOODNESS-OF-FIT STATISTICS:\n")
        stat            = round(cbind(iT_min:iT_max,rep(iM,num_iT),ll_seq,AIC_seq,BIClow_seq,BIChigh_seq,ICL_BIClow_seq,ICL_BIChigh_seq),2)
        stat[is.na(stat)] = "-"
        rownames(stat)  = rep("",num_iT)
        colnames(stat)    = c("iT","iM","LL","AIC","BIClow","BIChigh","ICL_BIClow","ICL_BIChigh")
        if(verbose)print(noquote(stat),right=TRUE)
        if(verbose)cat("\nBest iT|iM =",iM,"based on BIClow:",iT_best,"\n")
        
      }
      
    }
    
  }
  
  # Model selection with respect to <iM>
  if(approach=="model selection on high"){
    
    # Non-mixture specification: one-low multilevel model
    if(iT==1){
      
      iM_max = max(iM)
      iM_min = min(iM)
      num_iM = iM_max-iM_min+1
      
      if(verbose)cat("\nCALL:\n")
      if(verbose)print(match.call())
      if(verbose)cat("\n---------------------------\n")
      if(verbose)cat("\nGOODNESS-OF-FIT STATISTICS:\n")
      stat              = cbind(rep(1,num_iM),iM_min:iM_max,matrix("-",num_iM,6))
      rownames(stat)    = rep("",num_iM)
      colnames(stat)    = c("iT","iM","LL","AIC","BIClow","BIChigh","ICL_BIClow","ICL_BIChigh")
      if(verbose)print(noquote(stat),right=TRUE)
      if(verbose)cat("\n")
      
      stop("Non-mixture specification.",call.=FALSE)
      
    }
    
    # Mixture specification
    if(iT>1){
      
      iM_max = max(iM)
      iM_min = min(iM)
      num_iM = iM_max-iM_min+1
      
      ll_seq          = rep(NA,num_iM)
      AIC_seq         = rep(NA,num_iM)
      BIClow_seq      = rep(NA,num_iM)
      BIChigh_seq     = rep(NA,num_iM)
      ICL_BIClow_seq  = rep(NA,num_iM)
      ICL_BIChigh_seq = rep(NA,num_iM)
      
      LCAseq = list()
      
      for(i in iM_min:iM_max){
        
        if(i==1){
          
          LCAseq[[1]] = estlca(data=data,Y=Y,iT=iT,id_high=NULL,iM=NULL,
                               Z=NULL,Zh=NULL,
                               extout=TRUE,dataout=FALSE,kmea=TRUE,
                               maxIter=1e3,tol=1e-8,fixed=TRUE,reord=TRUE,
                               fixedpars=1,NRtol=1e-6,NRmaxit=100,freqout=TRUE)
          
          iH    = ncol(data[,Y])
          iJ    = length(table(data[,id_high])[table(data[,id_high])>0])
          npar  = iT*iH+iT-1
          
          ll = tail(LCAseq[[1]]$LLKSeries,1)
          
          ll_seq[1]          = ll
          AIC_seq[1]         = -2*ll+2*npar
          BIClow_seq[1]      = LCAseq[[1]]$BIC
          BIChigh_seq[1]     = -2*ll+npar*log(iJ)
          ICL_BIClow_seq[1]  = LCAseq[[1]]$BIC+2*sum(LCAseq[[1]]$freq*apply(-LCAseq[[1]]$mU*log(LCAseq[[1]]$mU),1,sum))
          ICL_BIChigh_seq[1] = NA
          
        } else{
          
          LCAseq[[1+i-iM_min]] = estlca(data=data,Y=Y,iT=iT,id_high=id_high,iM=i,
                                        Z=NULL,Zh=NULL,
                                        extout=FALSE,dataout=FALSE,kmea=TRUE,
                                        maxIter=1e3,tol=1e-8,fixed=TRUE,
                                        reord=TRUE,fixedpars=1,
                                        NRtol=1e-6,NRmaxit=100)
          
          ll_seq[1+i-iM_min]          = tail(LCAseq[[1+i-iM_min]]$LLKSeries,1)
          AIC_seq[1+i-iM_min]         = LCAseq[[1+i-iM_min]]$AIC
          BIClow_seq[1+i-iM_min]      = LCAseq[[1+i-iM_min]]$BIClow
          BIChigh_seq[1+i-iM_min]     = LCAseq[[1+i-iM_min]]$BIChigh
          ICL_BIClow_seq[1+i-iM_min]  = LCAseq[[1+i-iM_min]]$ICL_BIClow
          ICL_BIChigh_seq[1+i-iM_min] = LCAseq[[1+i-iM_min]]$ICL_BIChigh
          
        }
        
      }
      
      iM_best   = which.min(BIChigh_seq)+iM_min-1
      est       = LCAseq[[1+iM_best-iM_min]]
      est$call  = match.call()
      
      if(verbose)cat("\nCALL:\n")
      if(verbose)print(match.call())
      if(verbose)cat("\n---------------------------\n")
      if(verbose)cat("\nGOODNESS-OF-FIT STATISTICS:\n")
      stat            = round(cbind(rep(iT,num_iM),iM_min:iM_max,ll_seq,AIC_seq,BIClow_seq,BIChigh_seq,ICL_BIClow_seq,ICL_BIChigh_seq),2)
      stat[is.na(stat)] = "-"
      rownames(stat)  = rep("",num_iM)
      colnames(stat)  = c("iT","iM","LL","AIC","BIClow","BIChigh","ICL_BIClow","ICL_BIChigh")
      if(verbose)print(noquote(stat),right=TRUE)
      if(verbose)cat("\nBest iM|iT =",iT,"based on BIChigh:",iM_best,"\n")
      
    }
    
  }
  
  # Model selection with respect to <iT> and <iM>
  if(approach=="model selection on low and high"){
    
    # Sequential model selection
    if(sequential==TRUE){
      
      iT_max = max(iT)
      iT_min = min(iT)
      num_iT = iT_max-iT_min+1
      
      iM_max = max(iM)
      iM_min = min(iM)
      num_iM = iM_max-iM_min+1
      
      # Step 1
      ll_step1          = rep(NA,num_iT)
      AIC_step1         = rep(NA,num_iT)
      BIClow_step1      = rep(NA,num_iT)
      BIChigh_step1     = rep(NA,num_iT)
      ICL_BIClow_step1  = rep(NA,num_iT)
      ICL_BIChigh_step1 = rep(NA,num_iT)
      
      for(i in iT_min:iT_max){
        
        if(i==1){
          
          iH = ncol(data[,Y])
          iN = nrow(data)
          iJ = length(table(data[,id_high])[table(data[,id_high])>0])
          
          ll = sum(apply(data[,Y],2,function(x){dbinom(x,1,mean(x),log=TRUE)}))
          
          ll_step1[1]          = ll
          AIC_step1[1]         = -2*ll+2*iH
          BIClow_step1[1]      = -2*ll+iH*log(iN)
          BIChigh_step1[1]     = -2*ll+iH*log(iJ)
          ICL_BIClow_step1[1]  = NA
          ICL_BIChigh_step1[1] = NA
          
        } else{
          
          LCAstep1 = estlca(data=data,Y=Y,iT=i,id_high=NULL,iM=NULL,
                            Z=NULL,Zh=NULL,
                            extout=TRUE,dataout=FALSE,kmea=TRUE,maxIter=1e3,
                            tol=1e-8,fixed=TRUE,reord=TRUE,fixedpars=1,
                            NRtol=1e-6,NRmaxit=100,freqout=TRUE)
          
          iH    = ncol(data[,Y])
          iJ    = length(table(data[,id_high])[table(data[,id_high])>0])
          npar  = i*iH+i-1
          
          ll = tail(LCAstep1$LLKSeries,1)
          
          ll_step1[1+i-iT_min]          = ll
          AIC_step1[1+i-iT_min]         = -2*ll+2*npar
          BIClow_step1[1+i-iT_min]      = LCAstep1$BIC
          BIChigh_step1[1+i-iT_min]     = -2*ll+npar*log(iJ)
          ICL_BIClow_step1[1+i-iT_min]  = LCAstep1$BIC+2*sum(LCAstep1$freq*apply(-LCAstep1$mU*log(LCAstep1$mU),1,sum))
          ICL_BIChigh_step1[1+i-iT_min] = NA
          
        }
        
      }
      
      iT_currbest = which.min(BIClow_step1)+iT_min-1
      
      # Step 2
      ll_step2          = rep(NA,num_iM)
      AIC_step2         = rep(NA,num_iM)
      BIClow_step2      = rep(NA,num_iM)
      BIChigh_step2     = rep(NA,num_iM)
      ICL_BIClow_step2  = rep(NA,num_iM)
      ICL_BIChigh_step2 = rep(NA,num_iM)
      
      LCAstep2 = list()
      
      if(iT_currbest==1){
        
        ll_step2[1]          = ll_step1[1+iT_currbest-iT_min]
        AIC_step2[1]         = AIC_step1[1+iT_currbest-iT_min]
        BIClow_step2[1]      = BIClow_step1[1+iT_currbest-iT_min]
        BIChigh_step2[1]     = BIChigh_step1[1+iT_currbest-iT_min]
        ICL_BIClow_step2[1]  = ICL_BIClow_step1[1+iT_currbest-iT_min]
        ICL_BIChigh_step2[1] = ICL_BIChigh_step1[1+iT_currbest-iT_min]
        
        iM_best = 1
        
      } else if(iT_currbest>1){
        
        for(i in iM_min:iM_max){
          
          if(i==1){
            
            ll_step2[1]          = ll_step1[1+iT_currbest-iT_min]
            AIC_step2[1]         = AIC_step1[1+iT_currbest-iT_min]
            BIClow_step2[1]      = BIClow_step1[1+iT_currbest-iT_min]
            BIChigh_step2[1]     = BIChigh_step1[1+iT_currbest-iT_min]
            ICL_BIClow_step2[1]  = ICL_BIClow_step1[1+iT_currbest-iT_min]
            ICL_BIChigh_step2[1] = ICL_BIChigh_step1[1+iT_currbest-iT_min]
            
          } else{
            
            LCAstep2[[1+i-iM_min]] = suppressWarnings(estlca(data=data,Y=Y,iT=iT_currbest,id_high=id_high,iM=i,
                                                             Z=NULL,Zh=NULL,
                                                             extout=FALSE,dataout=FALSE,kmea=TRUE,
                                                             maxIter=1e3,tol=1e-8,fixed=TRUE,
                                                             reord=TRUE,fixedpars=1,
                                                             NRtol=1e-6,NRmaxit=100))
            
            ll_step2[1+i-iM_min]          = tail(LCAstep2[[1+i-iM_min]]$LLKSeries,1)
            AIC_step2[1+i-iM_min]         = LCAstep2[[1+i-iM_min]]$AIC
            BIClow_step2[1+i-iM_min]      = LCAstep2[[1+i-iM_min]]$BIClow
            BIChigh_step2[1+i-iM_min]     = LCAstep2[[1+i-iM_min]]$BIChigh
            ICL_BIClow_step2[1+i-iM_min]  = LCAstep2[[1+i-iM_min]]$ICL_BIClow
            ICL_BIChigh_step2[1+i-iM_min] = LCAstep2[[1+i-iM_min]]$ICL_BIChigh
            
          }
          
        }
        
        iM_best = which.min(BIChigh_step2)+iM_min-1
        
      }
      
      # Step 3
      ll_step3          = rep(NA,num_iT)
      AIC_step3         = rep(NA,num_iT)
      BIClow_step3      = rep(NA,num_iT)
      BIChigh_step3     = rep(NA,num_iT)
      ICL_BIClow_step3  = rep(NA,num_iT)
      ICL_BIChigh_step3 = rep(NA,num_iT)
      
      LCAstep3 = list()
      
      if(iM_best==1){
        
        iT_best   = iT_currbest
        est       = estlca(data=data,Y=Y,iT=iT_best,id_high=NULL,iM=NULL,
                           Z=NULL,Zh=NULL,
                           extout=FALSE,dataout=FALSE,kmea=TRUE,maxIter=1e3,
                           tol=1e-8,fixed=TRUE,reord=TRUE,fixedpars=1,
                           NRtol=1e-6,NRmaxit=100)
        est$call  = match.call()
        
      } else if(iM_best>1){
        
        for(i in iT_min:iT_max){
          
          if(i==1){
            
            LCAstep3[[1]] = list()
            
          } else{
            
            LCAstep3[[1+i-iT_min]] = suppressWarnings(estlca(data=data,Y=Y,iT=i,id_high=id_high,iM=iM_best,
                                                             Z=NULL,Zh=NULL,
                                                             extout=FALSE,dataout=FALSE,kmea=TRUE,
                                                             maxIter=1e3,tol=1e-8,
                                                             fixed=TRUE,reord=TRUE,fixedpars=1,
                                                             NRtol=1e-6,NRmaxit=100))
            
            ll_step3[1+i-iT_min]          = tail(LCAstep3[[1+i-iT_min]]$LLKSeries,1)
            AIC_step3[1+i-iT_min]         = LCAstep3[[1+i-iT_min]]$AIC
            BIClow_step3[1+i-iT_min]      = LCAstep3[[1+i-iT_min]]$BIClow
            BIChigh_step3[1+i-iT_min]     = LCAstep3[[1+i-iT_min]]$BIChigh
            ICL_BIClow_step3[1+i-iT_min]  = LCAstep3[[1+i-iT_min]]$ICL_BIClow
            ICL_BIChigh_step3[1+i-iT_min] = LCAstep3[[1+i-iT_min]]$ICL_BIChigh
            
          }
          
        }
        
        iT_best   = which.min(BIClow_step3)+iT_min-1
        est       = LCAstep3[[1+iT_best-iT_min]]
        est$call  = match.call()
        
      }
      
      if(verbose)cat("\nCALL:\n")
      if(verbose)print(match.call())
      if(verbose)cat("\n------------------------------\n")
      if(verbose)cat("\nSTEP 1 - LOW-LEVEL CLUSTERING:\n")
      stat_step1                    = round(cbind(iT_min:iT_max,rep(NA,num_iT),ll_step1,AIC_step1,BIClow_step1,BIChigh_step1,ICL_BIClow_step1,ICL_BIChigh_step1),2)
      stat_step1[is.na(stat_step1)] = "-"
      rownames(stat_step1)  = rep("",num_iT)
      colnames(stat_step1)  = c("iT","iM","LL","AIC","BIClow","BIChigh","ICL_BIClow","ICL_BIChigh")
      if(verbose)print(noquote(stat_step1),right=TRUE)
      if(verbose)cat("\nCurrent best iT based on low-level BIC:",iT_currbest,"\n")
      if(verbose)cat("\n-------------------------------\n")
      if(verbose)cat("\nSTEP 2 - HIGH-LEVEL CLUSTERING:\n")
      stat_step2                    = round(cbind(rep(iT_currbest,num_iM),iM_min:iM_max,ll_step2,AIC_step2,BIClow_step2,BIChigh_step2,ICL_BIClow_step2,ICL_BIChigh_step2),2)
      stat_step2[is.na(stat_step2)] = "-"
      rownames(stat_step2)  = rep("",num_iM)
      colnames(stat_step2)  = c("iT","iM","LL","AIC","BIClow","BIChigh","ICL_BIClow","ICL_BIChigh")
      if(verbose)print(noquote(stat_step2),right=TRUE)
      if(verbose)cat("\nBest iM based on high-level BIC:",iM_best,"\n")
      if(verbose)cat("\n--------------------------------\n")
      if(verbose)cat("\nSTEP 3 - REVISITING LOWER LEVEL:\n")
      stat_step3                    = round(cbind(iT_min:iT_max,rep(iM_best,num_iT),ll_step3,AIC_step3,BIClow_step3,BIChigh_step3,ICL_BIClow_step3,ICL_BIChigh_step3),2)
      stat_step3[is.na(stat_step3)] = "-"
      rownames(stat_step3)  = rep("",num_iT)
      colnames(stat_step3)  = c("iT","iM","LL","AIC","BIClow","BIChigh","ICL_BIClow","ICL_BIChigh")
      if(verbose)print(noquote(stat_step3),right=TRUE)
      if(verbose)cat("\nBest iT based on low-level BIC:",iT_best,"\n")
      if(verbose)cat("\n---------------------------\n")
      if(verbose)cat("\nSUMMARY:\n")
      if(verbose)cat("\n", paste0("iT = ",iT_best,", iM = ",iM_best),"identified as the optimal model.\n")
      
    }
    
    # Parallelized model selection
    if(sequential==FALSE){
      
      datafoo     = data[,c(Y,id_high)]
      Yfoo        = Y
      id_highfoo  = id_high
      select_mat  = as.matrix(tidyr::expand_grid(iT,iM))
      nmod        = nrow(select_mat)
      iV          = parallel::detectCores()-numFreeCores
      cluster     = parallel::makeCluster(iV)
      parallel::clusterExport(cluster,c("datafoo","Yfoo","id_highfoo","select_mat","nmod"), envir = environment())
      parallel::clusterExport(cluster, 
                              unclass(lsf.str(envir = asNamespace("multilevLCA"), 
                                              all = T)),
                              envir = as.environment(asNamespace("multilevLCA"))
      )
      parallel::clusterEvalQ(cluster, {
        library(clustMixType)
        library(multilevLCA)
        library(mclust)
        library(tictoc)
        library(numDeriv)
        library(klaR)
        library(tidyverse)
        library(MASS)})
      simultaneous_out = parallel::parLapply(cluster,1:nmod,
                                             function(x){
                                               iFoo       = select_mat[x,]
                                               iT_curr    = iFoo[1]
                                               iM_curr    = iFoo[2]
                                               out_simult = simultsel(datafoo,Yfoo,iT_curr,id_highfoo,iM_curr)
                                               ll         = out_simult$ll
                                               AIC        = out_simult$AIC
                                               BIClow     = out_simult$BIClow
                                               BIChigh    = out_simult$BIChigh
                                               ICLlow     = out_simult$ICL_BIClow
                                               ICLhigh    = out_simult$ICL_BIChigh
                                               nout       = c("LL","AIC","BIClow","BIChigh","ICL_BIClow","ICL_BIChigh")
                                               out        = c(ll,AIC,BIClow,BIChigh,ICLlow,ICLhigh)
                                               names(out) = nout
                                               return(list(out=out,fit=out_simult$modFIT))
                                             }
      )
      parallel::stopCluster(cluster)
      stat = cbind(select_mat,apply(t(sapply(simultaneous_out,function(x){x$out})),2,function(x){round(as.numeric(x),2)}))
      
      mod_best  = which.min(stat[,5])
      est       = simultaneous_out[[mod_best]]$fit
      est$call  = match.call()
      
      if(verbose)cat("\nCALL:\n")
      if(verbose)print(match.call())
      if(verbose)cat("\n---------------------------\n")
      if(verbose)cat("\nGOODNESS-OF-FIT STATISTICS:\n")
      stat[is.na(stat)] = "-"
      rownames(stat)  = rep("",nrow(stat))
      if(verbose)print(noquote(stat),right=TRUE)
      mod_best_print = paste0("iT = ",stat[mod_best,1],",iM = ",stat[mod_best,2])
      if(verbose)cat("\nBest iT,iM based on BIClow:",mod_best_print,"\n")
      
    }
    
  }
  
  #################
  ### 6. Return ###
  #################
  
  class(est) = "multiLCA"
  return(est)
  
}

print.multiLCA = function(x,...){
  
  out = NULL
  
  if(x$spec == "Single-level LC measurement model"){
    
    ######################################
    ### Single-level measurement model ###
    ######################################
    
    cat("\nCALL:\n")
    print(x$call)
    cat("\nSPECIFICATION:\n")
    print(noquote(x$spec))
    
    iter            = t(matrix(c(x$iter, x$LLKSeries[1], tail(x$LLKSeries, 1))))
    colnames(iter)  = c("EMiter", "LLfirst", "LLlast")
    rownames(iter)  = ""
    
    cat("\nESTIMATION DETAILS:\n\n")
    print(iter)
    cat("\n---------------------------\n")
    cat("\nCLASS PROPORTIONS:\n")
    print(round(x$vPg, 4))
    cat("\nRESPONSE PROBABILITIES:\n\n")
    print(round(x$mPhi, 4))
    cat("\n---------------------------\n")
    cat("\nMODEL AND CLASSIFICATION STATISTICS:\n")
    
    stat            = as.matrix(c(as.character(round(x$AIC, 4)), as.character(round(x$BIC, 4)), as.character(round(x$AvgClassErrProb, 4)), as.character(round(x$R2entr, 4))))
    rownames(stat)  = c("AIC", "BIC", "ClassErr", "EntR-sqr")
    colnames(stat)  = ""
    
    print(noquote(stat), right = TRUE)
    cat("\n")
    
  } else if(x$spec == "Single-level LC structural model"){
    
    #####################################
    ### Single-level structural model ###
    #####################################
    
    cat("\nCALL:\n")
    print(x$call)
    cat("\nSPECIFICATION:\n")
    print(noquote(x$spec))
    
    iter            = t(matrix(c(x$iter, x$LLKSeries[1], tail(x$LLKSeries, 1))))
    colnames(iter)  = c("EMiter", "LLfirst", "LLlast")
    rownames(iter)  = ""
    
    cat("\nESTIMATION DETAILS:\n\n")
    print(iter)
    cat("\n---------------------------\n")
    cat("\nCLASS PROPORTIONS (SAMPLE MEAN):\n")
    print(round(x$vPg_avg, 4))
    cat("\nRESPONSE PROBABILITIES:\n\n")
    print(round(x$mPhi, 4))
    cat("\n---------------------------\n")
    cat("\nMODEL AND CLASSIFICATION STATISTICS:\n")
    
    stat             = as.matrix(c(as.character(round(x$AIC, 4)), as.character(round(x$BIC, 4)), as.character(round(x$AvgClassErrProb, 4)), as.character(round(x$R2entr, 4))))
    rownames(stat)   = c("AIC", "BIC", "ClassErr", "EntR-sqr")
    colnames(stat)   = ""
    
    print(noquote(stat), right = TRUE)
    cat("\n")
    
    cat("\n---------------------------\n")
    cat("\nLOGISTIC MODEL FOR CLASS MEMBERSHIP:\n\n")
    
    betas   = format(round(as.matrix(x$betas), 4), nsmall = 4)
    SE      = format(round(as.matrix(x$SEs_cor_beta), 4), nsmall = 4)
    Zscore  = format(round(as.matrix(x$betas)/as.matrix(x$SEs_cor_beta), 4), nsmall = 4)
    pval    = format(round(2*(1 - pnorm(abs(as.matrix(x$betas)/as.matrix(x$SEs_cor_beta)))), 4), nsmall = 4)
    psignf  = matrix("   ", nrow(betas), ncol(betas))
    psignf[pval < 0.1 & pval >= 0.05]   = "*  "
    psignf[pval < 0.05 & pval >= 0.01]  = "** "
    psignf[pval < 0.01]                 = "***"
    
    for (i in 1:ncol(betas)){
      
      C = paste0("C", 1 + i)
      logit_params = noquote(cbind(betas[,i], SE[,i], Zscore[,i], matrix(paste0(pval[,i], psignf[,i]))))
      colnames(logit_params) = c("Beta", "S.E.", "Z-score", "p-value")
      rownames(logit_params) = paste0(substr(rownames(as.matrix(x$betas)), 1, nchar(rownames(as.matrix(x$betas))) - 1), 1 + i, ")")
      
      cat("\nMODEL FOR", C, "(BASE C1)\n\n")
      print(logit_params, right = TRUE)
      cat("\n")
      cat("\n", "*** p < 0.01, ** p < 0.05, * p < 0.1")
      cat("\n")
    }
    
  } else if(x$spec == "Multilevel LC measurement model"){
    
    ####################################
    ### Multilevel measurement model ###
    ####################################
    
    cat("\nCALL:\n")
    print(x$call)
    cat("\nSPECIFICATION:\n")
    print(noquote(x$spec))
    
    iter            = t(matrix(c(x$iter, x$LLKSeries[1], tail(x$LLKSeries, 1))))
    colnames(iter)  = c("EMiter", "LLfirst", "LLlast")
    rownames(iter)  = ""
    
    cat("\nESTIMATION DETAILS:\n\n")
    print(iter)
    cat("\n---------------------------\n")
    cat("\nGROUP PROPORTIONS:\n")
    print(round(x$vOmega, 4))
    cat("\nCLASS PROPORTIONS:\n")
    print(round(x$mPi, 4))
    cat("\nRESPONSE PROBABILITIES:\n\n")
    print(round(x$mPhi, 4))
    cat("\n---------------------------\n")
    cat("\nMODEL AND CLASSIFICATION STATISTICS:\n")
    
    stat = as.matrix(c(as.character(round(x$AIC, 4)), as.character(round(x$BIClow, 4)), as.character(round(x$BIChigh, 4)), as.character(round(x$ICL_BIClow, 4)), as.character(round(x$ICL_BIChigh, 4)), as.character(round(x$R2entr_low, 4)), as.character(round(x$R2entr_high, 4))))
    rownames(stat) = c("AIC", "BIClow", "BIChigh", "ICLBIClow", "ICLBIChigh", "R2entrlow", "R2entrhigh")
    colnames(stat) = ""
    
    print(noquote(stat), right = TRUE)
    
  } else if(x$spec == "Multilevel structural LC model with low-level covariates"){
    
    #############################################################
    ### Multilevel structural model with low-level covariates ###
    #############################################################
    
    cat("\nCALL:\n")
    print(x$call)
    cat("\nSPECIFICATION:\n")
    print(noquote(x$spec))
    
    iter            = t(matrix(c(x$iter, x$LLKSeries[1], tail(x$LLKSeries, 1))))
    colnames(iter)  = c("EMiter", "LLfirst", "LLlast")
    rownames(iter)  = ""
    
    cat("\nESTIMATION DETAILS:\n\n")
    print(iter)
    cat("\n---------------------------\n")
    cat("\nGROUP PROPORTIONS:\n")
    print(round(x$vOmega, 4))
    cat("\nCLASS PROPORTIONS (SAMPLE MEAN):\n\n")
    print(round(x$mPi_avg, 4))
    cat("\nRESPONSE PROBABILITIES:\n\n")
    print(round(x$mPhi, 4))
    cat("\n---------------------------\n")
    cat("\nMODEL AND CLASSIFICATION STATISTICS:\n")
    
    stat = as.matrix(c(as.character(round(x$AIC, 4)), as.character(round(x$BIClow, 4)), as.character(round(x$BIChigh, 4)), as.character(round(x$ICL_BIClow, 4)), as.character(round(x$ICL_BIChigh, 4)), as.character(round(x$R2entr_low, 4)), as.character(round(x$R2entr_high, 4))))
    rownames(stat) = c("AIC", "BIClow", "BIChigh", "ICLBIClow", "ICLBIChigh", "R2entrlow", "R2entrhigh")
    colnames(stat) = ""
    
    print(noquote(stat), right = TRUE)
    cat("\n")
    
    cat("\n---------------------------\n")
    cat("\nLOGISTIC MODEL FOR LOW-LEVEL CLASS MEMBERSHIP:\n\n")
    
    for (i in 1:dim(x$cGamma)[3]){
      
      gammas  = format(round(as.matrix(x$cGamma[,,i]), 4), nsmall = 4)
      SE      = format(round(as.matrix(x$SEs_cor_gamma[,,i]), 4), nsmall = 4)
      Zscore  = format(round(as.matrix(x$cGamma[,,i])/as.matrix(x$SEs_cor_gamma[,,i]), 4), nsmall = 4)
      pval    = format(round(2*(1 - pnorm(abs(as.matrix(x$cGamma[,,i])/as.matrix(x$SEs_cor_gamma[,,i])))), 4), nsmall = 4)
      psignf  = matrix("   ", nrow(gammas), ncol(gammas))
      psignf[pval < 0.1 & pval >= 0.05]   = "*  "
      psignf[pval < 0.05 & pval >= 0.01]  = "** "
      psignf[pval < 0.01]                 = "***"
      
      for (j in 1:ncol(gammas)){
        
        G = paste0("G", i)
        C = paste0("C", 1 + j)
        logit_params = noquote(cbind(gammas[,j], SE[,j], Zscore[,j], matrix(paste0(pval[,j], psignf[,j]))))
        colnames(logit_params) = c("Gamma", "S.E.", "Z-score", "p-value")
        rownames(logit_params) = paste0(substr(rownames(as.matrix(x$cGamma[,,i])), 1, nchar(rownames(as.matrix(x$cGamma[,,i]))) - 1), j + 1, "|G", i, ")")
        
        cat("\nMODEL FOR", C, "(BASE C1) GIVEN", G, "\n\n")
        print(logit_params, right = TRUE)
        cat("\n")
        cat("\n", "*** p < 0.01, ** p < 0.05, * p < 0.1")
        cat("\n")
        
      }
      
    }
    
  } else if(x$spec == "Multilevel structural LC model with low- and high-level covariates"){
    
    #######################################################################
    ### Multilevel structural model with low- and high-level covariates ###
    #######################################################################
    
    cat("\nCALL:\n")
    print(x$call)
    cat("\nSPECIFICATION:\n")
    print(noquote(x$spec))
    
    iter            = t(matrix(c(x$iter, x$LLKSeries[1], tail(x$LLKSeries, 1))))
    colnames(iter)  = c("EMiter", "LLfirst", "LLlast")
    rownames(iter)  = ""
    
    cat("\nESTIMATION DETAILS:\n\n")
    print(iter)
    cat("\n---------------------------\n")
    cat("\nGROUP PROPORTIONS (SAMPLE MEAN):\n")
    print(round(x$vOmega_avg, 4))
    cat("\nCLASS PROPORTIONS (SAMPLE MEAN):\n\n")
    print(round(x$mPi_avg, 4))
    cat("\nRESPONSE PROBABILITIES:\n\n")
    print(round(x$mPhi, 4))
    cat("\n---------------------------\n")
    cat("\nMODEL AND CLASSIFICATION STATISTICS:\n")
    
    stat = as.matrix(c(as.character(round(x$AIC, 4)), as.character(round(x$BIClow, 4)), as.character(round(x$BIChigh, 4)), as.character(round(x$ICL_BIClow, 4)), as.character(round(x$ICL_BIChigh, 4)), as.character(round(x$R2entr_low, 4)), as.character(round(x$R2entr_high, 4))))
    rownames(stat) = c("AIC", "BIClow", "BIChigh", "ICLBIClow", "ICLBIChigh", "R2entrlow", "R2entrhigh")
    colnames(stat) = ""
    
    print(noquote(stat), right = TRUE)
    cat("\n")
    
    cat("\n---------------------------\n")
    cat("\nLOGISTIC MODEL FOR HIGH-LEVEL CLASS MEMBERSHIP:\n\n")
    
    deltas  = format(round(as.matrix(x$mDelta), 4), nsmall = 4)
    SE      = format(round(as.matrix(x$SEs_cor_delta), 4), nsmall = 4)
    Zscore  = format(round(as.matrix(x$mDelta)/as.matrix(x$SEs_cor_delta), 4), nsmall = 4)
    pval    = format(round(2*(1 - pnorm(abs(as.matrix(x$mDelta)/as.matrix(x$SEs_cor_delta)))), 4), nsmall = 4)
    psignf  = matrix("   ", nrow(deltas), ncol(deltas))
    psignf[pval < 0.1 & pval >= 0.05]   = "*  "
    psignf[pval < 0.05 & pval >= 0.01]  = "** "
    psignf[pval < 0.01]                 = "***"
    
    for (i in 1:ncol(deltas)){
      
      G = paste0("G", 1 + i)
      logit_params = noquote(cbind(deltas[,i], SE[,i], Zscore[,i], matrix(paste0(pval[,i], psignf[,i]))))
      colnames(logit_params) = c("Delta", "S.E.", "Z-score", "p-value")
      rownames(logit_params) = paste0(substr(rownames(as.matrix(x$mDelta)), 1, nchar(rownames(as.matrix(x$mDelta))) - 1), 1 + i, ")")
      
      cat("\nMODEL FOR", G, "(BASE G1)\n\n")
      print(logit_params, right = TRUE)
      cat("\n")
      cat("\n", "*** p < 0.01, ** p < 0.05, * p < 0.1")
      cat("\n")
      
    }
    
    cat("\n---------------------------\n")
    cat("\nLOGISTIC MODEL FOR LOW-LEVEL CLASS MEMBERSHIP:\n\n")
    
    for (i in 1:dim(x$cGamma)[3]){
      
      gammas  = format(round(as.matrix(x$cGamma[,,i]), 4), nsmall = 4)
      SE      = format(round(as.matrix(x$SEs_cor_gamma[,,i]), 4), nsmall = 4)
      Zscore  = format(round(as.matrix(x$cGamma[,,i])/as.matrix(x$SEs_cor_gamma[,,i]), 4), nsmall = 4)
      pval    = format(round(2*(1 - pnorm(abs(as.matrix(x$cGamma[,,i])/as.matrix(x$SEs_cor_gamma[,,i])))), 4), nsmall = 4)
      psignf  = matrix("   ", nrow(gammas), ncol(gammas))
      psignf[pval < 0.1 & pval >= 0.05]   = "*  "
      psignf[pval < 0.05 & pval >= 0.01]  = "** "
      psignf[pval < 0.01]                 = "***"
      
      for (j in 1:ncol(gammas)){
        
        G = paste0("G", i)
        C = paste0("C", 1 + j)
        logit_params = noquote(cbind(gammas[,j], SE[,j], Zscore[,j], matrix(paste0(pval[,j], psignf[,j]))))
        colnames(logit_params) = c("Gamma", "S.E.", "Z-score", "p-value")
        rownames(logit_params) = paste0(substr(rownames(as.matrix(x$cGamma[,,i])), 1, nchar(rownames(as.matrix(x$cGamma[,,i]))) - 1), j + 1, "|G", i, ")")
        
        cat("\nMODEL FOR", C, "(BASE C1) GIVEN", G, "\n\n")
        print(logit_params, right = TRUE)
        cat("\n")
        cat("\n", "*** p < 0.01, ** p < 0.05, * p < 0.1")
        cat("\n")
        
      }
      
    }
    
  }
  
}
plot.multiLCA = function(x, horiz = TRUE, clab = NULL, ...){
  
  out = NULL
  
  iT        = ncol(x$mPhi)
  items     = 1:nrow(x$mPhi)
  itemnames = substr(rownames(x$mPhi), 3, nchar(rownames(x$mPhi)) - 3)
  
  if(horiz == TRUE){
    
    las = 1
    
  } else{
    
    las = 3
    
  }
  if(is.null(clab)){
    
    legend = paste0("Class ", 1:iT)
    
  } else{
    
    if(!is.vector(clab) | is.list(clab)){
      
      stop("Invalid clab.",call.=FALSE)
      
    } else if(length(clab) != iT | !is.character(clab)){
      
      stop("Invalid clab.",call.=FALSE)
      
    }
    
    legend = clab
    
  }
  
  oldpar = par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mar = c(5, 5, 3, 8), xpd = TRUE)
  plot(items, x$mPhi[,1], type = "b", xaxt = "n", pch = 1, 
       xlab = "", ylab = "Response probability", frame.plot = FALSE, ylim = c(0,1), lty = 1, ...)
  axis(1, at = items, labels = itemnames, las = las)
  for(h in 2:iT){
    lines(items, x$mPhi[,h], pch = h, type = "b", lty = h)
  }
  legend("topright", legend = legend,
         lty = (1:iT), inset = c(-0.115, 0), pch = 1:iT, cex = 0.8)
  
  
}
