MLCAsimulate = function(n_idlow = 1000, K = 5, nj = 3, J = 2, n_idhigh = 10,
                        ncov_ind = 2,
                        beta0 = log(0.9/0.1), beta1 = matrix(0,K,nj),
                        cGamma = array(0,c(nj-1,ncov_ind+1,J)),
                        delta0 =rep(0,J-1), delta1=rep(0,J-1),mLog_sep=NULL){
  # N is the total sample size
  # n_idhigh is the number of high lev units
  # K is the number of indicators
  # nj is the number of Classes
  # J is the number of GClasses
  # beta0 is the intercept in the item equations
  # beta1 are the DEs of the covariate on items (should be of length K. put zero for non dif items)
  # cGamma is the array of regression coefficients (including intercept)
  # of all Class | Gclass equations (nj-1 x J)
  # delta0 is the intercept in the GClass equation
  #  Z is the covariate
  # wpath should be written as wpath="D:\\Dropbox\\robust 2_3 step\\multilevel LCA\\"

  N = n_idlow*n_idhigh
  # Z = matrix(rnorm(N),n_idlow,n_idhigh)
  Z=matrix(rnorm(N*ncov_ind),N,ncov_ind)
  Z=cbind(1,Z)
  Z_high = rnorm(n_idhigh)
  Wprob = matrix(1,n_idhigh,J)
  for(ih in 1:n_idhigh){
    Wprob[ih,-1] = exp(delta0 + delta1*Z_high[ih])
    Wprob[ih,] = Wprob[ih,]/sum(Wprob[ih,])
  }

  w = numeric(n_idhigh)
  if(sum(delta1==0)==0){
    w = sample(J,n_idhigh,prob=Wprob[1,],replace=T)
  }
  else{
    for(ih in 1:n_idhigh){
      w[ih] = sample(J,1,prob=Wprob[ih,])
    }
  }
  control = (1:J) %in% w
  if(sum(control) != J){
    classchng = sample(1:n_idhigh,size=sum(!control),replace=F)
    foo=0
    for(l in classchng){
      foo = foo + 1
      w[l] = foo
    }
  }
  idhigh_long   = rep(1:n_idhigh,each=n_idlow)
  W = w[idhigh_long]
  X = numeric(N)
  Zf_high = rep(Z_high,each=n_idlow)

  probXfoo = array(1,c(N,nj,J))
  X = numeric(N)
  for(j in 1:J){
    for(n in 2:nj){
      probXfoo[,n,j] = exp(Z%*%cGamma[n-1,,j]) #
    }
    probXfoo[,,j] = sweep(probXfoo[,,j],1,apply(probXfoo[,,j],1,sum),"/")
  }
  for(i in 1:N){
    X[i] = sample(nj,1,prob = probXfoo[i,,W[i]])
  }

  if(is.null(mLog_sep)){
    mLog_sep = matrix(beta0,K,nj)
    if(nj == 2){
      mLog_sep[,2] = -beta0
    }
    if(nj>2){
      for(n in 2:nj){
        foo = floor(((n+1)/nj)*(K-2))
        mLog_sep[1:foo,n] = -mLog_sep[1:foo,n]
      }
    }
  }
  Y = matrix(0,N,K)
  resfull = Y
  for(k in 1:K){
    for(i in 1:N){
      resprob = exp(mLog_sep[k,X[i]] + beta1[k,X[i]]*Z[i])
      # resprob = exp(mLog_sep[k,Xvec[i]] + beta1[k,Xvec[i]]*Z[i])
      resprob = resprob/(1+resprob)
      Y[i,k] = rbinom(1,1,resprob)
      resfull[i,k] = resprob;
    }
  }
  data = data.frame(Y=Y,Z=Z,Zf_high = Zf_high, W=W,
                    idhigh=idhigh_long,X=X)
  parvec = c(delta0,c(cGamma),c(mLog_sep))
  return(list(idhigh_long = idhigh_long,
              Y=Y,
              Z=Z,
              Z_high =Z_high,
              Wprob = Wprob,
              Xprob = probXfoo,
              resfullY = resfull,
              w=w,
              W=W,
              X=X,
              mLog_sep=mLog_sep,
              beta1=beta1,
              cGamma = cGamma,
              delta0=delta0, delta1=delta1,
              parvec = parvec)
  )

}



# 
# MLCAsimulate = function(n_idlow = 1000, K = 5, nj = 3, J = 2, n_idhigh = 10,
#                         ncov_ind = 2,
#                         beta0 = log(0.9/0.1), beta1 = matrix(0,K,nj),
#                         cGamma = array(0,c(nj-1,ncov_ind+1,J)),
#                         delta0 =rep(0,J-1), delta1=rep(0,J-1),mLog_sep=NULL){
#   # N is the total sample size
#   # n_idhigh is the number of high lev units
#   # K is the number of indicators
#   # nj is the number of Classes
#   # J is the number of GClasses
#   # beta0 is the intercept in the item equations
#   # beta1 are the DEs of the covariate on items (should be of length K. put zero for non dif items)
#   # cGamma is the array of regression coefficients (including intercept)
#   # of all Class | Gclass equations (nj-1 x J)
#   # delta0 is the intercept in the GClass equation
#   #  Z is the covariate
#   # wpath should be written as wpath="D:\\Dropbox\\robust 2_3 step\\multilevel LCA\\"
# 
#   N = n_idlow*n_idhigh
#   # Z = matrix(rnorm(N),n_idlow,n_idhigh)
#   Z=matrix(rnorm(N*ncov_ind),N,ncov_ind)
#   Z=cbind(1,Z)
#   Z_high = rnorm(n_idhigh)
#   cGamma_foo = array(0,c(nj,ncov_ind+1,J))
#   cGamma_foo[-1,,] = cGamma
#   cGamma = cGamma_foo
#   delta0 = c(0,delta0)
#   delta1 = c(0,delta1)
#   Wprob = matrix(1,n_idhigh,J)
#   for(ih in 1:n_idhigh){
#     Wprob[ih,] = exp(delta0 + delta1*Z_high[ih])
#     Wprob[ih,] = Wprob[ih,]/sum(Wprob[ih,])
#   }
#   #
#   w = numeric(n_idhigh)
#   if(sum(delta1==0)==0){
#     w = sample(J,n_idhigh,prob=Wprob[1,],rep=T)
#   }
#   else{
#     for(ih in 1:n_idhigh){
#       w[ih] = sample(J,1,prob=Wprob[ih,])
#     }
#   }
#   control = (1:J) %in% w
#   if(sum(control) != J){
#     classchng = sample(1:n_idhigh,size=sum(!control),rep=F)
#     foo=0
#     for(l in classchng){
#       foo = foo + 1
#       w[l] = foo
#     }
#   }
#   idhigh_long   = rep(1:n_idhigh,each=n_idlow)
#   W = w[idhigh_long]
#   X = numeric(N)
#   Zf_high = rep(Z_high,each=n_idlow)
# 
#   probXfoo = array(1,c(N,nj,J))
#   X = numeric(N)
#   for(j in 1:J){
#     for(n in 1:nj){
#       probXfoo[,n,j] = exp(Z%*%cGamma[n,,j]) #
#     }
#     probXfoo[,,j] = sweep(probXfoo[,,j],1,apply(probXfoo[,,j],1,sum),"/")
#   }
#   for(i in 1:N){
#     X[i] = sample(nj,1,prob = probXfoo[i,,W[i]])
#   }
# 
#   if(is.null(mLog_sep)){
#     mLog_sep = matrix(beta0,K,nj)
#     if(nj == 2){
#       mLog_sep[,2] = -beta0
#     }
#     if(nj>2){
#       for(n in 2:nj){
#         foo = floor(((n+1)/nj)*(K-2))
#         mLog_sep[1:foo,n] = -mLog_sep[1:foo,n]
#       }
#     }
#   }
#   Y = matrix(0,N,K)
#   resfull = Y
#   for(k in 1:K){
#     for(i in 1:N){
#       resprob = exp(mLog_sep[k,X[i]] + beta1[k,X[i]]*Z[i])
#       # resprob = exp(mLog_sep[k,Xvec[i]] + beta1[k,Xvec[i]]*Z[i])
#       resprob = resprob/(1+resprob)
#       Y[i,k] = rbinom(1,1,resprob)
#       resfull[i,k] = resprob;
#     }
#   }
#   data = data.frame(Y=Y,Z=Z,Zf_high = Zf_high, W=W,
#                     idhigh=idhigh_long,X=X)
#   parvec = c(delta0[-1],as.vector(cGamma_foo[-1,,]),as.vector(mLog_sep))
#   return(list(idhigh_long = idhigh_long,
#               Y=Y,
#               Z=Z,
#               Z_high =Z_high,
#               Wprob = Wprob,
#               Xprob = probXfoo,
#               resfullY = resfull,
#               w=w,
#               W=W,
#               X=X,
#               mLog_sep=mLog_sep,
#               beta1=beta1,
#               cGamma = cGamma_foo[-1,,],
#               delta0=delta0, delta1=delta1,
#               parvec = parvec)
#   )
# 
# }


LCAsimulate = function(N = 1000, K = 5, J = 2, 
                       beta0 = log(0.9/0.1), beta1 = matrix(0,K,J), 
                       gamma0 = rep(0,J-1), gamma1 = rep(0,J-1)){
  # N is the sample size
  # K is the number of indicators
  # J is the number of Classes
  # beta0 is the intercept in the item equations
  # beta1 are the DEs of the covariate on items (should be of length K. put zero for non dif items)
  # gamma0 is the intercept in the Class equation (J x (nj -1) matrix)
  # gamma1 is the covariate effect in the Class equation
  #  Z is the covariate
  Z=rnorm(N)
  # gamma0 = cbind(0,gamma0)
  # gamma1 = c(0,gamma1)
  X = numeric(N)
  
  probXfoo = matrix(1,N,J)
  X = numeric(N)
  for(j in 2:J){
    probXfoo[,j] = exp(gamma0[j-1] + gamma1[j-1]*Z) 
  }
  probXfoo = sweep(probXfoo,1,apply(probXfoo,1,sum),"/")
  
  for(i in 1:N){
    X[i] = sample(J,1,prob = probXfoo[i,])
  }
  
  mLog_sep = matrix(beta0,K,J)
  if(J == 2){
    mLog_sep[,2] = -beta0
  }
  if(J>2){
    for(n in 2:J){
      foo = floor(((n+1)/J)*(K-2))
      mLog_sep[1:foo,n] = -mLog_sep[1:foo,n]
    }
  }
  Y = matrix(0,N,K)
  resfull = Y
  for(k in 1:K){
    for(i in 1:N){
      resprob = exp(mLog_sep[k,X[i]] + beta1[k,X[i]]*Z[i])
      # resprob = exp(mLog_sep[k,Xvec[i]] + beta1[k,Xvec[i]]*Z[i])
      resprob = resprob/(1+resprob)
      Y[i,k] = rbinom(1,1,resprob)
      resfull[i,k] = resprob;
    }
  }
  data = data.frame(Y=Y,X=Z,Cluster=X)
  return(list(Y=Y,
              Z=Z,
              data=data,
              Xprob = probXfoo,
              resfullY = resfull,
              X=X,
              mLog_sep=mLog_sep,
              beta1=beta1,
              gamma0=gamma0, gamma1=gamma1)
  )
  
}