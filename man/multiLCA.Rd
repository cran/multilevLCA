\name{multiLCA}
\alias{multiLCA}

\title{
Estimates and plots single- and multilevel latent class models
}
\description{
The \code{multiLCA} function in the \code{multilevLCA} package estimates single- and multilevel measurement and structural latent class models. Moreover, the function performs two different strategies for model selection. Methodological details can be found in Bakk et al. (2022), Bakk and Kuha (2018), and Di Mari et al. (2023).

Different output visualization tools are available for all model specifications. See, e.g., \code{\link{plot.multiLCA}}.
}
\usage{
multiLCA(data, 
Y, 
iT, 
id_high = NULL,
iM = NULL,
Z = NULL,
Zh = NULL,
extout = FALSE,
dataout = FALSE,
kmea = TRUE,
sequential = TRUE,
numFreeCores = 2,
maxIter = 1e3,
tol = 1e-8,
fixed = TRUE,
reord = TRUE,
fixedpars = 1,
NRtol = 1e-6,
NRmaxit = 100,
verbose = TRUE)
}

\arguments{
  \item{data}{
     Input matrix or dataframe
}
  \item{Y}{
     Names of \code{data} columns with 0-1 coded items
}
  \item{iT}{
     Number of low-level clusters
}
  \item{id_high}{
     Name of \code{data} column with high-level id. Default \code{NULL}
}
  \item{iM}{
     Number of high-level clusters. Default \code{NULL}
}
  \item{Z}{
     Names of \code{data} columns with low-level covariates. Default \code{NULL}
}
  \item{Zh}{
     Names of \code{data} columns with high-level covariates. Default \code{NULL}
}
  \item{extout}{
     Whether to output extensive model and estimation statistics. Default \code{FALSE}
}
  \item{dataout}{
     Whether to output the cleaned dataset on which estimation was performed. Default \code{FALSE}
}
  \item{kmea}{
     Whether to compute starting values for single-level estimation using K-means (\code{TRUE}), recommended for stability, or K-modes (\code{FALSE}). Default \code{TRUE}
}
  \item{sequential}{
     Whether to perform sequential (\code{TRUE}) or parallelized (\code{FALSE}) model selection. Default \code{TRUE}
}
  \item{numFreeCores}{
     Number of CPU cores to keep free in parallelized model selection. Default 2
}
  \item{maxIter}{
     Maximum number of iterations for EM algorithm in single-level model estimation. Default 1e3
}
  \item{tol}{
     EM tolerance for measurement model estimation. Default 1e-8
}
  \item{fixed}{
     Whether to fix (\code{TRUE}) or update (\code{FALSE}) measurement model parameters in regression stage of single-level structural model estimation. Default \code{TRUE}
}
  \item{reord}{
     Whether to (re)order classes in decreasing order according to probability of scoring yes on all items. Default \code{TRUE}
}
  \item{fixedpars}{
     Estimator of multilevel model; one-step estimator (0), two-step estimator (1) or three-step estimator (2). Default 1
}
  \item{NRtol}{
     Newton-Raphson tolerance for structural model estimation. Default 1e-6
}
  \item{NRmaxit}{
     Maximum number of iterations for Newton-Raphson algorithm in multilevel model estimation. Default 100
}
  \item{verbose}{
     Whether to print model selection results. Default TRUE
}

}
\details{
To directly estimate a latent class model, \code{iT} and (optionally) \code{iM} should be specified as a single positive integer. To perform model selection over a set of number of classes, \code{iT} and/or (optionally) \code{iM} can be specified as a range of consecutive positive integers of the form \code{iT_min:iT_max} and \code{iM_min:iM_max}, respectively. It is possible to specify \code{iT = iT_min:iT_max} with either \code{iM = NULL} or \code{iM} equal to a single positive integer, \code{iM = iM_min:iM_max} with \code{iT} equal to a single positive integer, or \code{iT = iT_min:iT_max} with \code{iM = iM_min:iM_max}. All model selection procedures print the results and returns the output of the optimal model based on BIC.

In the case where both \code{iT} and \code{iM} are defined as a range of consecutive positive integers, model selection can be performed using the sequential three-stage approach (Lukociene et al., 2010) or a simultaneous approach. The choice between the two is indicated using the \code{sequential} argument, with default \code{TRUE} for sequential model selection. The sequential approach involves (first step:) estimating \code{iT_min:iT_max} single-level models and identifying the optimal alternative \code{iT_opt1} based on BIC, (second step:) estimating \code{iM_min:iM_max|iT = iT_opt1} multilevel models and identifying the optimal alternative \code{iM_opt2} based on high-level BIC, and finally (third step:) estimating \code{iT_min:iT_max|iM = iM_opt2} multilevel models and identifying the optimal alternative \code{iT_opt3} based on low-level BIC. The simultaneous approach involves devoting multiple CPU cores on the local machine to estimate all combinations in \code{iT = iT_min:iT_max, iM = iM_min:iM_max} and identifying the optimal alternative based on low-level BIC. The number of CPU cores to keep free during estimation is indicated using the \code{numFreeCores} argument with default 2.
}
\value{
Single-level measurement model estimation returns (if \code{extout = TRUE}, a subset):
  \item{vPg}{Class proportions}
  \item{mPhi}{Conditional response probabilities}
  \item{alphas}{Intercept parameters in logistic models for class proportions}
  \item{gammas}{Intercept parameters in logistic models for response probabilities}
  \item{parvec}{Vector of all model parameters}
  \item{Infomat}{Expected information matrix}
  \item{Varmat}{Variance-covariance matrix}
  \item{SEs}{Standard errors}
  \item{eps}{Epsilon: difference between the last two elements of the EM log-likelihood sequence}
  \item{mU}{Posterior class membership probabilities}
  \item{mU_modal}{Modal class assignments as matrix}
  \item{vU_modal}{Modal class assignments as vector}
  \item{mClassErr}{Expected number of classification errors}
  \item{mClassErrProb}{Expected proportion of classification errors}
  \item{AvgClassErrProb}{Average expected proportion of classification errors}
  \item{R2entr}{Entropy R\eqn{^2}}
  \item{AIC}{Akaike Information Criterion}
  \item{BIC}{Bayesian Information Criterion}
  \item{iter}{Number of iterations of EM algorithm}
  \item{LLKSeries}{Full log-likelihood series of EM algorithm}
  \item{mScore}{Individual contributions to log-likelihood score}
  \item{spec}{Specification}
Single-level structural model estimation returns (if \code{extout = TRUE}, a subset):
  \item{mPg}{Class proportions conditional on covariates}
  \item{vPg_avg}{Conditional class proportions averaged over units}
  \item{mPhi}{Conditional response probabilities}
  \item{gammas}{Intercept parameters in logistic models for response probabilities}
  \item{betas}{Intercept and slope parameters in logistic models for conditional class membership}
  \item{parvec}{Vector of all model parameters}
  \item{mV2}{Inverse of the information matrix from the regression stage}
  \item{mQ}{Matrix of cross-derivatives for the asymptotic standard error correction in two-step estimation (see Bakk & Kuha, 2018; Di Mari et al., 2023)}
  \item{Varmat_unc}{Uncorrected variance-covariance matrix}
  \item{Varmat_cor}{Corrected variance-covariance matrix}
  \item{SEs_unc}{Uncorrected standard errors}
  \item{SEs_cor}{Corrected standard errors}
  \item{SEs_cor_beta}{Corrected standard errors for beta}
  \item{eps}{Epsilon: difference between the last two elements of the EM log-likelihood sequence}
  \item{mU}{Posterior class membership probabilities}
  \item{mClassErr}{Expected number of classification errors}
  \item{mClassErrProb}{Expected proportion of classification errors}
  \item{AvgClassErrProb}{Average expected proportion of classification errors}
  \item{R2entr}{Entropy R\eqn{^2}}
  \item{AIC}{Akaike Information Criterion}
  \item{BIC}{Bayesian Information Criterion}
  \item{iter}{Number of iterations of EM algorithm}
  \item{LLKSeries}{Full log-likelihood series of EM algorithm}
  \item{spec}{Specification}
Multilevel measurement model estimation returns (if \code{extout = TRUE}, a subset):
  \item{vOmega}{High-level class proportions}
  \item{mPi}{Conditional low-level class proportions}
  \item{mPhi}{Conditional response probabilities}
  \item{vDelta}{Intercept parameters in logistic models for high-level class proportions}
  \item{mGamma}{Intercept parameters in logistic models for conditional low-level class proportions}
  \item{mBeta}{Intercept parameters in logistic models for response probabilities}
  \item{parvec}{Vector of all model parameters}
  \item{Infomat}{Expected information matrix}
  \item{Varmat}{Variance-covariance matrix}
  \item{SEs}{Standard errors}
  \item{eps}{Epsilon: difference between the last two elements of the EM log-likelihood sequence}
  \item{cPMX}{Joint posterior class membership probabilities}
  \item{cLogPMX}{Log of joint posterior class membership probabilities}
  \item{cPX}{Conditional posterior low-level class membership probabilities}
  \item{cLogPX}{Log of conditional posterior low-level class membership probabilities}
  \item{mSumPX}{Posterior high-level class membership probabilities for low-level units after marginalizing over low-level classes}
  \item{mPW}{Posterior high-level class membership probabilities for high-level units}
  \item{mlogPW}{Log of posterior high-level class membership probabilities for high-level units}
  \item{mPW_N}{Posterior high-level class membership probabilities for low-level units}
  \item{mPMsumX}{Posterior low-level class membership probabilities for low-level units after marginalizing over high-level classes}
  \item{R2entr_low}{Low-level entropy R\eqn{^2}}
  \item{R2entr_high}{High-level entropy R\eqn{^2}}
  \item{AIC}{Akaike Information Criterion}
  \item{BIClow}{Low-level Bayesian Information Criterion}
  \item{BIChigh}{High-level Bayesian Information Criterion}
  \item{ICL_BIClow}{Low-level BIC-type approximation the integrated complete likelihood}
  \item{ICL_BIChigh}{High-level BIC-type approximation the integrated complete likelihood}
  \item{iter}{Number of iterations of EM algorithm}
  \item{LLKSeries}{Full log-likelihood series of EM algorithm}
  \item{vLLK}{Current log-likelihood for high-level units}
  \item{mScore}{Individual contributions to log-likelihood score}
  \item{spec}{Specification}
Multilevel structural model with low-level covariates estimation returns (if \code{extout = TRUE}, a subset):
  \item{vOmega}{High-level class proportions}
  \item{mPi}{Low-level class proportions conditional on high-level class membership and covariates}
  \item{mPi_avg}{Conditional low-level class proportions averaged over low-level units}
  \item{mPhi}{Conditional response probabilities}
  \item{vDelta}{Intercept parameters in logistic models for high-level class proportions}
  \item{cGamma}{Intercept and slope parameters in logistic models for conditional low-level class membership}
  \item{mBeta}{Intercept parameters in logistic models for response probabilities}
  \item{parvec}{Vector of all model parameters}
  \item{Infomat}{Expected information matrix}
  \item{cGamma_Info}{Expected information matrix for gamma}
  \item{mV2}{Inverse of the information matrix from the regression stage}
  \item{mQ}{Matrix of cross-derivatives for the asymptotic standard error correction in two-step estimation (see Bakk & Kuha, 2018; Di Mari et al., 2023)}
  \item{Varmat_unc}{Uncorrected variance-covariance matrix}
  \item{Varmat_cor}{Corrected variance-covariance matrix}
  \item{SEs_unc}{Uncorrected standard errors}
  \item{SEs_cor}{Corrected standard errors}
  \item{SEs_cor_gamma}{Corrected standard errors for gamma}
  \item{eps}{Epsilon: difference between the last two elements of the EM log-likelihood sequence}
  \item{cPMX}{Joint posterior class membership probabilities}
  \item{cLogPMX}{Log of joint posterior class membership probabilities}
  \item{cPX}{Conditional posterior low-level class membership probabilities}
  \item{cLogPX}{Log of conditional posterior low-level class membership probabilities}
  \item{mSumPX}{Posterior high-level class membership probabilities for low-level units after marginalizing over low-level classes}
  \item{mPW}{Posterior high-level class membership probabilities for high-level units}
  \item{mlogPW}{Log of posterior high-level class membership probabilities for high-level units}
  \item{mPW_N}{Posterior high-level class membership probabilities for low-level units}
  \item{mPMsumX}{Posterior low-level class membership probabilities for low-level units after marginalizing over high-level classes}
  \item{R2entr_low}{Low-level entropy R\eqn{^2}}
  \item{R2entr_high}{High-level entropy R\eqn{^2}}
  \item{AIC}{Akaike Information Criterion}
  \item{BIClow}{Low-level Bayesian Information Criterion}
  \item{BIChigh}{High-level Bayesian Information Criterion}
  \item{ICL_BIClow}{Low-level BIC-type approximation the integrated complete likelihood}
  \item{ICL_BIChigh}{High-level BIC-type approximation the integrated complete likelihood}
  \item{iter}{Number of iterations of EM algorithm}
  \item{LLKSeries}{Full log-likelihood series of EM algorithm}
  \item{vLLK}{Current log-likelihood for high-level units}
  \item{mScore}{Individual contributions to log-likelihood score}
  \item{mGamma_Score}{Individual contributions to log-likelihood score for gamma}
  \item{spec}{Specification}
Multilevel structural model with low- and high-level covariates estimation returns (if \code{extout = TRUE}, a subset):
  \item{mOmega}{High-level class proportions conditional on covariates}
  \item{vOmega_avg}{High-level class proportions averaged over high-level units}
  \item{mPi}{Low-level class proportions conditional on high-level class membership and covariates}
  \item{mPi_avg}{Conditional low-level class proportions averaged over low-level units}
  \item{mPhi}{Conditional response probabilities}
  \item{mDelta}{Intercept and slope parameters in logistic models for conditional high-level class membership}
  \item{cGamma}{Intercept and slope parameters in logistic models for conditional low-level class membership}
  \item{mBeta}{Intercept parameters in logistic models for response probabilities}
  \item{parvec}{Vector of all model parameters}
  \item{Infomat}{Expected information matrix}
  \item{cDelta_Info}{Expected information matrix for delta}
  \item{cGamma_Info}{Expected information matrix for gamma}
  \item{mV2}{Inverse of the information matrix from the regression stage}
  \item{mQ}{Matrix of cross-derivatives for the asymptotic standard error correction in two-step estimation (see Bakk & Kuha, 2018; Di Mari et al., 2023)}
  \item{Varmat_unc}{Uncorrected variance-covariance matrix}
  \item{Varmat_cor}{Corrected variance-covariance matrix}
  \item{SEs_unc}{Uncorrected standard errors}
  \item{SEs_cor}{Corrected standard errors}
  \item{SEs_cor_delta}{Corrected standard errors for delta}
  \item{SEs_cor_gamma}{Corrected standard errors for gamma}
  \item{eps}{Epsilon: difference between the last two elements of the EM log-likelihood sequence}
  \item{cPMX}{Joint posterior class membership probabilities}
  \item{cLogPMX}{Log of joint posterior class membership probabilities}
  \item{cPX}{Conditional posterior low-level class membership probabilities}
  \item{cLogPX}{Log of conditional posterior low-level class membership probabilities}
  \item{mSumPX}{Posterior high-level class membership probabilities for low-level units after marginalizing over low-level classes}
  \item{mPW}{Posterior high-level class membership probabilities for high-level units}
  \item{mlogPW}{Log of posterior high-level class membership probabilities for high-level units}
  \item{mPW_N}{Posterior high-level class membership probabilities for low-level units}
  \item{mPMsumX}{Posterior low-level class membership probabilities for low-level units after marginalizing over high-level classes}
  \item{R2entr_low}{Low-level entropy R\eqn{^2}}
  \item{R2entr_high}{High-level entropy R\eqn{^2}}
  \item{AIC}{Akaike Information Criterion}
  \item{BIClow}{Low-level Bayesian Information Criterion}
  \item{BIChigh}{High-level Bayesian Information Criterion}
  \item{ICL_BIClow}{Low-level BIC-type approximation the integrated complete likelihood}
  \item{ICL_BIChigh}{High-level BIC-type approximation the integrated complete likelihood}
  \item{iter}{Number of iterations of EM algorithm}
  \item{LLKSeries}{Full log-likelihood series of EM algorithm}
  \item{vLLK}{Current log-likelihood for high-level units}
  \item{mScore}{Individual contributions to log-likelihood score}
  \item{mDelta_Score}{Individual contributions to log-likelihood score for delta}
  \item{mGamma_Score}{Individual contributions to log-likelihood score for gamma}
  \item{spec}{Specification}
}
\references{
Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between latent classes and external variables. \emph{Psychometrika, 83, 871-892}.

Bakk, Z., Di Mari, R., Oser, J., & Kuha, J. (2022). Two-stage multilevel latent class analysis with covariates in the presence of direct effects. \emph{Structural Equation Modeling: A Multidisciplinary Journal, 29(2), 267-277}.

Di Mari, Bakk, Z., R., Oser, J., & Kuha, J. (2023). A two-step estimator for multilevel latent class analysis with covariates. Under review. Available from \url{https://arxiv.org/abs/2303.06091}.

Lukociene, O., Varriale, R., & Vermunt, J. K. (2010). 6. The simultaneous decision (s) about the number of lower-and higher-level classes in multilevel latent class analysis. Sociological Methodology, 40(1), 247-283.
}
\examples{
\donttest{
# Use the artificial data set
data = dataTOY

# Define vector with names of columns with items
Y = colnames(data)[1+1:10]

# Define name of column with high-level id
id_high = "id_high"

# Define vector with names of columns with low-level covariates
Z = c("Z_low")

# Define vector with names of columns with high-level covariates
Zh = c("Z_high")

# Single-level 3-class LC model with covariates
out = multiLCA(data, Y, 3, Z = Z, verbose = FALSE)
out

# Multilevel LC model
out = multiLCA(data, Y, 3, id_high, 2, verbose = FALSE)
out

# Multilevel LC model low-level covariates
out = multiLCA(data, Y, 3, id_high, 2, Z, verbose = FALSE)
out

# Multilevel LC model low- and high-level covariates
out = multiLCA(data, Y, 3, id_high, 2, Z, Zh, verbose = FALSE)
out

# Model selection over single-level models with 1-3 classes
out = multiLCA(data, Y, iT = 1:3, verbose = FALSE)
out

# Model selection over multilevel models with 1-3 low-level classes and
# 2 high-level classes
out = multiLCA(data, Y, iT = 1:3, id_high, 2, verbose = FALSE)
out

# Model selection over multilevel models with 3 low-level classes and 
# 1-2 high-level classes
out = multiLCA(data, Y, 3, id_high, 1:2, verbose = FALSE)
out

# Model selection over multilevel models with 1-3 low-level classes and 
# 1-2 high-level classes using the default sequential approach
out = multiLCA(data, Y, 1:3, id_high, 1:2, verbose = FALSE)
out
}
}