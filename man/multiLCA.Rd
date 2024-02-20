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
multiLCA(
data,
Y,
iT,
id_high = NULL,
iM = NULL,
Z = NULL,
Zh = NULL,
extout = FALSE,
dataout = TRUE,
kmea = TRUE,
sequential = TRUE,
numFreeCores = 2,
maxIter = 1e3,
tol = 1e-8,
reord = TRUE,
fixedpars = 1,
NRmaxit = 100,
NRtol = 1e-6,
verbose = TRUE
)
}

\arguments{
  \item{data}{
     Input matrix or dataframe.
}
  \item{Y}{
     Names of \code{data} columns with indicators.
}
  \item{iT}{
     Number of lower-level latent classes.
}
  \item{id_high}{
     Name of \code{data} column with higher-level id. Default: \code{NULL}.
}
  \item{iM}{
     Number of higher-level latent classes. Default: \code{NULL}.
}
  \item{Z}{
     Names of \code{data} columns with lower-level covariates (non-numeric covariates are treated as nominal). Default: \code{NULL}.
}
  \item{Zh}{
     Names of \code{data} columns with higher-level covariates (non-numeric covariates are treated as nominal). Default \code{NULL}
}
  \item{extout}{
     Whether to output extensive model and estimation information. Default: \code{FALSE}.
}
  \item{dataout}{
     Whether to match class predictions to the observed data. Default: \code{TRUE}.
}
  \item{kmea}{
     Whether to compute starting values for single-level model using \eqn{K}-means (\code{TRUE}), which is recommended for algorithmic stability, or \eqn{K}-modes (\code{FALSE}). Default: \code{TRUE}.
}
  \item{sequential}{
     Whether to perform sequential model selection (\code{TRUE}) or parallelized model selection (\code{FALSE}). Default: \code{TRUE}.
}
  \item{numFreeCores}{
     If performing parallelized model selection, the number of CPU cores to keep free. Default: \code{2}.
}
  \item{maxIter}{
     Maximum number of iterations for EM algorithm. Default: \code{1e3}.
}
  \item{tol}{
     Tolerance for EM algorithm. Default: \code{1e-8}.
}
  \item{reord}{
     Whether to (re)order classes in decreasing order according to probability of scoring yes on all items. Default: \code{TRUE}.
}
  \item{fixedpars}{
     One-step estimator (\code{0}), two-step estimator (\code{1}) or two-stage estimator (\code{2}). Default: \code{1}.
}
  \item{NRmaxit}{
     Maximum number of iterations for Newton-Raphson algorithm. Default: \code{100}.
}
  \item{NRtol}{
     Tolerance for Newton-Raphson algorithm. Default: \code{1e-6}.
}
  \item{verbose}{
     Whether to print estimation progress. Default: \code{TRUE}.
}

}
\details{
The indicator columns may be coded as as consecutive sequence of integers from 0, or as characters.

To directly estimate a latent class model, \code{iT} and (optionally) \code{iM} should be specified as a single positive integer. To perform model selection over range of consecutive positive integers as the number of latent classes, \code{iT} and/or \code{iM} may be specified in the form \code{iT_min:iT_max} and/or \code{iM_min:iM_max}. It is possible to specify \code{iT = iT_min:iT_max} with either \code{iM = NULL} or \code{iM} equal to a single positive integer, \code{iM = iM_min:iM_max} with \code{iT} equal to a single positive integer, or \code{iT = iT_min:iT_max} with \code{iM = iM_min:iM_max}. All model selection procedures return the output of the optimal model based on the BIC.

In the case where both \code{iT} and \code{iM} are defined as a range of consecutive positive integers, model selection can be performed using the sequential three-stage approach (Lukociene et al., 2010) or a simultaneous approach. The sequential approach involves (first step) estimating \code{iT_min:iT_max} single-level models and identifying the optimal alternative \code{iT_opt1} based on the BIC, (second step) estimating \code{iM_min:iM_max|iT = iT_opt1} multilevel models and identifying the optimal alternative \code{iM_opt2} based on the higher-level BIC, and (third step) estimating \code{iT_min:iT_max|iM = iM_opt2} multilevel models and identifying the optimal alternative \code{iT_opt3} based on the lower-level BIC. The simultaneous approach involves devoting multiple CPU cores on the local machine to estimate all combinations in \code{iT = iT_min:iT_max, iM = iM_min:iM_max} and identifying the optimal alternative based on the lower-level BIC.
}
\value{
Single-level model estimation returns (if \code{extout = TRUE}, a subset):
  \item{vPi}{Class proportions}
  \item{mPhi}{Response probabilities given the latent classes}
  \item{mU}{Matrix of posterior class assignment (proportional assignment)}
  \item{mU_modal}{Matrix of posterior class assignment (modal assignment)}
  \item{vU_modal}{Vector of posterior class assignment (modal assignment)}
  \item{mClassErr}{Expected number of classification errors}
  \item{mClassErrProb}{Expected proportion of classification errors}
  \item{AvgClassErrProb}{Average of \code{mClassErrProb}}
  \item{R2entr}{Entropy-based R\eqn{^2}}
  \item{BIC}{Bayesian Information Criterion (BIC)}
  \item{AIC}{Akaike Information Criterion (AIC)}
  \item{vGamma}{Intercepts in logistic parametrization for class proportions}
  \item{mBeta}{Intercepts in logistic parametrization for response probabilities}
  \item{parvec}{Vector of logistic parameters}
  \item{SEs}{Standard errors}
  \item{Varmat}{Variance-covariance matrix}
  \item{iter}{Number of iterations for EM algorithm}
  \item{eps}{Difference between last two elements of log-likelihood sequence for EM algorithm}
  \item{LLKSeries}{Full log-likelihood series for EM algorithm}
  \item{mScore}{Contributions to log-likelihood score}
  \item{spec}{Model specification}
Single-level model estimation with covariates returns (if \code{extout = TRUE}, a subset):
  \item{mPi}{Class proportions given the covariates}
  \item{vPi_avg}{Sample average of \code{mPi}}
  \item{mPhi}{Response probabilities given the latent classes}
  \item{mU}{Matrix of posterior class assignment (proportional assignment)}
  \item{mClassErr}{Expected number of classification errors}
  \item{mClassErrProb}{Expected proportion of classification errors}
  \item{AvgClassErrProb}{Average of \code{mClassErrProb}}
  \item{R2entr}{Entropy-based R\eqn{^2}}
  \item{BIC}{Bayesian Information Criterion (BIC)}
  \item{AIC}{Akaike Information Criterion (AIC)}
  \item{cGamma}{Intercept and slope parameters in logistic models for conditional class membership}
  \item{mBeta}{Intercepts in logistic parametrization for response probabilities}
  \item{parvec}{Vector of logistic parameters}
  \item{SEs_unc}{Uncorrected standard errors}
  \item{SEs_cor}{Corrected standard errors}
  \item{SEs_cor_gamma}{Corrected standard errors only for the gammas}
  \item{mQ}{Cross-derivatives for asymptotic standard error correction in two-step estimation (see Bakk & Kuha, 2018; Di Mari et al., 2023)}
  \item{Varmat_unc}{Uncorrected variance-covariance matrix}
  \item{Varmat_cor}{Corrected variance-covariance matrix}
  \item{mV2}{Inverse of information matrix for structural model}
  \item{iter}{Number of iterations for EM algorithm}
  \item{eps}{Difference between last two elements of log-likelihood sequence for EM algorithm}
  \item{LLKSeries}{Full log-likelihood series for EM algorithm}
  \item{spec}{Model specification}
Multilevel model estimation returns (if \code{extout = TRUE}, a subset):
  \item{vOmega}{Higher-level class proportions}
  \item{mPi}{Lower-level class proportions given the higher-level latent classes}
  \item{mPhi}{Response probabilities given the lower-level latent classes}
  \item{cPMX}{Posterior joint class assignment (proportional assignment)}
  \item{cLogPMX}{Log of \code{cPMX}}
  \item{cPX}{Posterior lower-level class assignment given high-level class membership (proportional assignment)}
  \item{cLogPX}{Log of \code{cPX}}
  \item{mSumPX}{Posterior higher-level class assignment for lower-level units after marginalization over the lower-level classes (proportional assignment)}
  \item{mPW}{Posterior higher-level class assignment for higher-level units (proportional assignment)}
  \item{mlogPW}{Log of \code{mPW}}
  \item{mPW_N}{Posterior higher-level class assignment for lower-level units (proportional assignment)}
  \item{mPMsumX}{Posterior lower-level class assignment for lower-level units after marginalization over the higher-level classes (proportional assignment)}
  \item{R2entr_low}{Lower-level entropy-based R\eqn{^2}}
  \item{R2entr_high}{Higher-level entropy-based R\eqn{^2}}
  \item{BIClow}{Lower-level Bayesian Information Criterion (BIC)}
  \item{BIChigh}{Higher-level Bayesian Information Criterion (BIC)}
  \item{ICL_BIClow}{Lower-level BIC-type approximation the integrated complete likelihood}
  \item{ICL_BIChigh}{Higher-level BIC-type approximation the integrated complete likelihood}
  \item{AIC}{Akaike Information Criterion (AIC)}
  \item{vAlpha}{Intercepts in logistic parametrization for higher-level class proportions}
  \item{mGamma}{Intercepts in logistic parametrization for conditional lower-level class proportions}
  \item{mBeta}{Intercepts in logistic parametrization for response probabilities}
  \item{parvec}{Vector of logistic parameters}
  \item{SEs}{Standard errors}
  \item{Varmat}{Variance-covariance matrix}
  \item{Infomat}{Expected information matrix}
  \item{iter}{Number of iterations for EM algorithm}
  \item{eps}{Difference between last two elements of log-likelihood sequence for EM algorithm}
  \item{LLKSeries}{Full log-likelihood series for EM algorithm}
  \item{vLLK}{Current log-likelihood for higher-level units}
  \item{mScore}{Contributions to log-likelihood score}
  \item{spec}{Model specification}
Multilevel model estimation with lower-level covariates returns (if \code{extout = TRUE}, a subset):
  \item{vOmega}{Higher-level class proportions}
  \item{mPi}{Lower-level class proportions given the higher-level latent classes and the covariates}
  \item{mPi_avg}{Sample average of \code{mPi}}
  \item{mPhi}{Response probabilities given the lower-level latent classes}
  \item{cPMX}{Posterior joint class assignment (proportional assignment)}
  \item{cLogPMX}{Log of \code{cPMX}}
  \item{cPX}{Posterior lower-level class assignment given high-level class membership (proportional assignment)}
  \item{cLogPX}{Log of \code{cPX}}
  \item{mSumPX}{Posterior higher-level class assignment for lower-level units after marginalization over the lower-level classes (proportional assignment)}
  \item{mPW}{Posterior higher-level class assignment for higher-level units (proportional assignment)}
  \item{mlogPW}{Log of \code{mPW}}
  \item{mPW_N}{Posterior higher-level class assignment for lower-level units (proportional assignment)}
  \item{mPMsumX}{Posterior lower-level class assignment for lower-level units after marginalization over the higher-level classes (proportional assignment)}
  \item{R2entr_low}{Lower-level entropy-based R\eqn{^2}}
  \item{R2entr_high}{Higher-level entropy-based R\eqn{^2}}
  \item{BIClow}{Lower-level Bayesian Information Criterion (BIC)}
  \item{BIChigh}{Higher-level Bayesian Information Criterion (BIC)}
  \item{ICL_BIClow}{Lower-level BIC-type approximation the integrated complete likelihood}
  \item{ICL_BIChigh}{Higher-level BIC-type approximation the integrated complete likelihood}
  \item{AIC}{Akaike Information Criterion (AIC)}
  \item{vAlpha}{Intercepts in logistic parametrization for higher-level class proportions}
  \item{cGamma}{Intercept and slope parameters in logistic models for conditional lower-level class membership}
  \item{mBeta}{Intercepts in logistic parametrization for response probabilities}
  \item{parvec}{Vector of logistic parameters}
  \item{SEs_unc}{Uncorrected standard errors}
  \item{SEs_cor}{Corrected standard errors}
  \item{SEs_cor_gamma}{Corrected standard errors only for the gammas}
  \item{mQ}{Cross-derivatives for asymptotic standard error correction in two-step estimation (see Bakk & Kuha, 2018; Di Mari et al., 2023)}
  \item{Varmat_unc}{Uncorrected variance-covariance matrix}
  \item{Varmat_cor}{Corrected variance-covariance matrix}
  \item{Infomat}{Expected information matrix}
  \item{cGamma_Info}{Expected information matrix only for the gammas}
  \item{mV2}{Inverse of information matrix for structural model}
  \item{iter}{Number of iterations for EM algorithm}
  \item{eps}{Difference between last two elements of log-likelihood sequence for EM algorithm}
  \item{LLKSeries}{Full log-likelihood series for EM algorithm}
  \item{vLLK}{Current log-likelihood for higher-level units}
  \item{mScore}{Contributions to log-likelihood score}
  \item{mGamma_Score}{Contributions to log-likelihood score only for the gammas}
  \item{spec}{Model specification}
Multilevel model estimation with lower- and higher-level covariates returns (if \code{extout = TRUE}, a subset):
  \item{mOmega}{Higher-level class proportions given the covariates}
  \item{vOmega_avg}{Higher-level class proportions averaged over higher-level units}
  \item{mPi}{Lower-level class proportions given the higher-level latent classes and the covariates}
  \item{mPi_avg}{Sample average of \code{mPi}}
  \item{mPhi}{Response probabilities given the lower-level latent classes}
  \item{cPMX}{Posterior joint class assignment (proportional assignment)}
  \item{cLogPMX}{Log of \code{cPMX}}
  \item{cPX}{Posterior lower-level class assignment given high-level class membership (proportional assignment)}
  \item{cLogPX}{Log of \code{cPX}}
  \item{mSumPX}{Posterior higher-level class assignment for lower-level units after marginalization over the lower-level classes (proportional assignment)}
  \item{mPW}{Posterior higher-level class assignment for higher-level units (proportional assignment)}
  \item{mlogPW}{Log of \code{mPW}}
  \item{mPW_N}{Posterior higher-level class assignment for lower-level units (proportional assignment)}
  \item{mPMsumX}{Posterior lower-level class assignment for lower-level units after marginalization over the higher-level classes (proportional assignment)}
  \item{R2entr_low}{Lower-level entropy-based R\eqn{^2}}
  \item{R2entr_high}{Higher-level entropy-based R\eqn{^2}}
  \item{BIClow}{Lower-level Bayesian Information Criterion (BIC)}
  \item{BIChigh}{Higher-level Bayesian Information Criterion (BIC)}
  \item{ICL_BIClow}{Lower-level BIC-type approximation the integrated complete likelihood}
  \item{ICL_BIChigh}{Higher-level BIC-type approximation the integrated complete likelihood}
  \item{AIC}{Akaike Information Criterion (AIC)}
  \item{mAlpha}{Intercept and slope parameters in logistic models for conditional higher-level class membership}
  \item{cGamma}{Intercept and slope parameters in logistic models for conditional lower-level class membership}
  \item{mBeta}{Intercepts in logistic parametrization for response probabilities}
  \item{parvec}{Vector of logistic parameters}
  \item{SEs_unc}{Uncorrected standard errors}
  \item{SEs_cor}{Corrected standard errors}
  \item{SEs_cor_alpha}{Corrected standard errors only for the alphas}
  \item{SEs_cor_gamma}{Corrected standard errors only for the gammas}
  \item{mQ}{Cross-derivatives for asymptotic standard error correction in two-step estimation (see Bakk & Kuha, 2018; Di Mari et al., 2023)}
  \item{Varmat_unc}{Uncorrected variance-covariance matrix}
  \item{Varmat_cor}{Corrected variance-covariance matrix}
  \item{Infomat}{Expected information matrix}
  \item{cAlpha_Info}{Expected information matrix only for the alphas}
  \item{cGamma_Info}{Expected information matrix only for the gammas}
  \item{mV2}{Inverse of information matrix for structural model}
  \item{iter}{Number of iterations for EM algorithm}
  \item{eps}{Difference between last two elements of log-likelihood sequence for EM algorithm}
  \item{LLKSeries}{Full log-likelihood series for EM algorithm}
  \item{vLLK}{Current log-likelihood for higher-level units}
  \item{mScore}{Contributions to log-likelihood score}
  \item{mAlpha_Score}{Contributions to log-likelihood score only for the alphas}
  \item{mGamma_Score}{Contributions to log-likelihood score only for the gammas}
  \item{spec}{Model specification}
}
\references{
Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between latent classes and external variables. \emph{Psychometrika, 83, 871-892}.

Bakk, Z., Di Mari, R., Oser, J., & Kuha, J. (2022). Two-stage multilevel latent class analysis with covariates in the presence of direct effects. \emph{Structural Equation Modeling: A Multidisciplinary Journal, 29(2), 267-277}.

Di Mari, Bakk, Z., R., Oser, J., & Kuha, J. (2023). A two-step estimator for multilevel latent class analysis with covariates. Psychometrika.

Lukociene, O., Varriale, R., & Vermunt, J. K. (2010). The simultaneous decision(s) about the number of lower-and higher-level classes in multilevel latent class analysis. Sociological Methodology, 40(1), 247-283.
}
\examples{
\donttest{
# Use the artificial data set
data = dataTOY

# Define vector with names of columns with items
Y = colnames(data)[1+1:10]

# Define name of column with higher-level id
id_high = "id_high"

# Define vector with names of columns with lower-level covariates
Z = c("Z_low")

# Define vector with names of columns with higher-level covariates
Zh = c("Z_high")

# Single-level 3-class LC model with covariates
out = multiLCA(data, Y, 3, Z = Z, verbose = FALSE)
out

# Multilevel LC model
out = multiLCA(data, Y, 3, id_high, 2, verbose = FALSE)
out

# Multilevel LC model lower-level covariates
out = multiLCA(data, Y, 3, id_high, 2, Z, verbose = FALSE)
out

# Multilevel LC model lower- and higher-level covariates
out = multiLCA(data, Y, 3, id_high, 2, Z, Zh, verbose = FALSE)
out

# Model selection over single-level models with 1-3 classes
out = multiLCA(data, Y, 1:3, verbose = FALSE)
out

# Model selection over multilevel models with 1-3 lower-level classes and
# 2 higher-level classes
out = multiLCA(data, Y, 1:3, id_high, 2, verbose = FALSE)
out

# Model selection over multilevel models with 3 lower-level classes and 
# 1-2 higher-level classes
out = multiLCA(data, Y, 3, id_high, 1:2, verbose = FALSE)
out

# Model selection over multilevel models with 1-3 lower-level classes and 
# 1-2 higher-level classes using the default sequential approach
out = multiLCA(data, Y, 1:3, id_high, 1:2, verbose = FALSE)
out
}
}