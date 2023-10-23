\name{multilevLCA-package}
\alias{multilevLCA-package}
\docType{package}
\title{
\packageTitle{multilevLCA}
}
\description{
\packageDescription{multilevLCA}
}
\details{
For estimating latent class models, see \code{\link{multiLCA}}.

For plotting latent class models, see \code{\link{plot.multiLCA}}
}
\author{
Roberto Di Mari and Johan Lyrvall.

Maintainer: \packageMaintainer{multilevLCA}
}
\references{
Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between latent classes and external variables. \emph{Psychometrika, 83, 871-892}.

Bakk, Z., Di Mari, R., Oser, J., & Kuha, J. (2022). Two-stage multilevel latent class analysis with covariates in the presence of direct effects. \emph{Structural Equation Modeling: A Multidisciplinary Journal, 29(2), 267-277}.

Di Mari, Bakk, Z., R., Oser, J., & Kuha, J. (2023). A two-step estimator for multilevel latent class analysis with covariates. Psychometrika.

Lukociene, O., Varriale, R., & Vermunt, J. K. (2010). The simultaneous decision(s) about the number of lower-and higher-level classes in multilevel latent class analysis. Sociological Methodology, 40(1), 247-283.
}

\keyword{ package }
\examples{
\donttest{
data = dataIEA
Y = colnames(dataIEA)[4+1:12]

out = multiLCA(data = data, Y = Y, iT = 2)
out
plot(out, horiz = FALSE)
}
}