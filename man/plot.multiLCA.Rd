\name{plot.multiLCA}
\alias{plot.multiLCA}

\title{
Plots conditional response probabilities
}
\description{
Visualizes conditional response probabilities estimated by the \code{\link{multiLCA}} function. The method works for both single- and multilevel models.

Let \code{out} denote the list object returned by the \code{\link{multiLCA}} function. Executing \code{plot(out)} visualizes the conditional response probabilities given by the \code{mPhi} matrix in \code{out}.
}
\usage{
\method{plot}{multiLCA}(x, horiz = FALSE, clab = NULL, \dots)
}
\arguments{
  \item{x}{
    The object returned by the \code{\link{multiLCA}} function
}
  \item{horiz}{
    Whether item labels should be oriented horizontally (\code{TRUE}) or vertically (\code{FALSE}). Default \code{FALSE}
}
  \item{clab}{
    A character vector with user-specified class labels, if available, in the order "Class 1", "Class 2", ... under the default settings, i.e. top-to-bottom. Default \code{NULL}
}
  \item{\dots}{Additional plotting arguments}
}
\value{
  No return value
}
\examples{
\donttest{
# Use IEA data
data = dataIEA

# Define vector with names of columns with items
Y = colnames(data)[4+1:12]

# Define number of (low-level) classes
iT = 3

# Estimate single-level measurement model
out = multiLCA(data = data, Y = Y, iT = iT)
out

# Plot conditional response probabilities with default settings
plot(out)

# Plot with vertical item labels and custom class labels
plot(out, horiz = FALSE, clab = c("Maximal", "Engaged", "Subject"))
}
}