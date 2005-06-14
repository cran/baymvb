\name{summary.baymvb}
\alias{summary.baymvb}
\alias{print.baymvb}
\title{Printing summary of the `baymvb' mcmc output}
\description{
    `print.baymvb' and  `summary.baymvb'provide an efficient way of 
     summarizing the mcmc from a Bayesian Multivariate Binary Regression
}
\usage{
 print.baymvb(x, digits = max(3, getOption("digits") - 3),...)

 summary.baymvb(object, digits = max(3, getOption("digits") - 3),...)

}


\arguments{
  \item{x}{an object from the `baymvb' fit}
  \item{object}{an object from the `baymvb' fit}
  \item{digits}{number of significant decimal points.}  
  \item{...}{other options `summary' arguments.}  

}


\seealso{
\code{\link{tracePlot}},\code{\link{densPlot}}.
}

\keyword{models}
\author{S. M. Mwalili <samuelmwalili@yahoo.com>}
