\name{baymvb}
\alias{baymvb}
\title{Bayesian analysis of multivariate binary data}
\description{
   `baymvb' is used to fit a multivariate binary data
    regression model in a Bayesian way. In particular, 
    either multivariate probit or multivariate t-link model
    is fitted.
}
\usage{
baymvb(formula, data = parent.frame(), nvars=8,burnin = 1000, mcmc = 2000,
       thin=1, seed = NA, beta.start = NA, b0 = 0, B0 = 0, R0 = 0.0, G0 = 1, 
       R.start=NA, refresh=100, sd=0.3, distr=c("mvprobit", "mvt"),...)

}

\arguments{
    \item{formula}{Model formula.}

    \item{data}{Data frame. NOTE the response is (N*nvar) by 1 vector of 
    response variables,  stacked across all subjects, i.e. data is ordered 
    by subjects and  variables: \eqn{Y_{11},Y_{12},...,Y_{1,nvars},...,Y_{N,nvars}.}}
    
    \item{nvars}{Number of variables.}
    
    \item{distr}{Choice of the link, either ``mvprobit" or ``mvt".}

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of Gibbs iterations for the sampler.}

    \item{thin}{The thinning interval used in the simulation.  The number of
    Gibbs iterations must be divisible by this value.}

    \item{refresh}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If \code{refresh} is greater
    than 0 the iteration number and the MCMC iterations and 
    acceptace rate are printed to the screen every \code{refresh}-th iteration.}
  
    \item{seed}{The seed for the random number generator.  If NA, the Mersenne
      Twister generator is used with default seed 12345; if an integer is 
      passed it is used to seed the Mersenne twister.  The user can also
      pass a list of length two to use the L'Ecuyer random number generator,
      which is suitable for parallel computation.  The first element of the
      list is the L'Ecuyer seed, which is a vector of length six or NA (if NA 
      a default seed of \code{rep(12345,6)} is used).  The second element of 
      list is a positive substream number.}
  
    \item{beta.start}{The starting value for the \eqn{\beta}{beta} vector.    
    This can either be a scalar or a column vector with dimension equal 
    to the number of betas.  If this takes a scalar value, then that value 
    will serve as the starting value for all of the betas. The default 
    value of NA will use the maximum likelihood estimate of \eqn{\beta}{beta} 
    as the starting value.}

    \item{b0}{The prior mean of \eqn{\beta}{beta}.  This can either be a 
    scalar or a column vector with dimension equal to the number of betas. 
    If this takes a scalar value, then that value will serve as the prior mean 
    for all of the betas.}
  
  \item{B0}{The prior precision of \eqn{\beta}{beta}.  This can either
    be a scalar  or a square matrix with dimensions equal to the number of betas.  If this
    takes a scalar value, then that value times an identity matrix serves
    as the prior precision of \eqn{\beta}{beta}. Default value of 0 is 
    equivalent to an improper uniform prior on \eqn{\beta}{beta}.}

  \item{R0}{The prior of the correlation matrix \eqn{R}{R}.  This can either
    be a scalar  or a square matrix with dimensions equal to the number of nvars. 
    Default value of 0 is equivalent to an improper uniform prior on \eqn{R}{R}.}

  \item{G0}{The prior precision of \eqn{R}{R}. This can either
    be a scalar  or a square matrix with dimensions equal to the number of nvars*(nvars-1)/2.  
    If this takes a scalar value, then that value times an identity matrix serves
    as the prior precision of \eqn{R}{R}. Default value of 1 is 
    equivalent to an improper uniform prior on \eqn{R}{R}.}
    
    \item{R.start}{The starting value for the \eqn{R}{R} matrix.    
    This can either be a scalar or a matrix with dimension equal 
    to nvars by nvars.  If this takes a scalar value, then that value 
    will serve as the starting value for all of the corelation. The default 
    value of NA will use the data augmented estimate of \eqn{R}{R} 
    as the starting value.}
 
   \item{sd}{ The hit-and-run MH tuning parmaeter in sampling R}
    
   \item{...}{further arguments to be passed}       
}

\value{
   An mcmc object that contains the posterior density sample.  This 
   object can be summarized by functions provided by the coda package.
}

\references{
  Mwalili, S. M. and Lesaffre, E. 2005. ``Bayesian analysis of 
  misclassified multivariate binary data.'' \emph{In Press} xx, xx--xx.
  
  Chib, S. and Greenberg, E. 1998. ``Bayesian analysis of 
  multivariate probit models.'' \emph{Biometrika} 85, 47--361.

  Chen, M.-H. and Dipak, K. D. 1998. ``Bayesian modelling of 
  correlated binary responses via scale mixture of multivariate
  normal link functions.'' \emph{Sankhy{\={a}}} 60, 322-343.
      

  {O'B}rien, S. M. and Dunson, D. B. 2004. ``Bayesian Multivariate 
  Logistic Regression.'' \emph{Biometrics} 60, 739--746.
 
   Andrew D. Martin, Kevin M. Quinn, and Daniel Pemstein.  2004.  
   \emph{Scythe Statistical Library 1.0.} \url{http://scythe.wustl.edu}.
 
   Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}.
}
 
 \seealso{
 \code{\link{densPlot}},\code{\link{tracePlot}}.
} 

\examples{
library(baymvb) # load library
data(teeth4)   # simulated 4-teeth data
names(teeth4)

# a: Multivariate probit
out.mvprobit <- baymvb(formula=y~ factor(tooth) + gender + age,data=teeth4, 
			nvars=4,burnin = 1000, mcmc = 1000,distr="mvprobit")
out.mvprobit
summary(out.mvprobit)
# b: Multivariate t-link
out.mvt      <- baymvb(formula=y~ factor(tooth) + gender + age,data=teeth4, 
			nvars=4,burnin = 1000, mcmc = 1000,distr="mvt")
out.mvt
summary(out.mvt)			
fit <- glm(y~ factor(tooth) + gender + age,family=binomial(logit),data=teeth4)

# compare with GLM. NOTE: increase the value of burnin and mcmc for better results 
cbind(mvt=colMeans(out.mvt$samples)[1:6],mvprobit=colMeans(out.mvprobit$samples)[1:6],glm=fit$coeff)

# correlation matrix from a vector of unique elements of R
R.mvprobit <- xpndCor(out.mvprobit)
R.mvprobit

R.mvt <- xpndCor(out.mvt)
R.mvt



# Density plots
par(mfrow=c(2,2))
densPlot(out.mvprobit,5,main="Gender(mvprobit)")
densPlot(out.mvprobit,6,main="Age(mvprobit)")
densPlot(out.mvt,5,main="Gender(mvt)")
densPlot(out.mvt,6,main="Age(mvt)")

# Trace plots
par(mfrow=c(2,2))
tracePlot(out.mvprobit,5,main="Gender(mvprobit)")
tracePlot(out.mvprobit,6,main="Age(mvprobit)")
tracePlot(out.mvt,5,main="Gender(mvt)")
tracePlot(out.mvt,6,main="Age(mvt)")

}

\keyword{models}

\author{S. M. Mwalili <samuelmwalili@yahoo.com>}
