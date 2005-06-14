\name{hidden}
\alias{cor.start}		 
\alias{form.R.prior}		
\alias{check.bayes.parm}	
\alias{calling.function}	
\alias{form.seeds}		
\alias{parse.formula}		
\alias{form.mvn.prior}	
\alias{xpndCor}	

\title{Helper/utility  functions required by 'baymvb'}
\description{
  Provides helper/utility functions required by 'baymvb'
}
\usage{
        cor.start(Y,n, nvars) 
	form.R.prior(R0, G0, nvars)
        check.bayes.parm(burnin, mcmc, thin) 
  	calling.function(parentheses=TRUE)   
  	form.seeds(seed)   
	parse.formula(formula, data, intercept=TRUE, justX=FALSE)
  	form.mvn.prior(b0, B0, K)
  	xpndCor(obj)
}


\arguments{
    \item{formula}{Model formula.}

    \item{data}{Data frame.}
    
    \item{nvars}{Number of variables.}
    
   \item{Y}{Response.}
        
    \item{n}{Length of Y.}
    
    \item{K}{nvars.}

    \item{posterior.object}{posterior object.}
    
    \item{names}{names}
         
          
    \item{intercept}{intercept}
            
    \item{title}{title}
   
    \item{justX}{justX}
    
    \item{title}{title}
        
        
    \item{parentheses}{TRUE or FALSE} 

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of Gibbs iterations for the sampler.}

    \item{thin}{The thinning interval used in the simulation.  The number of
     Gibbs iterations must be divisible by this value.}

    \item{x}{An object that may be coerced to an mcmc object}

    \item{start}{the iteration number of the first observation}

    \item{end}{the iteration number of the last observation}
    
 
    \item{seed}{The seed for the random number generator.  If NA, the Mersenne
      Twister generator is used with default seed 12345; if an integer is 
      passed it is used to seed the Mersenne twister.  The user can also
      pass a list of length two to use the L'Ecuyer random number generator,
      which is suitable for parallel computation.  The first element of the
      list is the L'Ecuyer seed, which is a vector of length six or NA (if NA 
      a default seed of \code{rep(12345,6)} is used).  The second element of 
      list is a positive substream number.}
  

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
    
  \item{obj}{Posterior samples from `baymvb'}  
  
  \item{...}{a list of mcmc objects}
    
}

\keyword{internal}
