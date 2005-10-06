"baymvb" <-
  function(formula, data = parent.frame(), nvars=8,burnin = 1000, mcmc = 2000,
           thin=1, seed = NA, beta.start = NA, b0 = 0, B0 = 0, R0 = 0.0, G0 = 1, 
           R.start=NA, sd=0.6, distr=c("mvprobit", "mvt"),...)
           {

    check.bayes.parm(burnin, mcmc, thin)
    
    # choice of distribution function
     distr <- match.arg(distr)
     distr <- ifelse(distr=="mvt",1,0)
 

    ## seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## form response and model matrices
    holder <- parse.formula(formula, data)
    Y <- holder[[1]]
    X <- holder[[2]]
    
      
    xnames <- holder[[3]]    
    P  <- ncol(X)  # number of covariates
    N  <- nrow(Y)/nvars  # number of covariates
    
   ## y \in {0, 1,-999} error checking
    if (sum(Y!=0 & Y!=1 & Y!=-999) > 0) {
       cat("Elements of Y equal to something other than 0 or 1 or -999.\n")
       stop("Check data and call baymvb() again. \n") 
    }
    
   ## starting values and priors
    Y2 <- Y + (999 + rbinom(1, 1, 0.5)) * (Y == -999)
    glm.beta <- glm(Y2 ~ X - 1, family = binomial(logit))$coeff
    if (is.na(beta.start)) {
        beta.start <- glm.beta
    }
    if (is.null(dim(beta.start))) {
        beta.start <- beta.start * matrix(1, P, 1)
    }
    if (!is.na(R.start) & is.null(dim(R.start))) {
        R.start <- R.start  + (1 - R.start) *  diag(nvars)
    }
    if (is.null(dim(R.start))) {
        R.start <-  diag(nvars)
    }

    
    if(length(beta.start) != P) {
          cat(paste("Error: Starting value for  beta not
          	comformable [ length(beta) =",length(beta)," != ",P,"].",sep=""),"\n")
          stop("Please respecify and call baymvb() again.\n", call.=FALSE)
    }
    
    if((dim(R.start)[1] != nvars) || (dim(R.start)[2] != nvars)) {
      cat(paste("Error: Starting value for  R not comformable [",nvars," times ",nvars,"].",sep=""),"\n")
      stop("Please respecify and call baymvb() again.\n", call.=FALSE)
    }
    
    mvn.prior <- form.mvn.prior(b0, B0, P )
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]
    
    R.prior <- form.R.prior(R0, G0,nvars)
    R0 <- R.prior[[1]]
    G0 <- R.prior[[2]]
   

    
 ## define holder for posterior density sample
 sample <- matrix(data=0, mcmc/thin, dim(X)[2] + nvars*(nvars-1)/2)
 
 posterior <- .C("baymvb",
            sampledata  = as.double(sample), 
            samplerow   = as.integer(nrow(sample)),
            samplecol   = as.integer(ncol(sample)), 
            Y           = as.double((Y)),
            Yrow        = as.integer(nrow(Y)), 
            Ycol        = as.integer(ncol(Y)), 
            X           = as.double(X),
            Xrow        = as.integer(nrow(X)), 
            Xcol        = as.integer(ncol(X)), 
            burnin      = as.integer(burnin),   
            mcmc        = as.integer(mcmc), 
            thin        = as.integer(thin), 
            lecuyer     = as.integer(lecuyer), 
            seedarray   = as.integer(seed.array), 
            lstream     = as.integer(lecuyer.stream), 
            betastart   = as.double(beta.start), 
            betaow      = as.integer(nrow(beta.start)), 
            betacol     = as.integer(ncol(beta.start)),
            Rstart      = as.double(R.start), 
            Rrow        = as.integer(nrow(R.start)), 
            Rcol        = as.integer(ncol(R.start)), 
            b0          = as.double(b0), 
            b0row       = as.integer(nrow(b0)), 
            b0col       = as.integer(ncol(b0)), 
            B0          = as.double(B0), 
            B0row       = as.integer(nrow(B0)), 
            B0col       = as.integer(ncol(B0)), 
            R0          = as.double(R0), 
            R0row       = as.integer(nrow(R0)), 
            R0col       = as.integer(ncol(R0)), 
            G0          = as.double(G0), 
            G0row       = as.integer(nrow(G0)), 
            G0col       = as.integer(ncol(G0)), 
            N           = as.integer(N), 
            nvars       = as.integer(nvars), 
            P           = as.integer(P),  
            distr       = as.integer(distr), 
            SD          = as.double(sd),
            PACKAGE="baymvb")
    
   ## put together matrix and build MCMC object to return
   rvars  <- numeric(2)
   rvars[1] <- length(xnames) + 1
   Rnames <- NULL
   for(j in 1:(nvars-1))
    for(l in (j+1):nvars)
     Rnames <- c(Rnames,paste("rho_",j,l,sep=""))
     xnames  <- c(xnames, Rnames)
  
   rvars[2] <- length(xnames)

   
    outmat <- matrix(posterior$sampledata,
                     posterior$samplerow,
                     posterior$samplecol,
                     byrow=TRUE)
   
    colnames(outmat) <-  xnames
    
    distrname <- ifelse(distr,"Multivariate t- link","Multivariate probit")
    
    output <-  list(samples = outmat, distrname = distrname,
     coeff =colMeans(outmat, na.rm = TRUE),call=formula,burnin=burnin,
     mcmc=mcmc,thin=thin,rvars =rvars)
     
    class(output) <- c("baymvb")
    return(output) 
  }


print.baymvb <- function(x, digits = max(3, getOption("digits") - 3),...)
{
   cat("\n\nBayesian Analysis of Multivariate Binary data\n")
   cat("\nModel:", x$distrname, "\n")
   print(x$call)
   cat("\nCoefficients:\n")
   print.default(format(x$coeff, digits = digits), print.gap = 2,
            quote = FALSE)
 
   invisible(x)
}


summary.baymvb <- function(object, digits = max(3, getOption("digits") - 3),...)
{
    
  
    statnames <- c("Mean", "SD", "2.5%", "97.5%")
    statsumma <- matrix(nrow = ncol(object$samples), ncol = length(statnames), 
                     dimnames = list(colnames(object$samples), statnames))
                     
    xmean <- colMeans(object$samples, na.rm = TRUE)
    xvar <- diag(var(object$samples, na.rm = TRUE))
    varquant <- apply(object$samples,2, function(x) quantile(x, c(0.025,0.975)))
   
    statsumma[, 1] <- xmean
    statsumma[, 2] <- sqrt(xvar)
    statsumma[, 3:4] <- t(varquant)

    ans <- list()
    ans$samples<- object$samples
    ans$coeff <- statsumma
    ans$distrname <- object$distrname
    ans$call <- object$call
   
     class(ans) <- c("baymvb")
    ans
}


 .First.lib <- function(lib, pkg)
{
   cat("##\n## Bayesian Multivariate Binary Model (baymvb)\n")
   cat("## Copyright (C) 2005, S . M. Mwalili\n")
   cat("##\n## Biostatistical Centre\n")
   cat("## Katholieke Universiteit Leuven\n##\n")
  
    library.dynam("baymvb", pkg, lib)
}
