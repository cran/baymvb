"cor.start" <-
  function(Y,n, nvars) {
  
   Z <- Y
   
    for(i in 1:(n*nvars))
    {
        if(Y[i]==0) Z[i] <- -abs(rnorm(1,1,1))
        if(Y[i]==1) Z[i] <-  abs(rnorm(1,-1,1))
        if(Y[i]==0 && Y[i]==1) Z[i] <- rnorm(1)
    }

    zz      <- matrix(Z,ncol=nvars,byrow=TRUE)
    muz     <-  zz - t(colMeans(zz)%x%matrix(1,1,n))
    varZ    <-  t(muz)%*%muz
    stdZ    <-  apply(zz,2,sd)
    R.start <- ((1.0/(n-1))*varZ/(stdZ%x%t(stdZ)))

    return(R.start)
   }
    
# form multivariate Normal prior
"form.R.prior" <-
   function(R0, G0, nvars) {
   
    nvar2 <- nvars*(nvars-1)/2
    
    # prior precision of R
     if(is.null(dim(R0))) {
       R0 <- R0*matrix(1,nvars,nvars) + (1-R0) * diag(nvars)    
     }
     if((dim(R0)[1] != nvars) || (dim(R0)[2] != nvars)) {
       cat("Error: pi(R|R0,G0^-1) prior R0 not conformable.\n")
       stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
     }
     if(is.null(dim(G0))) {
       G0 <- G0*diag(nvar2)    
     }
     if((dim(G0)[1] != nvar2) || (dim(G0)[2] != nvar2)) {
       cat("Error: pi(R|R0,G0^-1) prior G0 not conformable.\n")
       stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
     }
     
     return(list(R0,G0))
   }
   
"check.bayes.parm" <-
  function(burnin, mcmc, thin) {
  
    if(mcmc %% thin != 0) {
      cat("Error: MCMC iterations not evenly divisible by thinning interval.\n")
      stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
    }
    if(mcmc < 0) {
      cat("Error: MCMC iterations negative.\n")
      stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE) 
    }
    if(burnin < 0) {
      cat("Error: Burnin iterations negative.\n")
      stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
    }
    if(thin < 0) {
      cat("Error: Thinning interval negative.\n")
      stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
    }    
    return(0)
  }
  
"calling.function" <-
   function(parentheses=TRUE) {
     calling.function <- strsplit(toString(sys.call(which=-3)),",")[[1]][1]
     if (parentheses){
       calling.function <- paste(calling.function, "()", sep="")
     }
     return(calling.function)
   }
  
  
# parse the passed seeds
# 1] if a scalar is passed, it is used by Mersennse twister
# 2] if a list of length two is passed, a parallel-friendly stream is
#    created using L'Ecuyer
"form.seeds" <-
   function(seed) {
      if(length(seed)==1) {
         if(is.na(seed)) seed <- 12345
         seed <- as.integer(seed)
         if(seed < 0) {
            cat("Error: Mersenne seed negative.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)                       
         }
         seeds <- list(0, rep(seed,6), 0)
      }
      if(length(seed)==2) {
         if(!is.list(seed)) {
            cat("Error: List must be passed to use L'Ecuyer.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)          
         }
         lec.seed <- seed[[1]]
         lec.substream <- as.integer(seed[[2]])
         if(is.na(lec.seed[1])) lec.seed <- rep(12345, 6)
         if(length(lec.seed) != 6) {
            cat("Error: L'Ecuyer seed not of length six.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)          
         }
         if(!all(lec.seed >= 0))  {
             cat("Error: At least one L'Ecuyer seed negative.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)          
         }
         if( max(lec.seed[1:3]) >= 4294967087){
           cat("Error: At least one of first three L'Ecuyer seeds\n")
           cat("  greater than or equal to 4294967087\n")
           stop("Please respecify and call ", calling.function(), " again.",
                call.=FALSE)          
         }
         if( all(lec.seed[1:3]) == 0 ){
           cat("Error: first three L'Ecuyer seeds == 0\n")
           stop("Please respecify and call ", calling.function(), " again.",
                call.=FALSE)          
         }
         if( max(lec.seed[4:6]) >= 4294944443){
           cat("Error: At least one of last three L'Ecuyer seeds\n")
           cat("  greater than or equal to 4294944443\n")
           stop("Please respecify and call ", calling.function(), " again.",
                call.=FALSE)          
         }         
         if( all(lec.seed[4:6]) == 0 ){
           cat("Error: last three L'Ecuyer seeds == 0\n")
           stop("Please respecify and call ", calling.function(), " again.",
                call.=FALSE)          
         }
         if(lec.substream < 1) {
            cat("Error: L'Ecuyer substream number not positive.\n")
            stop("Please respecify and call ", calling.function(), " again.",
                 call.=FALSE)               
         }
         seeds <- list(1, lec.seed, lec.substream) 
      }
      if(length(seed)>2) {
            cat("Error: Seed passed as length greater than two.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)        
      }
      return(seeds)
   }


"parse.formula" <- 
   function(formula, data, intercept=TRUE, justX=FALSE) {

    # extract Y, X, and variable names for model formula and frame
    mt <- terms(formula, data=data)
    if(missing(data)) data <- sys.frame(sys.parent())
    mf <- match.call(expand.dots = FALSE)
    mf$intercept <- mf$justX <- NULL
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
    if (!intercept){
      attributes(mt)$intercept <- 0
    }

    # null model support
    X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
    X <- as.matrix(X)         # X matrix
    xvars <- dimnames(X)[[2]] # X variable names
    xobs  <- dimnames(X)[[1]] # X observation names
    if (justX){
      Y <- NULL
    }
    else {
      Y <- as.matrix(model.response(mf, "numeric")) # Y matrix
    }
    return(list(Y, X, xvars, xobs))
   }

   
# form multivariate Normal prior
"form.mvn.prior" <-
   function(b0, B0, K) {
  
     # prior mean
     if(is.null(dim(b0))) {
       b0 <- b0 * matrix(1,K,1)  
     } 
     if((dim(b0)[1] != K) || (dim(b0)[2] != 1)) {
       cat("Error: N(b0,B0^-1) prior b0 not conformable.\n")
       stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
     }
     
     # prior precision
     if(is.null(dim(B0))) {
       B0 <- B0 * diag(K)    
     }
     if((dim(B0)[1] != K) || (dim(B0)[2] != K)) {
       cat("Error: N(b0,B0^-1) prior B0 not conformable.\n")
       stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
     }
     return(list(b0,B0))
   }



# trace plot for each parameter
"tracePlot" <-
function (obj, parm=1,smooth = TRUE, col = 1:6, type = "l", ylab = "", 
    ...) 
{
        j <- parm
        
        x <- obj$samples
     
        xp <- obj$burnin + obj$thin - 1 + seq(1, obj$mcmc, by = obj$thin)
        yp <- x[, j, drop = TRUE]
        plot(xp, yp, xlab = "Iterations", ylab = ylab, type = type, 
            col = col, ...)
       
        if (smooth) {
           lines(lowess(xp, yp))
        }
}

# density plot for each parameter

"densPlot" <- 
function (obj, parm=1, show.obs = TRUE, bwf, main = "", ylim, ...) 
{
         j <- parm
        
        x <- obj$samples
     
    #    xp <- obj$burnin + obj$thin - 1 + seq(1, obj$mcmc, by = obj$thin)
        y <- x[, j, drop = TRUE]
        
        if (missing(bwf)) 
            bwf <- function(x) {
                x <- x[!is.na(as.vector(x))]
                return(1.06 * min(sd(x), IQR(x)/1.34) * length(x)^-0.2)
            }
        bw <- bwf(y)
        width <- 4 * bw
        if (max(abs(y - floor(y))) == 0 || bw == 0) 
            hist(y, prob = TRUE, main = main, ...)
        else {
            scale <- "open"
            if (max(y) <= 1 && 1 - max(y) < 2 * bw) {
                if (min(y) >= 0 && min(y) < 2 * bw) {
                  scale <- "proportion"
                  y <- c(y, -y, 2 - y)
                }
            }
            else if (min(y) >= 0 && min(y) < 2 * bw) {
                scale <- "positive"
                y <- c(y, -y)
            }
            else scale <- "open"
            dens <- density(y, width = width)
            if (scale == "proportion") {
                dens$y <- 3 * dens$y[dens$x >= 0 & dens$x <= 
                  1]
                dens$x <- dens$x[dens$x >= 0 & dens$x <= 1]
            }
            else if (scale == "positive") {
                dens$y <- 2 * dens$y[dens$x >= 0]
                dens$x <- dens$x[dens$x >= 0]
            }
            if (missing(ylim)) 
                ylim <- c(0, max(dens$y))
            plot(dens, ylab = "", main = main, type = "l", xlab = paste("N =", 
                length(y), "  Bandwidth =", formatC(dens$bw)), 
                ylim = ylim, ...)
            if (show.obs) 
                lines(y, rep(max(dens$y)/100, length(y)), 
                  type = "h")
        }

    return(invisible(obj))
}
 

# correlation matrix from a vector of unique elements of R
"xpndCor" <- function(obj)
{
        unikR <- colMeans(obj$samples)[obj$rvars[1]:obj$rvars[2]]
	x   <- length(unikR)
  	n   <- (1 + sqrt(1+ 8*x))/2
  
  	R   <- matrix(1,n,n)
  	cnt <- 0
  	for(i in 1:(n-1)){
  		for(j in (i+1):n){
   			R[i,j] <- R[j,i] <- unikR[cnt <- cnt + 1]   			
   		}
   	}
    colnames(R) <- paste("rho_", 1:n, sep="")
    rownames(R) <- paste("rho_", 1:n, sep="")
   
   return(R)
}
 