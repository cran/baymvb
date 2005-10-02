\name{densPlot}
\alias{densPlot}

\title{Provides density plot for individual parameters}
\description{
   `densPlot' offers a flexibility of
    plotting density of each parameter individually
}
\usage{
	densPlot(obj,parm=1, show.obs = TRUE, bwf, main = "", ylim, ...)

}

\arguments{
	
	\item{obj}{posterior samppels from 'baymvb'}
	
	\item{parm}{indexes the parameter to plot}
	
	\item{show.obs}{Show observations along the x-axis}

     	\item{bwf}{Function for calculating the bandwidth.  If omitted, the
          bandwidth is calculate by 1.06 times the minimum of the
          standard deviation and the interquartile range divided by
          1.34 times the sample size to the negative one fifth power}

    	\item{main}{Title. See 'par()'}

    	\item{ylim}{Limits on y axis.  See 'par()'}

     	\item{...}{Further graphical parameters}
     
}

\value{
Density plot
}

\seealso{
\code{\link{tracePlot}},\code{\link{baymvb}}.
}  

\keyword{models}

 \author{S. M. Mwalili \email{samuel\_mwalili@yahoo.com}}
