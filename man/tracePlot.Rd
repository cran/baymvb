\name{tracePlot}
\alias{tracePlot}

\title{Provides trace for an individual parameter}
\description{
   `ptraceplot' offers a flexibility of
    plotting trace for each parameter alone
}
\usage{
 	tracePlot(obj, parm=1,smooth = TRUE, col, type, ylab, ...)
}

\arguments{
	
	\item{obj}{posterior samppels from 'baymvb'}
	
	\item{parm}{indexes the parameter to plot}
	
	\item{smooth}{draw smooth line through trace plot}

     	\item{col}{graphical parameter (see 'par')}

    	\item{ylab}{graphical parameter (see 'par')}

    	\item{type}{graphical parameter (see 'par')}

     	\item{...}{Further graphical parameters}
     
}

\value{
Trace plot
}
 
\seealso{
\code{\link{densPlot}},\code{\link{baymvb}}.
} 

\keyword{models}

\author{S. M. Mwalili <samuelmwalili@yahoo.com>}
