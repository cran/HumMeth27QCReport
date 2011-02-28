\name{QCCheck}
\alias{QCCheck}
\title{Quality Control of samples}
\description{
Creates all the plots relatve to the quality of the samples.
}
\usage{
QCCheck(Dir, pval)
}
\arguments{
  \item{Dir}{A character string containing the location of the directory in which the input file are }
  \item{pval}{A p-value threshold number to define which sample keep for the normalization and the following analyses }
}
\details{
 Creates 3 different plots to asses the quality of the samples:	
\itemize{ 
\item an Intensity Graph plot foreach sample recalling the "plotSampleIntensities" function of methylumi package.
\item a histogram with the percentage of non dectected CPG (that is the CPGs tha have a detection p-value bigger than 0.05 or 0.01.
\item a histogram with the average p-value for each sample.
}}

\author{ F.M. Mancuso <francesco.mancuso@crg.es> }
\examples{
\dontrun{
	Dir <- system.file("extdata/",package="HumMeth27QCReport")
	QCCheck(Dir,pval=0.05)
}
}
\keyword{ hplot }
\keyword{ manip }
