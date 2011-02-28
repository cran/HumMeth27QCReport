\name{QCCheck}
\alias{QCCheck}
\title{Quality Control of samples}
\description{
Creates all the plots relatve to the quality of the samples.
}
\usage{
QCCheck(ImportDataR, pval)
}
\arguments{
  \item{ImportDataR}{the object that comes from ImportData function}
  \item{pval}{A p-value threshold number to define which sample keep for the normalization and the following analyses }
}
\details{
 Creates 3 different plots to asses the quality of the samples:	
\itemize{ 
\item an Intensity Graph plot foreach sample recalling the "plotSampleIntensities" function of methylumi package.
\item a histogram with the percentage of non dectected CPG (that is the CPGs that have a detection p-value bigger than 0.05 or 0.01.
\item a histogram with the average p-value for each sample.
}}
\value{
 list composed by three data.frame (1. summary of the analysis; 2. list of CPGs that have a detection p-value bigger than 0.05 or 3. than 0.01; a pdf file called 
 QualityCheck.pdf with all the plots.
 }
\author{ F.M. Mancuso <francesco.mancuso@crg.es> }
\examples{
\dontrun{
	Dir <- system.file("extdata/",package="HumMeth27QCReport")
    ImportDataR <- ImportData(Dir)
	QCresults <- QCCheck(ImportDataR,pval=0.05)
}
}
\keyword{ hplot }
\keyword{ manip }
