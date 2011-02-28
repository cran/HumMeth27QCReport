\name{ImportData}
\alias{ImportData}
\title{
ImportData
}
\description{
The function imports the data into a list 
}

\usage{
ImportData(Dir)
}

\arguments{
  \item{Dir}{A character string containing the location of the directory in which the input file are. It will be also the place were the outputs will be}
}
\value{
  A list containing three or four data.frame, one per each file, plus the number of samples and the name of the control file. Moreover it creates a pdf file with the sample names. 
}
\author{ F.M. Mancuso <francesco.mancuso@crg.es> }
\examples{
\dontrun{
	Dir <- system.file("extdata/",package="HumMeth27QCReport")
	ImportDataR <- ImportData(Dir)
}
}
\keyword{ hplot }
\keyword{ manip }
