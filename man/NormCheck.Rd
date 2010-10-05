\name{NormCheck}
\alias{NormCheck}
\title{
Normalization
}
\description{
Normalize the Beta Values and plot a PCA and a hierarchical Clustering of the samples using the noralized data
}

\usage{
NormCheck(Dir)
}

\arguments{
  \item{Dir}{A character string containing the location of the directory 
    in which the input file are
  }
}
\details{
Creates:
\itemize{
\item  a PCA of normalized Beta values
\item  a Cluster of normalized Beta values
\item  an Excel file with the normalized Betas values.
}}
\author{ F.M. Mancuso <francesco.mancuso@crg.es> }

\keyword{ hplot }
\keyword{ manip }
